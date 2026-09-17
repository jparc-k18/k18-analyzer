#!/usr/bin/env python3
"""
K18 Geant4 Simulation & Tracking Pipeline
=========================================
初見のユーザーでも1行でシミュレーションからミッシングマス解析まで
一気通貫で実行できるワンストップパイプラインスクリプト。

使用例:
  # 12C(K-, pi-) 12_LambdaC 反応を 1000 イベント実行
  python3 scripts/run_k18_pipeline.py --reaction kpi_c12lambda -n 1000 -o out/kpi_test

  # 対応反応一覧を表示
  python3 scripts/run_k18_pipeline.py --list-reactions
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# デフォルトのリポジトリパス解決
SCRIPT_DIR = Path(__file__).resolve().parent
ANALYZER_DIR = SCRIPT_DIR.parent
GEANT4_DIR = Path(os.environ.get("K18GEANT4_DIR", ANALYZER_DIR.parent / "k18-geant4"))

SUPPORTED_REACTIONS = {
    "kpi_c12lambda": {
        "description": "12C(K-, pi-)12_LambdaC (1.4 GeV/c, B_Lambda=10.76 MeV)",
        "g4_conf": "param/conf/k18_6381_kpi_c12lambda_example.conf",
        "k18_conf": "runmanager/runlist/k18_6381_kpi_c12lambda_k18.conf",
        "s2s_conf": "runmanager/runlist/k18_6381_kpi_c12lambda_s2s.conf",
        "preset": "kpi_c12lambda",
    },
    "e90_cusp": {
        "description": "d(K-, pi-)Sigma-N cusp (1.4 GeV/c, liquid deuterium cylinder)",
        "g4_conf": "param/conf/k18_6381_e90_cusp_example.conf",
        "k18_conf": "runmanager/runlist/k18_6381_e90_cusp_k18.conf",
        "s2s_conf": "runmanager/runlist/k18_6381_e90_cusp_s2s.conf",
        "preset": "e90_cusp",
    },
    "e75_phase1_li7": {
        "description": "7Li(K-, K+)7_XiH (1.8 GeV/c, enriched Li7 slab)",
        "g4_conf": "param/conf/k18_6381_e75_phase1_li7_example.conf",
        "k18_conf": "runmanager/runlist/k18_6381_e75_phase1_li7_k18.conf",
        "s2s_conf": "runmanager/runlist/k18_6381_e75_phase1_li7_s2s.conf",
        "preset": "e75_phase1_li7",
    },
    "kk_c12xibe": {
        "description": "12C(K-, K+)12_XiBe (1.8 GeV/c, thick C slab)",
        "g4_conf": "param/conf/k18_6381_kk_c12xibe_example.conf",
        "k18_conf": "runmanager/runlist/k18_6381_kk_c12xibe_k18.conf",
        "s2s_conf": "runmanager/runlist/k18_6381_kk_c12xibe_s2s.conf",
        "preset": "kk_c12xibe",
    },
    "pik_c12lambda": {
        "description": "12C(pi+, K+)12_LambdaC (1.4 GeV/c, C slab)",
        "g4_conf": "param/conf/k18_6381_pik_c12lambda_example.conf",
        "k18_conf": "runmanager/runlist/k18_6381_pik_c12lambda_k18.conf",
        "s2s_conf": "runmanager/runlist/k18_6381_pik_c12lambda_s2s.conf",
        "preset": "pik_c12lambda",
    },
}


def run_command(cmd: list[str], cwd: Path, log_file: Path | None = None) -> None:
    """コマンドを実行し、エラー時は詳細を出力して中断する"""
    cmd_str = " ".join(str(c) for c in cmd)
    print(f"  [RUN] {cmd_str}")
    start = time.time()
    if log_file:
        with open(log_file, "w") as f:
            res = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    else:
        res = subprocess.run(cmd, cwd=cwd)
    elapsed = time.time() - start
    if res.returncode != 0:
        print(f"\n  [ERROR] コマンドが失敗しました (終了コード: {res.returncode})", file=sys.stderr)
        if log_file and log_file.exists():
            print(f"  ログ末尾 (詳細は {log_file} を参照):", file=sys.stderr)
            lines = log_file.read_text().splitlines()
            for line in lines[-20:]:
                print(f"    {line}", file=sys.stderr)
        sys.exit(res.returncode)
    print(f"  [DONE] 完了 ({elapsed:.1f} 秒)")


def main():
    parser = argparse.ArgumentParser(
        description="K18 Geant4 Simulation & Tracking Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--reaction", "-r",
        choices=list(SUPPORTED_REACTIONS.keys()),
        help="実行する反応名 (一覧は --list-reactions で確認)",
    )
    parser.add_argument(
        "--events", "-n",
        type=int,
        default=1000,
        help="シミュレーションイベント数 (デフォルト: 1000)",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        help="出力先ディレクトリ (デフォルト: out/<reaction>_<timestamp>)",
    )
    parser.add_argument(
        "--g4-dir",
        type=Path,
        default=GEANT4_DIR,
        help=f"k18-geant4 ディレクトリパス (デフォルト: {GEANT4_DIR})",
    )
    parser.add_argument(
        "--skip-g4",
        action="store_true",
        help="Geant4 シミュレーションをスキップし、既存の G4 ROOT ファイルを使用",
    )
    parser.add_argument(
        "--g4-root",
        type=Path,
        help="--skip-g4 時に使用する既存の G4 ROOT ファイルパス",
    )
    parser.add_argument(
        "--list-reactions",
        action="store_true",
        help="利用可能な反応プリセット一覧を表示して終了",
    )

    args = parser.parse_args()

    if args.list_reactions:
        print("\n利用可能な反応プリセット一覧:")
        print("=" * 60)
        for tag, info in SUPPORTED_REACTIONS.items():
            print(f"  {tag:<18} : {info['description']}")
        print("=" * 60)
        print("使用例:")
        print("  python3 scripts/run_k18_pipeline.py --reaction kpi_c12lambda -n 1000\n")
        return 0

    if not args.reaction:
        parser.error("--reaction (または -r) を指定してください。一覧は --list-reactions で確認できます。")

    rx = SUPPORTED_REACTIONS[args.reaction]

    # 出力先ディレクトリの決定
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    out_dir = args.output_dir or ANALYZER_DIR / "out" / f"{args.reaction}_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = out_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    print("\n" + "=" * 60)
    print(f" K18 Simulation & Tracking Pipeline")
    print("=" * 60)
    print(f" 反応:       {args.reaction} ({rx['description']})")
    print(f" イベント数: {args.events}")
    print(f" 出力先:     {out_dir}")
    print("=" * 60 + "\n")

    # バイナリパスの確認
    g4_bin = args.g4_dir / "bin" / "G4S2S"
    k18_bin = ANALYZER_DIR / "bin" / "DstK18TrackingGeant4"
    s2s_bin = ANALYZER_DIR / "bin" / "DstS2sTrackingGeant4"
    mm_bin = ANALYZER_DIR / "bin" / "DstMissingMassGeant4"

    for name, bin_path in [
        ("G4S2S", g4_bin),
        ("DstK18TrackingGeant4", k18_bin),
        ("DstS2sTrackingGeant4", s2s_bin),
        ("DstMissingMassGeant4", mm_bin),
    ]:
        if not bin_path.exists():
            print(f"[ERROR] 必要なバイナリ {name} が見つかりません: {bin_path}", file=sys.stderr)
            print(f"先に各リポジトリでビルドを実行してください。", file=sys.stderr)
            sys.exit(1)

    # 中間ファイル・出力ファイルパス
    g4_root = out_dir / f"{args.reaction}_g4.root"
    k18_root = out_dir / f"{args.reaction}_k18.root"
    s2s_root = out_dir / f"{args.reaction}_s2s.root"
    mm_root = out_dir / f"{args.reaction}_missing_mass.root"
    mm_csv = out_dir / f"{args.reaction}_missing_mass.csv"

    # Step 1: Geant4 シミュレーション
    print(">>> Step 1/4: Geant4 シミュレーションの実行 (G4S2S)")
    if args.skip_g4:
        if not args.g4_root or not args.g4_root.exists():
            print("[ERROR] --skip-g4 を使用する場合は実在する --g4-root を指定してください。", file=sys.stderr)
            sys.exit(1)
        g4_root = args.g4_root
        print(f"  [SKIP] 既存の G4 ファイルを使用: {g4_root}")
    else:
        # 一時マクロ作成
        with tempfile.NamedTemporaryFile(mode="w", suffix=".mac", delete=False) as mac:
            mac.write(f"/run/beamOn {args.events}\n")
            mac_path = Path(mac.name)
        try:
            g4_conf = args.g4_dir / rx["g4_conf"]
            run_command(
                [str(g4_bin), str(g4_conf), str(g4_root), str(mac_path)],
                cwd=args.g4_dir,
                log_file=logs_dir / "step1_g4.log",
            )
        finally:
            mac_path.unlink(missing_ok=True)

    # Step 2: K18 ビームライントラッキング
    print("\n>>> Step 2/4: K18 ビームライントラッキングの実行 (DstK18TrackingGeant4)")
    k18_conf = ANALYZER_DIR / rx["k18_conf"]
    run_command(
        [str(k18_bin), str(k18_conf), str(g4_root), str(k18_root)],
        cwd=ANALYZER_DIR,
        log_file=logs_dir / "step2_k18.log",
    )

    # Step 3: S-2S スペクトロメータトラッキング
    print("\n>>> Step 3/4: S-2S スペクトロメータトラッキングの実行 (DstS2sTrackingGeant4)")
    s2s_conf = ANALYZER_DIR / rx["s2s_conf"]
    run_command(
        [str(s2s_bin), str(s2s_conf), str(g4_root), str(s2s_root)],
        cwd=ANALYZER_DIR,
        log_file=logs_dir / "step3_s2s.log",
    )

    # Step 4: ミッシングマス & 反応点解析
    print("\n>>> Step 4/4: ミッシングマス・反応点解析の実行 (DstMissingMassGeant4)")
    run_command(
        [
            str(mm_bin),
            f"--preset={rx['preset']}",
            f"--g4={g4_root}",
            f"--k18={k18_root}",
            f"--s2s={s2s_root}",
            f"--root={mm_root}",
            f"--csv={mm_csv}",
        ],
        cwd=ANALYZER_DIR,
        log_file=logs_dir / "step4_missing_mass.log",
    )

    print("\n" + "=" * 60)
    print(" パイプライン処理が正常に完了しました！")
    print("=" * 60)
    print("生成された主要成果物:")
    print(f"  Geant4 RAW データ : {g4_root}")
    print(f"  K18 トラック ROOT : {k18_root}")
    print(f"  S-2S トラック ROOT: {s2s_root}")
    print(f"  ミッシングマス ROOT: {mm_root}")
    print(f"  解析サマリ CSV     : {mm_csv}")
    print(f"  実行ログ一覧       : {logs_dir}")
    print("=" * 60 + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
