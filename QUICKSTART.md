# K18 Geant4 & Analyzer クイックスタートガイド

本ガイドは、J-PARC K1.8 ビームラインおよび S-2S スペクトロメータの Geant4 シミュレーション（`k18-geant4`）と解析ツール（`k18-analyzer`）を初めて使う方のための解説書です。

---

## 1. 一番簡単な動かし方 (ワンストップ実行)

付属のパイプラインスクリプトを使用すると、**1行のコマンドで「シミュレーション → K18トラッキング → S-2Sトラッキング → ミッシングマス解析」まで自動で一括実行**できます。

```bash
# 例: 12C(K-, pi-) 12_LambdaC 反応を 1000 イベント実行
python3 scripts/run_k18_pipeline.py --reaction kpi_c12lambda -n 1000 -o out/kpi_test
```

### 実行結果
指定した出力ディレクトリ（上記例では `out/kpi_test/`）に以下のファイルが生成されます：
- `kpi_c12lambda_g4.root`: Geant4 シミュレーション RAW データ
- `kpi_c12lambda_k18.root`: K18 ビームライントラック再構成結果
- `kpi_c12lambda_s2s.root`: S-2S スペクトロメータトラック再構成結果
- `kpi_c12lambda_missing_mass.root`: ミッシングマスヒストグラムおよび `bbflow` ツリー
- `kpi_c12lambda_missing_mass.csv`: 統計量・カットフロー・分解能フィット結果サマリ
- `logs/`: 各ステージの詳細な実行ログ

---

## 2. 対応している反応プリセット一覧

以下の反応がプリセットとして登録されており、パラメータを個別に設定することなくすぐに実行できます。

| 反応名 (`--reaction` / `--preset`) | 反応内容 | ビーム運動量 | 標的 |
| :--- | :--- | :---: | :--- |
| `kpi_c12lambda` | $^{12}\mathrm{C}(K^-, \pi^-)^{12}_\Lambda\mathrm{C}$ | 1.4 GeV/c | 炭素 ($1.8\ \mathrm{g/cm^3}$, 厚さ 5.6 mm) |
| `e90_cusp` | $d(K^-, \pi^-)\Sigma N$ cusp | 1.4 GeV/c | 液体重水素 ($0.169\ \mathrm{g/cm^3}$, 円筒標的) |
| `e75_phase1_li7` | $^{7}\mathrm{Li}(K^-, K^+)^{7}_{\Xi}\mathrm{H}$ | 1.8 GeV/c | 濃縮リチウム7 ($0.534\ \mathrm{g/cm^3}$, 厚さ 16.9 mm) |
| `kk_c12xibe` | $^{12}\mathrm{C}(K^-, K^+)^{12}_{\Xi}\mathrm{Be}$ | 1.8 GeV/c | 厚い炭素標的 ($1.8\ \mathrm{g/cm^3}$, 厚さ 52 mm) |
| `pik_c12lambda` | $^{12}\mathrm{C}(\pi^+, K^+)^{12}_\Lambda\mathrm{C}$ | 1.4 GeV/c | 炭素 ($1.8\ \mathrm{g/cm^3}$, 厚さ 5.6 mm) |

```bash
# 登録されているプリセット一覧をコマンドで確認
python3 scripts/run_k18_pipeline.py --list-reactions
```

---

## 3. 手動でステップごとに動かす場合

シミュレーションデータが既にある場合や、個別の解析ステップをデバッグ・調整したい場合は、各バイナリを直接実行できます。

### Step 1: Geant4 シミュレーション
`k18-geant4` ディレクトリで実行します。
```bash
# 形式: bin/G4S2S [Confファイル] [出力ROOT] (マクロファイル)
bin/G4S2S param/conf/k18_6381_kpi_c12lambda_example.conf g4_output.root test.mac
```

### Step 2: K18 ビームライントラッキング
`k18-analyzer` ディレクトリで実行します。
```bash
# 形式: bin/DstK18TrackingGeant4 [Confファイル] [入力G4_ROOT] [出力ROOT]
bin/DstK18TrackingGeant4 runmanager/runlist/k18_6381_kpi_c12lambda_k18.conf g4_output.root k18_output.root
```

### Step 3: S-2S スペクトロメータトラッキング
`k18-analyzer` ディレクトリで実行します。
```bash
# 形式: bin/DstS2sTrackingGeant4 [Confファイル] [入力G4_ROOT] [出力ROOT]
bin/DstS2sTrackingGeant4 runmanager/runlist/k18_6381_kpi_c12lambda_s2s.conf g4_output.root s2s_output.root
```

### Step 4: ミッシングマス解析 (`DstMissingMassGeant4`)
プリセットフラグ `--preset` を指定するだけで、質量・標的形状・Bethe-Bloch 補正パラメータが自動適用されます。
```bash
bin/DstMissingMassGeant4 \
  --preset=kpi_c12lambda \
  --g4=g4_output.root \
  --k18=k18_output.root \
  --s2s=s2s_output.root \
  --root=missing_mass.root
```

> **Tips: 特定パラメータの上書き**
> プリセットを使用しながら、特定のパラメータ（例: 標的厚みやカット角度）だけ変更したい場合は、該当オプションを追加するだけで上書きされます：
> ```bash
> bin/DstMissingMassGeant4 --preset=kpi_c12lambda --target-thickness-mm=10.0 ...
> ```

---

## 4. ビルド方法

### `k18-geant4` のビルド
```bash
cd ~/pro/k18-geant4
./build.sh
```

### `k18-analyzer` のビルド
```bash
cd ~/pro/k18-analyzer
make -f Makefile.org dst
```

---

## 5. 設定方針（過剰設計の撤廃と現実的ジオメトリの固定）

本環境では、不用意な分解能悪化や非現実的なシミュレーションを防ぎ、誰でも迷わず使えるように以下の設計方針を採用しています：

- **物理的に実在する測定器要素は常時有効（固定）**:
  - ビームパイプ（`K18BeamPipe`）、真空窓（`K18VacuumWindows`）、BFTファイバー材質（`K18BFTRealisticMaterial`）、遮光シート（`K18BFTBlackSheet`）、ドリフトチェンバー薄膜（`K18BcOutFilms`）などは、現実の設計値で常に組み込まれます。
  - これらの ON/OFF トグルは廃止されており、古い設定ファイルに記述があっても警告を出して安全に既定の現実的ジオメトリが使用されます。
- **過剰な後方互換性やAI専用メタツールの排除**:
  - `config_lifecycle.yml` や移行スクリプト、未知キーによる即時クラッシュ例外は撤廃されました。
  - ユーザーが変更すべき意味のあるパラメータ（ビーム運動量、標的材質・寸法、磁場スケール等）のみを設定ファイルにすっきりと記述できます。
