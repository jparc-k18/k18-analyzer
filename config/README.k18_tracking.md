# K18 Geant4 trackingの最小設定

新しい実行設定は `config/k18_tracking_minimal.conf` から作る。BC3/BC4の位置分解能は
DCGEOの各層の `Res` 列（Gaussianのσ、mm）だけで指定する。hitのsmearと直線fitの
誤差が同じ値を使う。`Res` は有限かつ正でなければ起動を拒否する。

入力は現行K18 Geant4が保存する `K18BC` の面ローカル座標である。`Vx()` がすでに
wire方向の座標なので、追加の傾斜射影を選ぶ設定はない。wireをtruth位置で決め、
`DL = Gaus(abs(x-wire), Res)` を作り、従来の `-0.5 < DL < 1.8 mm` gateをかける。
負DLを一律に絶対値化する変更、gate、最小hit数、χ² cutの変更はしていない。

| 設定 | 役割 |
|---|---|
| `K18TrackingConfigVersion: 2` | 固定の設定形式識別子。無版の旧ideal configを新版で誤実行しないために必須 |
| `UNPACK`, `DIGIT`, `CMAP` | 既存の入出力・detector定義 |
| `DCGEO` | detector配置、wire情報、BC各層の位置分解能 |
| `K18TM`, `PK18` | G4の磁場条件に対応した輸送行列と中心運動量（GeV/c） |
| `USER` | BcOut trackの最小hit数などの共通解析条件 |
| `K18BFTAcceptedAbsPdg` | BFTとBH1/BH2で共通のbeam粒子種。正整数、必須 |
| `K18BFTPositionSmearSigma` | BFT cluster位置へ加えるσ（mm）。0以上、必須 |

最小例は10 keys。乱数seedを指定する場合は `K18HitSmearSeed` を加える（既定20260716、
0--4294967295の整数）。各入力entryごとにBCとBFTの独立streamを作る方式は従来通り。
`K18BFTClusterGap` は必要な場合だけ指定する（既定5 mm）。BFTのσはcluster位置の
応答であり、DCGEOにある個々のfiber面の分解能へ機械的に置き換えてはいけない。

K18の正本運動量はD2Uの `p_3rd`。native RKの物理設定は別の診断出力へ作用する。
数値計算の既定値や、使っていないfull-fit用設定をproduction configへ複写する必要はない。
`chisqrK18` はBcOut直線fitの換算χ²で、保存前に `>20` のtrackを除去する。

## 実行前確認

```bash
make -f Makefile.org -j8 bin/DstK18TrackingGeant4
bin/DstK18TrackingGeant4 --check-config config/k18_tracking_minimal.conf
python3 scripts/config_lifecycle.py --application DstK18TrackingGeant4 config/k18_tracking_minimal.conf
```

`--check-config` は設定の構文・キー・数値を検査する。ROOT入力や出力は開かない。
DCGEOの実ファイル内容と分解能は通常起動時に検査する。重複、未知キー、非有限値、
必須値の欠落を黙って既定値へ読み替えない。

## 古い設定とバイナリ

K18では `G4DCSmearResolutionScale`、`G4DCSmearSignedPosition`、
`G4DCUseTiltedReadout`、`G4DCUseTiltedReadoutAngleFactor`、
`G4DCSmearResolutionScaleSdcIn/Out` を拒否する。これらを消しても同じ応答になる
旧scale=1・追加射影なし・従来DL方式だけ、次のoffline移行ができる。

```bash
python3 scripts/migrate_k18_tracking_config.py old.conf --output new.conf
bin/DstK18TrackingGeant4 --check-config new.conf
```

移行は元ファイルを変更せず、既存の出力configも上書きしない。scale=0/別倍率、
signed-position、追加射影の比較は移行対象外。旧commitと入力を隔離して再現する。
S-2S側の同名設定は今回の廃止対象ではない。

新形式のconfigには必ず新版binaryを使う。旧binaryは形式識別子を検査せず、BCをsmearしない可能性がある。
