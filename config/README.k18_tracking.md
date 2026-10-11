# K18 Geant4 trackingの最小設定

新しい実行設定は `config/k18_tracking_minimal.conf` から作る。BC3/BC4の位置分解能は
DCGEOの各層の `Res` 列（Gaussianのσ、mm）だけで指定する。hitのsmearと直線fitの
誤差が同じ値を使う。`Res` は有限かつ正でなければ起動を拒否する。

S2Sもsmearingのσには各層のDCGEO `Res` を直接使う。
`G4DCSmearResolutionScale` と `G4DCSmearResolutionScaleSdcIn/Out` は全DCで廃止し、
設定例から削除した。倍率を1として掛ける処理も残していない。
従来倍率1の設定例の応答は変わらない。倍率省略でsmearingなしとなる旧既定動作は廃止した。
DCGEOの値、fit重み、DL gateは変更していない。BFTのcluster位置σは別の設定なので変更しない。

入力は非回転の敏感面を持つK18 Geant4が保存する `K18BC` のチェンバー局所座標である。
S2Sと共通に `s = Vx()*cos(tilt) + Vy()*sin(tilt)` を常に計算する。
枠は回さず、DCGEOのワイヤー角だけを読み出し投影へ使う。切替パラメータはない。
wireをtruth読み出し位置で決め、`DL = Gaus(abs(s-wire), Res)` を作り、
従来の `-0.5 < DL < 1.8 mm` gateをかける。
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

## 現行入力だけを使う

新規実行には Version 2 の config と対応する新版 binary を使う。
旧 `G4DC*` のK18応答switch・倍率・追加射影設定は削除済みであり、指定すると拒否する。
S2Sの `G4DCSmearSignedPosition` も廃止した。SDCもBCと同じく `DL = Gaus(abs(s-wire), Res)` でsmearする。

ROOTの `g4s2s` UserInfoには `G4DCReadoutFrame=chamber-local-v1` が必要。
識別情報がない入力は出力作成前に拒否する。旧座標のROOTは現行Geant4で再生成する。
`--check-config` はROOTを開かないため座標規約の検査は通常起動時に行う。
旧結果を再現するときだけ、当時のソフトウェア・設定・ROOTを一組で隔離して使う。
