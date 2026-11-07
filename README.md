# SINCHO-Gen-GUI

SINCHO-Gen 0.1.0用のStreamlit GUIアドインです。入力PDBの確認、条件YAML生成、MD・SINCHO・ChemTS・AA Score結果の可視化を行います。

GUIはSINCHO-Genの入出力アドインとして同一バージョンで管理します。配布版は両方とも`0.1.0`です。GUIイメージは同じバージョンタグのSINCHO-Genイメージを基盤にビルドします。

## 対応するSINCHO-Gen

- 公開モデル版: `sincho-gen:0.1.0`
- 計算ドライバ: `/SINCHO-Gen/extend_driver.py`

GUIは計算ドライバを複製しません。YAML生成後は、同じコンテナ内にあるSINCHO-Gen本体を実行します。

## GUIイメージの作成

```bash
cd /path/to/SINCHO-Gen-GUI
./build-image.sh
```

作成されるイメージ:

```text
sincho-gen-gui:0.1.0
```

```bash
SINCHO_GUI_BASE_IMAGE=sincho-gen:0.1.0 \
SINCHO_GUI_IMAGE=sincho-gen-gui:0.1.0 \
./build-image.sh
```

## 起動

計算ファイルを保存するディレクトリから起動します。

GUIのみの利用にはGPUは不要です。通常起動ではGPUをコンテナへ割り当てません。

```bash
cd /path/to/job
/path/to/SINCHO-Gen-GUI/run-gui.sh
```

ブラウザで次を開きます。

```text
http://localhost:8501
```

ポート、イメージは環境変数で変更できます。同じコンテナ内でGPUを使う計算も実行する場合は、`SINCHO_GEN_GPU` にGPU番号を指定して起動してください。その場合はホスト側のNVIDIAドライバとDockerのGPU実行環境が必要です。

```bash
SINCHO_GUI_HOST_PORT=8502 \
SINCHO_GEN_GPU=1 \
SINCHO_GUI_IMAGE=sincho-gen-gui:0.1.0-private \
/path/to/SINCHO-Gen-GUI/run-gui.sh
```

## 計算開始

GUIのSummary画面でYAMLを保存した後、コンテナへ入ります。
GPUを必要とする計算を実行する場合は、あらかじめ上記の `SINCHO_GEN_GPU` を指定してコンテナを起動してください。

```bash
docker exec -it sincho-gen-gui /usr/sbin/gosu dockeruser bash
python /SINCHO-Gen/extend_driver.py conditions.yaml
```

途中再開:

```bash
python /SINCHO-Gen/extend_driver.py conditions.yaml --skip-md
python /SINCHO-Gen/extend_driver.py conditions.yaml --skip-sincho
python /SINCHO-Gen/extend_driver.py conditions.yaml --skip-chemts
```

## 停止・削除

```bash
docker stop sincho-gen-gui
docker rm sincho-gen-gui
```
