# PRESTO用 ZIP 配布パッケージ（2026-11-19）

このZIPは、研究者が自身の方法でPRESTOへ転送・展開するためのものです。
ChatGPT はPRESTOへログイン・ジョブ投入していません。

## 中身

`latest_main_stage/` : ueda1984/ED2 の**配布生成時点**の main の全ソース（tar.gz）、
PR #36 の QuSpin/XDiag 比較コード、PR #37 のPRESTO PBSランチャーを、それぞれ
異なるコミットに固定したステージです。SHA256SUMSで検証できます。

`experimental_cubic_PR21/` : 未マージの立方格子専用高速化PR #21 の全ソースと、
既存のPRESTO立方格子ベンチマーク用PBSスクリプト。**mainとは別の実験版**で、
現行mainとして報告してはいけません。

`DELIVERY_MANIFEST.json` : 収録したmain/PR36/PR37/PR21の完全SHAと、
全ファイルのチェックサム。

## PRESTOでの展開

ZIPをLinuxの任意の新しいディレクトリに展開してください。
次の例ではそのディレクトリを `$BUNDLE` とします。

```bash
BUNDLE="$HOME/ED2_PRESTO_DELIVERY_20261119"
# ZIPの中身が展開されたルートをBUNDLEに指定
cd "$BUNDLE"
sha256sum -c DELIVERY_SHA256SUMS
(cd latest_main_stage && sha256sum -c SHA256SUMS)
(cd experimental_cubic_PR21 && sha256sum -c ed2-source.sha256)
python3 latest_main_stage/preflight.py --stage "$BUNDLE/latest_main_stage" --mode stage
cat latest_main_stage/main-commit.txt
cat experimental_cubic_PR21/source-commit.txt
```

**最新版mainの注意:** ZIPのmainは配布生成時点に固定されています。
別チャットでmainが後からマージ更新された場合は、このZIPの実行で
「更新後の最新版」を名乗れません。GitHubに接続可能な環境で
`git ls-remote origin refs/heads/main` と
`latest_main_stage/main-commit.txt` を照合してください。
不一致なら新規ZIPを作成してください。古いZIPを黙って上書きしないでください。

## 1000サイト・3マグノンの本実験を開始する前に

配布時点のmainは汎用CSRを含みますが、PR #13～#21由来の
非可換立方格子群証明・anchored基底・direct-hop表、
`ED2_LANCZOS_TWO_STAGE`、`ED2_PHASE_TIMING`、
`ED2_DIAGNOSTICS` などの専用実験インターフェースが未統合です。

したがって latest_main_stage/preflight.py が
`SOURCE_NOT_READY` を返す場合、**mainとQuSpinの
同一A1gセクター性能比較のPBSジョブを投入してはいけません。**
この失敗は正常な安全装置です。

## 実験的なPR21版で既存のPRESTO実験を継続する場合

PR21には最適化コードと既存PBSドライバーが入っています。
PR21は未マージの実験版として結果とソースを必ず区別してください。
PBSジョブ投入には従来のPRESTO管理者による独占利用・資源ポリシーの確認、
`qstat` の照合、実際の物理モデルと直接Ritz残差の確認が必要です。
ジョブはこのZIPの展開だけでは**投入されません**。

まず、以前のPR16/PR21 PRESTOテストと同様に、
`ED2_STAGE` に `experimental_cubic_PR21` の絶対パス、
`ED2_RESULTS_ROOT` に書き込み可能で空きのある結果ディレクトリ、
`ED2_EXCLUSIVE_CONFIRMED=1` に管理者確認済みの場合のみ設定してください。

例（確認前には**実行しない**）：
```bash
mkdir -p "$HOME/ED2_EXPERIMENT_RESULTS"
qsub -v "ED2_STAGE=$BUNDLE/experimental_cubic_PR21,ED2_RESULTS_ROOT=$HOME/ED2_EXPERIMENT_RESULTS,ED2_PRESET=cubic-main,ED2_SIZES=10,ED2_THREADS=1,ED2_SECTORS=oh-a1g,ED2_WARMUPS=0,ED2_REPEATS=1,ED2_TIMEOUT=1800,ED2_EXCLUSIVE_CONFIRMED=1" \
  "$BUNDLE/experimental_cubic_PR21/presto_cubic_exclusive.pbs"
```

上記は**PR21の実験的ED2**を計算する例で、QuSpinとの同一ノード
直接比較や、最新mainのベンチマークではありません。

## mainが統合後にREADYとなり、QuSpin 1.0.1環境を準備できた場合

まず `latest_main_stage/preflight.py --mode submit --repo ...` で
最新remote mainと一致することを確認します（ネットワークが必要）。
PRESTOのPBS割当権限・排他実行条件を確認してから、
`ED2_STAGE=$BUNDLE/latest_main_stage` と
`ED2_RESULTS_ROOT=...` を設定し、
`latest_main_stage/launch.pbs` を**手動で**投入してください。

QuSpin比較にはcompute nodeからアクセスできるPython 3.10以上の
QuSpin **1.0.1** 仮想環境が必要です。具体的なコマンド・
`ED2_PROGRAMS=ed2-cubic+quspin` の指定は、
`latest_main_stage/launch-tools.tar.gz` の
`tools/presto_cross_library/README.md` にあります。

## データの解釈

すべての数値結果を出力ハッシュ・プログラムcommit・CPU・BLAS・
OpenMPスレッド数・配置・直接残差と一緒に保存してください。
異なるセクター間の経過時間を直接「高速化倍率」として報告しないでください。
共同CIの値とPRESTOの排他計算資源での値を混在させないでください。
CPC投稿用の優位性主張には同一ノード・同一物理セクターでの反復測定と
管理者による独立した割当・排他状態の確認が必要です。

GitHub:
https://github.com/ueda1984/ED2/pull/36
https://github.com/ueda1984/ED2/pull/37
https://github.com/ueda1984/ED2/pull/21
