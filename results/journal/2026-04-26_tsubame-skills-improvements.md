# tsubame_skills wrapper 改善 (2026-04-26)

## 依頼

> 今、自分の TSUBAME ポイントがどのくらい残っているのか、ストレージがどうかを確認したいです。
> グループは、ポイントは --group に指定するもの、ストレージはターゲットディレクトリの pwd -P から取得する必要があります。
> それも、tsubame skills に入れたいです。
> (続いて) home 領域にそれが配置されるのは良くないですね。workdir 側に配置できますか
> (続いて) workdir に隠しディレクトリを作成して配置
> (続いて) results を tsubame_skills ではなく docking_automation 側に戻してほしい

## 実施

### 1. `tsubame account` 動詞の追加
- **server (`templates/claude-wrapper.sh`)**: `account <group>` verb を追加。`/apps/t4/rhel9/hpe/ptl/info_com/user/bin/t4-user-info` を絶対パスで呼ぶ (forced command の制限 PATH に対応)
- **client (`bin/tsubame`)**: `tsubame account [--group <name>]` 動詞、`--group` 省略時は `TSUBAME_GROUP` 利用
- **出力**: `t4-user-info group point -g <group>` (ポイント残高) + `t4-user-info disk group/home` (ストレージ使用量)
- 反映: `scripts/init-remote.sh --ssh-host t4` で wrapper 再デプロイ

### 2. sif の home → workdir 移行
- **問題**: `$HOME/apptainer/tsubame-env.sif` (1.4 GB) が home 25 GB quota の 88% 圧迫の主因
- **方針**: workdir 上 `.sif/` (隠しディレクトリ) に移動 = group disk `/gs/bs` 上で home quota を消費しない
- **変更**:
  - `bin/tsubame push` の `--exclude` に `.sif/` 追加 (push の `--delete` 対象外)
  - `jobs/build-apptainer.sh` 出力先を `.sif/tsubame-env.sif` に
  - 全 dock/prep スクリプトの `SIF=` を `.sif/tsubame-env.sif` に統一
  - 一時ジョブ `migrate-sif.sh` で既存 sif を mv (再ビルド回避)、`$HOME/apptainer` 削除

### 3. pull 先のプロジェクト分離
- **問題**: `tsubame pull` が `tsubame_skills/results/` に書いていたが、結果は本来プロジェクト (`docking_automation/`) のもの
- **変更**:
  - `bin/tsubame` に `TSUBAME_RESULTS_DIR` env サポート + `--to <path>` フラグ追加
  - `tsubame_skills/.tsubame.conf` に `TSUBAME_RESULTS_DIR=../docking_automation/results` を設定
  - 既存 `tsubame_skills/results/` 内 1.1 GB を `docking_automation/results/` に mv

### 4. 結果ジャーナル設計 (本ファイル含む)
- 各タスク単位で `docking_automation/results/journal/<YYYY-MM-DD>_<slug>.md` を Claude が自動生成する規約
- ユーザの依頼プロンプト + ジョブ ID 表 + 主要ファイル + 学び を構造化
- index `RESULTS.md` に時系列で集約

## 主要ファイル変更

| ファイル | 変更点 |
|---|---|
| `bin/tsubame` | `account`/`pull --to` 動詞、`TSUBAME_RESULTS_DIR` env、push の `.sif/` exclude |
| `templates/claude-wrapper.sh` | `account` verb (t4-user-info 絶対パス) |
| `jobs/build-apptainer.sh` | sif 出力先 `.sif/tsubame-env.sif` |
| `jobs/{dock,prep}-*.sh` | `SIF=.sif/tsubame-env.sif` に統一 |
| `.tsubame.conf` | `TSUBAME_RESULTS_DIR=../docking_automation/results` |
| `.claude/skills/tsubame/SKILL.md` | 新動詞・運用ルール記載 |
| `.devcontainer/Dockerfile` (docking 側) | fpocket を conda 経由で追加 (CLAUDE.md 「install 後は dockerfile 追記」ルール準拠) |

## 学びと変更

- **forced command のセキュリティと利便性のバランス**: 動詞ベースの allowlist 構造は新機能追加が安全 (任意 shell 実行を許さない)。代わりに wrapper 改修 → init-remote.sh 再実行が必要なオーバーヘッドあり
- **t4-user-info は標準 PATH に無い** (TSUBAME 4 の謎): `/apps/t4/rhel9/hpe/ptl/info_com/user/bin/t4-user-info` という個別パッケージ。今後類似ユーティリティを叩く場合は **「ユーザに正確なパスを聞く」** のが最速
- **隠しディレクトリ慣例**: `.sif/`, `.git/` 等の hidden naming は「内部用」のシグナルとして有効。ただし rsync の `--delete` 対象外には**自動でならない** ので明示 `--exclude` 必須
- **疎結合 / プロジェクト境界**: tsubame_skills は wrapper のみ、結果はプロジェクト本体に置く設計が整理しやすい

## 開発予定 (未着手)

- `tsubame submit --after <jobid>` (ジョブ依存関係、qsub の `-hold_jid` 透過)。多段パイプライン (prep → dock) を 1 度の submit で連鎖実行できるようにする。forced command への引数許可拡張要
