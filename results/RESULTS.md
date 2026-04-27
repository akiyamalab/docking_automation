# Results Journal

TSUBAME 上の作業記録。各エントリは独立した md として `journal/` 配下に保存。
新しい依頼が完了したら Claude が自動的に追加する規約 (詳細は `claude_tsubame/.claude/skills/tsubame/SKILL.md`)。

| 日付 | タイトル | 概要 | 詳細 |
|---|---|---|---|
| 2026-04-26 | tsubame_skills wrapper 改善 | `tsubame account` 動詞、sif の `.sif/` 移行、pull 先のプロジェクト分離、結果ジャーナル設計 | [→](journal/2026-04-26_tsubame-skills-improvements.md) |
| 2026-04-26 | 1000×100 virtual screening | 1000 受容体 × 100 リガンド docking、999,509 poses 生成 (99.95%)。fpocket box / array 競合 fix の発見 | [→](journal/2026-04-26_1000x100-virtual-screen.md) |
| 2026-04-25 | MPS × apptainer 互換性調査 | TSUBAME 4 + apptainer 1.3.6 + MPS は post-RPC 段階で hang。実用解は MPS なし運用 | [→](journal/2026-04-25_mps-investigation.md) |
| 2026-04-25 | TSUBAME ドッキングパイプライン構築 | 10×10 / 100×100 を v1 → v2 へ移行、DMS キャッシュで 100×100 を 10 分化 | [→](journal/2026-04-25_docking-pipeline-establishment.md) |

## 主要成果物

- **dock-screen.7264678.results.tar.gz** (328 MB): 1000×100 docking の最終結果 999,509 poses
- **dock-100x100.7262439.results.tar.gz** (33 MB): 100×100 v2 + DMS cache 完走 (10 min)
- **dock-100x100.7261784.results.tar.gz** (29 MB): 100×100 v1 serial (15 min, baseline)
- TSUBAME workdir 永続: `~/workspace/.claude-workdir/results/receptors_dms/` に DMS 1100 件 (1.2 GB)、`.sif/tsubame-env.sif` (1.4 GB unidock2 + cuda 13)
