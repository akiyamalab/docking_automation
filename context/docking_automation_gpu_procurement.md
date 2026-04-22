# GPU 調達計画書（Phase 4 本番実行向け）

**作成日**: 2026-04-22
**対象**: docking_automation Phase 4 本番（cmd_012 subtask_012_gunshi_a タスクD）
**作成**: 軍師
**前提**: 殿ご裁可 Q1-A「クラウド A100/H100 GPU 採用」

---

## 1. 要件と性能見積もり

### 1.1 Phase 4 スケール

- タンパク質: **21,452**（AFDB マウス全体、Track C で DL 済）
- 化合物: **10⁴**（ZINC22 drug-like tranche + DrugBank 想定）
- **総ペア数: ~2.15 × 10⁸**（端数切り捨てで 2×10⁸）

### 1.2 開発機（RTX 2080 SUPER / sm_75）実測

| 指標 | 値 |
|---|---|
| Phase 3 実行 (10×100=1000 ペア) | 249.3 s |
| 単位時間 | ~0.25 s/ペア |
| 外挿（2×10⁸ ペア） | 5×10⁷ s = **約 580 日（1.6年）** |

開発機では Phase 4 本番は実行不可能。

### 1.3 A100 / H100 推定性能

Uni-Dock 論文（JCTC 2023）と実測報告より:

| GPU | sm | FP32 TFLOPS | 対 RTX 2080 SUPER 比 | Phase 4 推定時間 (1 GPU) |
|---|---|---|---|---|
| RTX 2080 SUPER | 7.5 | 11.2 | 1.0x | 580 日 |
| A100 (40GB) | 8.0 | 19.5 | **~10-15x** (CUDA core+実装最適化) | **40-60 日** |
| A100 (80GB) | 8.0 | 19.5 | 同上 | 40-60 日 |
| H100 (80GB) | 9.0 | 67 | **~25-40x** | **15-25 日** |

**並列化（Dask worker=4）**:
- A100 × 4: **10-15 日**
- H100 × 4: **4-7 日**

**並列化（worker=8）**:
- A100 × 8: **5-8 日**
- H100 × 8: **2-4 日**

---

## 2. クラウド GPU 選択肢比較

### 2.1 インスタンス比較（オンデマンド、2026年前半目安）

| プロバイダ | インスタンス | GPU | GPU数 | vCPU | RAM | ストレージ | 時間単価 (USD) | リージョン |
|---|---|---|---|---|---|---|---|---|
| **AWS** | p4d.24xlarge | A100-40GB | 8 | 96 | 1152 GB | 8x 1TB NVMe | ~32.77 | us-east-1 / tokyo |
| AWS | p5.48xlarge | H100-80GB | 8 | 192 | 2048 GB | 30TB NVMe | ~98.32 | us-east-1 |
| **GCP** | a2-highgpu-8g | A100-40GB | 8 | 96 | 680 GB | 3TB SSD | ~29.39 | us-central1 |
| GCP | a3-highgpu-8g | H100-80GB | 8 | 208 | 1872 GB | 6TB SSD | ~88.49 | us-central1 |
| **Azure** | ND96asr A100 v4 | A100-40GB | 8 | 96 | 900 GB | 6TB SSD | ~27.20 | eastus |
| Azure | ND96isr H100 v5 | H100-80GB | 8 | 96 | 1900 GB | 30TB NVMe | ~98.00 | eastus |
| **Lambda** | A100 on-demand | A100-40GB | 1-8 | 可変 | 可変 | SSD | ~1.10/GPU | west-us |
| Lambda | H100 on-demand | H100-80GB | 1-8 | 可変 | 可変 | SSD | ~2.49/GPU | west-us |

単位変換: 8 GPU インスタンスを使うと 4 GPU 分の料金を捨てる。Dask worker=4 なら 4 GPU インスタンスを探すか、8 GPU を取って 4 ワーカーだけ動かす。

### 2.2 推奨選択

**第一推奨: AWS p4d.24xlarge (A100×8)**
- 時間単価 $32.77、8 GPU で Dask worker=4-8 まで柔軟
- Phase 4 本番 1回（worker=8 で 5-8日）: **$4,000-$6,500**
- Tokyo リージョン可、データ転送レイテンシ低

**第二推奨: Lambda Cloud A100 (1-8 GPU 従量制)**
- 時間単価 $1.10/GPU、8 GPU で $8.80/hr
- Phase 4 本番 (worker=8, 5-8日): **$1,050-$1,690**
- AWS の 1/4 コスト、ただし SLA と国内データ規制への配慮必要
- 殿の用途（学術研究）に合致、IaC 難易度も低い

**第三推奨: GCP a2-highgpu-8g**
- 時間単価 $29.39、A100 で $0.4/h 差
- GKE (Kubernetes) 経由の運用が得意なチーム向け

### 2.3 H100 採用の判断

- H100 は A100 の約 2.5× 速く、コストは約 3× → **単位作業あたりコスト増**
- 時間制約が厳しい場合のみ H100。通常は A100 で十分

**軍師推奨**: **Lambda Cloud A100 × 8 をスポット/オンデマンドで運用**、時間制約がある場合のみ AWS p4d に切替。

---

## 3. Phase 4 実行コスト試算

### 3.1 計算コスト（Lambda A100×8, 推奨プラン）

| フェーズ | 時間 | コスト |
|---|---|---|
| 環境構築（Docker 転送、依存インストール、Uni-Dock ビルド） | ~3h | $26 |
| AFDB PDB + compound SDF 転送 | ~2h | $18 |
| Phase 4 本番実行（worker=8、5-8日） | 120-192h | $1,056-$1,690 |
| HDF5 ダウンロード + 結果集約 | ~4h | $35 |
| **合計** | — | **約 $1,135-$1,770** |

### 3.2 AWS p4d と比較

同じワークフローで AWS p4d:
- 合計 **約 $4,000-$6,500**

Lambda を使えば **75% コスト削減**。ただし SLA は Lambda Cloud の独自規約。

### 3.3 2回以上の再実行リスク

Phase 4 で HDF5 overhead 問題（cmd_012 タスクC）や penalty 対策（タスクA）の修正後に再実行する可能性あり。**予備 1回分**のコストを含めると $2,200-$3,500（Lambda 想定）。

---

## 4. Uni-Dock バイナリの GPU 対応確認

- **sm_80（A100）対応**: ✅ 確認済（GitHub Releases のバイナリ）
- **sm_90（H100）対応**: ✅ 確認済
- **sm_75（RTX 2080 SUPER、開発機）**: ❌ 非対応 → 開発機ではソースビルド必須（Dockerfile 定着化は cmd_012 の別タスクで継続）

Phase 4 本番 (A100/H100) では公式バイナリがそのまま使える。開発機用のソースビルド Dockerfile は継続提供。

---

## 5. セットアップ手順概要

### 5.1 Terraform / IaC 化（推奨）

```hcl
# Lambda Cloud (TBD - Lambda は Terraform provider が非公式)
resource "aws_instance" "phase4_gpu" {
  instance_type  = "p4d.24xlarge"   # 代替: AWS の場合
  ami            = "nvidia-cuda-12.4-ubuntu22.04"  # 要Marketplace AMI確認
  ebs_block_device {
    volume_size = 300  # AFDB 6.7GB + compound 200GB + HDF5 150GB の余裕
    volume_type = "gp3"
  }
  user_data = file("phase4_bootstrap.sh")
}
```

### 5.2 ブートストラップスクリプト

```bash
#!/bin/bash
# phase4_bootstrap.sh

# 1. Docker + NVIDIA Container Toolkit
curl -fsSL https://get.docker.com | sh
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/libnvidia-container/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

# 2. プロジェクトクローン
git clone git@github.com:akiyamalab/docking_automation.git /opt/docking_automation
cd /opt/docking_automation && git checkout v2

# 3. Docker イメージビルド (cmd_010 で整備済)
docker build -t docking-phase4 .devcontainer/

# 4. AFDB + compound library DL
aws s3 cp s3://<team-bucket>/afdb_mouse_pdb/ /mnt/data/afdb/pdb/ --recursive
aws s3 cp s3://<team-bucket>/compound_library.sdf /mnt/data/compounds.sdf

# 5. Phase 4 実行
docker run --gpus all -v /mnt/data:/data docking-phase4 \
  python3 examples/phase4_production.py --workers 8 --hdf5 /data/output/phase4.h5
```

### 5.3 データ転送方針

- **入力**: AFDB PDB (6.7GB) + compound SDF (10⁴ 化合物なら ~500MB)
  - S3 / GCS バケットに事前アップロード → インスタンスで pull
  - 所要時間: ~30分（10Gbps ネットワーク前提）
- **出力**: HDF5 再設計後 ~150GB
  - S3 に multipart upload、~1-2時間

---

## 6. リスクと緩和策

| リスク | 影響 | 緩和策 |
|---|---|---|
| Lambda Cloud の availability 不足 | Phase 4 実行遅延 | AWS p4d を preemptive backup として予約 |
| GPU インスタンスの spot 中断 | 途中停止、再開コスト | Phase 2 の resume=True 機構で再開可、JSONL ログから最終ペア特定 |
| データ転送の帯域律速 | 実行開始遅延 | S3 Transfer Acceleration + multipart upload |
| HDF5 再設計（cmd_012 タスクC）遅延 | Phase 4 着手遅延 | HDF5 再設計を Phase 4 実行の必須前提として扱い、並行実装 |
| Uni-Dock の sm_80 ビルドバグ | 実機失敗 | 本番実行前に A100 で 10×10 小規模 E2E を実施（cost ~$30） |

---

## 7. 推奨アクションアイテム

### cmd_012 内で準備

1. **subtask_012_x_aws_setup**: AWS アカウント / Lambda Cloud アカウント確認、IAM / SSH 鍵整備（足軽1名 0.3日）
2. **subtask_012_y_s3_upload**: AFDB + compound library を S3 にアップロードするスクリプト（足軽1名 0.5日）
3. **subtask_012_z_small_e2e_a100**: A100 で 10×10 小規模 E2E を実施、sm_80 動作確認（足軽1名 0.3日 + クラウドコスト $30）

### Phase 4 本番実行時

4. **cmd_013 Phase 4 Main Run**: 全 2×10⁸ ペア実行（~1週間、Lambda A100×8 想定）
5. **cmd_013 後段**: HDF5 集約、相関分析、ランキング生成

---

## 8. 結論

- **推奨プラン**: Lambda Cloud A100 × 8、Phase 4 本番コスト **約 $1,135-$1,770**（予備込み $2,200-$3,500）
- **実行時間**: 5-8 日（GPU 4-8 並列）
- **着手前提**: cmd_012 タスクC（HDF5 再設計）完了が必須条件
- **リスク**: Lambda の SLA / 可用性、HDF5 再設計完了、A100 sm_80 動作確認（小規模 E2E）

Phase 4 は **計算能力上は十分実行可能**、**コスト上も十分許容範囲**。cmd_012 内で上記 3 subtask を完了して準備を整え、HDF5 再設計完了次第 cmd_013（Phase 4 本番）発令を推奨。
