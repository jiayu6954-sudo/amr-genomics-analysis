# 工作流文档 — AMR 基因组学分析管道
> **Pipeline Workflow Reference** · 版本 1.3 · 2026-04-19

---

## 总览

```
[NCBI GenBank]
      │ Entrez API (esearch + esummary)
      ▼
┌─────────────────────────────────────────────────────────────┐
│  01_download.py   [P1修复分页bug: None→10,000 retmax]       │
│  获取中国临床K.pneu基因组 → data/raw/{accession}/           │
│  677 UIDs检索 → 492 GFF可用                                │
└───────────┬─────────────────────────────────────────────────┘
            │ download_status.tsv
            ▼
┌─────────────────────────────────────────────────────────────┐
│  00_data_audit.py   ★P1新增                                │
│  7项完整性检查：gzip/GFF3/FASTA/大小/N50/CDS/重复检测      │
│  677基因组全量审计 → 492 PASS / 181 MISSING_GFF / 4 QC失败 │
│  0注册表重复 + 0 BioSample重复 + 0 gzip损坏                │
└──────────────────────────┬──────────────────────────────────┘
                           │ audit_report.tsv + audit_summary.json
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  02_validate.py                                             │
│  质控门：大小/N50/CDS/可读性检查                            │
│  492 PASS                                                   │
└──────────────────────────┬──────────────────────────────────┘
                           │ genome_manifest.tsv (PASS only)
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  03_amr_scan.py                                             │
│  Tier1: AMRFinderPlus TSV → Tier2: GFF3正则扫描             │
│  217 carbapenem hits / 200基因组耐药 (40.7%)                │
└──────────────────────────┬──────────────────────────────────┘
                           │ amr_hits.tsv
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  04_is_context.py                                           │
│  IS元件侧翼结构分析（±10kb窗口）                           │
│  2,985 IS-AMR对 / 复合转座子率=94.6%                        │
└──────────────────────────┬──────────────────────────────────┘
                           │ is_context.tsv + context_summary.tsv
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  05_stats.py            │  06_figures.py  │  07_manuscript.py│
│  统计检验                │  图表           │  自动草稿        │
└─────────────────────────┴─────────────────┴─────────────────┘
      │ (全队列 IS 负担分析 — 独立分支，可并行运行)
      ▼
┌─────────────────────────────────────────────────────────────┐
│  08_is_burden_all.py                                        │
│  全492基因组IS26/IS6负担定量 → 耐药预测 AUC                 │
│  全队列AUC=0.807 / 2025队列AUC=0.976（见Step 11）           │
└──────────┬───────────────────────────┬───────────────────────┘
           │ IS_unknown蛋白质          │ AMR hit坐标
           ▼                           ▼
┌──────────────────────────┐ ┌─────────────────────────────────┐
│  09_is_hmmer_verify.py   │ │  10_amr_hmmer_verify.py         │
│  PFAM PF01527 HMM验证    │ │  phmmer vs CARD v3.3            │
│  IS_unknown→IS6重分类    │ │  AMR基因蛋白质序列确认          │
│  AUC: 0.952→0.949 ✓     │ │  217/217 (100%) CONFIRMED ✓    │
└──────────────────────────┘ └─────────────────────────────────┘
     ↑ data/db/pfam/PF01527.hmm      ↑ data/db/card/protein_fasta_*.fasta
           │ (提交偏差量化 + 时间分层)
           ▼
┌─────────────────────────────────────────────────────────────┐
│  11_subgroup_analysis.py   ★P1新增                          │
│  A. 年份分层耐药率+IS6负担                                  │
│  B. 2025队列敏感性分析 (n=188, AUC=0.976, δ=0.952)          │
│  C. Pre-2022 vs 2022+ 对比（Fisher OR=48.0, p=2×10⁻⁴⁶）   │
│  D. 复合转座子率跨时代稳定性（88–95%）                     │
└──────────────────────────┬──────────────────────────────────┘
                           │ subgroup_stats.json + subgroup_table.tsv
                           ▼
              [reports/manuscript_submission_v1.md]
              投稿级论文（Nature Communications）
```

---

## 从零开始：完整构建指南

### 前置条件

```bash
# 1. 创建conda环境（llama-env已有，跳过）
conda create -n amr_env python=3.12
conda activate amr_env

# 2. 安装依赖（pyhmmer为Step 09/10序列验证必需）
pip install pandas scipy numpy matplotlib requests pyhmmer

# 3. 确认环境
python -c "import pandas, scipy, matplotlib, requests, pyhmmer; print('OK')"
```

### 目录初始化

```bash
mkdir -p amr_project/{analysis,data/{raw,validated,processed,db/{pfam,card}},logs,figures,reports,docs}
cd amr_project

# 将以下文件放入 analysis/：
# config.py  01_download.py  02_validate.py  03_amr_scan.py
# 04_is_context.py  05_stats.py  06_figures.py  07_manuscript.py
# 08_is_burden_all.py  09_is_hmmer_verify.py  10_amr_hmmer_verify.py
```

### 数据库准备（Step 09/10需要）

```bash
# PFAM PF01527 — 脚本自动下载（无需手动操作）
# 下载到: data/db/pfam/PF01527.hmm

# CARD v3.3 — 需手动下载
# 1. 访问 https://card.mcmaster.ca/download
# 2. 下载 "Protein Homolog Model" (broadstreet-v3.3.0.tar.bz2)
# 3. 解压并复制:
tar -xjf broadstreet-v3.3.0.tar.bz2
cp protein_fasta_protein_homolog_model.fasta \
   amr_project/data/db/card/protein_fasta_protein_homolog_model.fasta
```

### Step 0 — 配置参数

编辑 `analysis/config.py`，最少需修改：

```python
NCBI_EMAIL   = 'your.email@example.com'   # 必须
NCBI_API_KEY = 'your_api_key_here'         # 可选，获取后速率提升10×

# 首次运行建议先用小批量测试
MAX_GENOMES_CHINA = 10   # 改为10测试，确认后改为300或500
```

---

## Step 0 — 数据完整性审计 ★P1新增

```bash
# 快速模式（下载进行中时使用，只读前64KB）
python analysis/00_data_audit.py --fast

# 完整模式（所有下载完成后，完整解压验证）
python analysis/00_data_audit.py
```

**7项检查内容**:

| 检查项 | 说明 |
|--------|------|
| A. 文件存在性 | gff.gz + fna.gz 必须存在 |
| B. gzip完整性 | 完整解压（fast模式仅读前64KB） |
| C. GFF3格式 | ≥9列 / start≤stop / 有效链向 / CDS≥100 |
| D. FASTA格式 | ≥1条序列 / 有效DNA字符 / 大小≥物种阈值×0.8 |
| E. 组装统计 | 基因组大小范围 + N50 + CDS计数 |
| F. 交叉验证 | accession一致性检查 |
| G. 重复检测 | 注册表重复（同accession）+ BioSample重复（同样本多组装）|

**P1审计结果**（677基因组）:

```
注册表重复:  0
BioSample重复: 0
gzip损坏:    0
MISSING_GFF: 181（GCF_前缀FTP路径问题 — 见D-012、L-002）
QC_N50_FAIL: 3
最终通过:    492/677 (72.7%)
```

**输出**:
- `data/processed/audit_report.tsv` — 每基因组一行，含所有检查结果
- `data/processed/audit_summary.json` — 汇总统计

---

## Step 1 — 数据下载

```bash
python analysis/01_download.py
```

**输出**: `data/raw/download_status.tsv`

| 列 | 含义 |
|----|------|
| assembly_accession | GCA_/GCF_ 编号 |
| organism_name | 物种全名 |
| assembly_level | Complete/Chromosome/Scaffold |
| status | OK / MISSING_REQUIRED / ERROR |
| local_dir | 本地存储路径 |
| has_amrfinder | True/False（AMRFinder TSV是否存在） |

**预期结果（300基因组）**:
- OK: ~271
- MISSING_REQUIRED: ~29（GCF_前缀FTP路径问题）
- 运行时间: ~45分钟

**常见问题**:
```
# 问题: 下载中断后重启
# 解决: 脚本自动读取已有的download_status.tsv跳过已完成项

# 问题: HTTP 429 Too Many Requests
# 解决: 在config.py中增大 RATE_LIMIT_DELAY_S (改为1.0)
# 或者: 申请NCBI API key (免费，速率提升至10req/s)

# 问题: 某个基因组的GFF下载失败
# 解决: 检查 logs/download.log，单独重试
```

---

## Step 2 — 数据质控

```bash
python analysis/02_validate.py
```

**输出**:
- `data/validated/genome_manifest.tsv` — 通过质控的基因组（PASS）
- `data/validated/qc_report.tsv` — 所有基因组的QC详情

**质控检查项**:
1. 本地目录存在
2. GFF3.gz 可解压、非空
3. FASTA.gz 可解压、非空
4. 基因组大小在物种阈值内（K.pneu: 4.8–6.5Mb）
5. N50 ≥ 50,000 bp
6. CDS数量 ≥ 3,000
7. GFF3 header中taxid交叉验证（非阻断性）

**预期结果**: 270/271 PASS

---

## Step 3 — AMR基因扫描

```bash
python analysis/03_amr_scan.py
```

**输出**: `data/processed/amr_hits.tsv`

**列说明**:

| 列 | 含义 |
|----|------|
| accession | 基因组编号 |
| gene_name | 标准化基因名（如KPC-2, NDM-5） |
| drug_class | 耐药类型（BETA-LACTAM） |
| drug_subclass | 亚类（CARBAPENEM） |
| contig | 所在contig名 |
| start/stop | 0-based坐标 |
| is_carbapenem | True/False |
| detection_tier | AMRFINDER / GFF_KEYWORD |

**Tier 2 正则模式（已修正版）**:
```python
r'\bNDM-\d'          # 匹配 NDM-1, NDM-5, NDM-13 等
r'\bKPC-\d|blaKPC'   # 匹配 KPC-2 或 blaKPC 基因属性
r'\bIMP-\d|blaIMP'   # 匹配 IMP-4 等（注意：不能用 IMP\b，会匹配rimP）
r'\bVIM-\d|blaVIM'   # 匹配 VIM-1, VIM-2 等
r'OXA-48|OXA-181|OXA-232'
r'carbapenem.{0,30}(resistance|beta-lactamase)'
r'bla[A-Z]{2,5}-\d'  # 通用bla命名（需等位号，排除blaZ等非碳青霉烯基因）
```

> **陷阱警告**: `r'IMP\b'` (case-insensitive) 会匹配 `rimP`、`purH`、`guaB` 等
> 代谢酶基因，产生大量假阳性。务必使用带等位号的模式。

---

## Step 4 — IS元件侧翼分析

```bash
python analysis/04_is_context.py
```

**输出**:
- `data/processed/is_context.tsv` — 每行为一个IS-AMR配对
- `data/processed/context_summary.tsv` — 每基因组汇总

**分类逻辑**:
```
AMR基因位置: ────[AMR_start ─── AMR_stop]────
IS元件检测:  寻找 AMR±10,000bp 内所有IS特征

gap_left  = max(0, AMR_start - IS_stop)   # IS在上游
gap_right = max(0, IS_start  - AMR_stop)  # IS在下游
distance  = min(gap_left, gap_right)       # 取较小值

if distance <= 10,000:
    position = UPSTREAM if gap_left <= gap_right else DOWNSTREAM

flank_class:
  上游≥1 AND 下游≥1  → COMPOSITE_TRANSPOSON
  只有上游            → SINGLE_IS_UPSTREAM
  只有下游            → SINGLE_IS_DOWNSTREAM
  无IS                → NO_IS
```

**IS家族提取优先级**:
1. Tn\d+ 模式（如Tn4401）
2. IS[A-Za-z0-9]+ 模式（如IS26, ISKpn26）
3. 无法识别 → IS_unknown

---

## Step 5 — 统计分析

```bash
python analysis/05_stats.py
```

**统计方法**:

| 指标 | 方法 |
|------|------|
| 耐药率置信区间 | Wilson score 95%CI |
| 复合转座子率置信区间 | Clopper-Pearson 精确 95%CI |
| IS家族分布均匀性 | chi-square 拟合优度检验 |
| IS负担中位数 | 描述统计（中位数+IQR） |

**输出**:
- `stats_summary.tsv` — 全部指标平铺表
- `stats_tables.json` — 结构化JSON（`07_manuscript.py`读取）

---

## Step 6 — 图表生成

```bash
python analysis/06_figures.py
```

**生成6张图**:

| 图 | 内容 |
|----|------|
| fig1_prevalence | A: 耐药率柱状图+CI；B: 基因家族分布 |
| fig2_is_context | IS侧翼分类堆积柱（按基因家族） |
| fig3_is_families | IS家族水平柱状图（Top 10） |
| fig4_is_burden | 耐药基因组IS负担箱线图 |
| fig5_forest | 复合转座子率森林图（按基因家族） |
| fig_combined | 6个面板合并图（论文投稿用） |

每张图同时输出 `.pdf`（向量，投稿用）和 `.png`（300 dpi，预览用）。

---

## Step 7 — 论文草稿生成

```bash
python analysis/07_manuscript.py
```

**输出**: `reports/manuscript_draft_v1.md`

草稿从 `stats_tables.json` 自动填充所有数字，包含：
- 标题 / 作者占位符 / 关键词
- 摘要（~250词）
- 引言 / 方法 / 结果 / 讨论 / 结论
- 声明 / 参考文献（15条）

> 每次重新运行 `05_stats.py` 后，重跑 `07_manuscript.py` 即可自动更新所有数字。

---

## Step 8 — 全队列IS负担分析

```bash
python analysis/08_is_burden_all.py
```

**输入**（自动读取）:
- `data/validated/genome_manifest.tsv` — 492个基因组路径（P1）
- `data/processed/amr_hits.tsv` — 区分耐药/敏感分组

**输出**:

| 文件 | 内容 |
|------|------|
| `data/processed/is_burden_all.tsv` | 492行，每基因组一行；包含n_is_total, n_is6, n_is26, assembly_year, is_resistant |
| `data/processed/is_burden_stats.json` | 全部统计结果（中位数、IQR、Mann-Whitney p、Cliff's δ、AUC、OR、预激株accessions） |
| `figures/fig_is_burden_violin.pdf/png` | IS26/IS6/Total IS负担violin图（耐药vs.敏感，三个子图） |
| `figures/fig_is_burden_combined.pdf/png` | 组合图（violin + ROC曲线 + IS26计数散点图） |

**预期关键结果（Q1–Q4）**:

| 问题 | 预期结果 |
|------|---------|
| Q1: IS26负担差异？ | 耐药株中位=7.5 vs 敏感株中位=0；Cliff's δ=0.614（全492队列） |
| Q2: 时间趋势？ | 2025队列(n=188, 15.4%耐药)中位=13 vs 0；AUC=0.976（最代表性） |
| Q3: 预测性能？ | 全队列AUC=0.807；2025队列AUC=0.976（提交偏差影响见Step 11） |
| Q4: 提交偏差？ | Pre-2022耐药率93.5% vs 2022+ 23.0%；Fisher OR=48.0, p=2×10⁻⁴⁶ |

**核心统计函数**:

```python
# Cliff's delta 效应量（非参数）
def cliffs_delta(a, b):
    mat = np.sign(a[:, None] - b[None, :])
    return float(mat.mean())
# δ=0.916 含义: 任取一对耐药/敏感基因组，95.8%的情况下耐药株IS26更多

# 从零实现逻辑回归 + ROC曲线（无sklearn依赖）
def logistic_regression_roc(X, y):
    # scipy.optimize.minimize → sigmoid → score → sort → trapezoid
    # 修复: np.trapz deprecated → getattr(np, 'trapezoid', getattr(np, 'trapz'))
    # 修复: AUC负值 → np.argsort(fprs) 保证单调排序

# Clopper-Pearson 精确置信区间
from scipy.stats import beta as beta_dist
ci_lo = beta_dist.ppf(0.025, k, n-k+1)
ci_hi = beta_dist.ppf(0.975, k+1, n-k)
```

**常见运行问题**:

```
# 错误: AttributeError: module 'numpy' has no attribute 'trapz'
# 原因: numpy 2.x 废弃 np.trapz，改名为 np.trapezoid
# 修正: trapz_fn = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)

# 错误: AUC = -0.952（负值）
# 原因: fprs数组降序排列，trapezoid积分方向错误
# 修正: order = np.argsort(fprs); fprs_s, tprs_s = fprs[order], tprs[order]

# 错误: fold_change格式化失败（NoneType）
# 原因: susceptible中位数=0，fold=None（不可定义）
# 修正: fold_str = f'{fold:.1f}×' if fold is not None and np.isfinite(fold) else '∞×'
```

---

## Step 9 — IS元件HMMER序列验证（L-005解决）✅

```bash
python analysis/09_is_hmmer_verify.py
```

**目的**: 验证GFF3注释偏差对IS26负担分析的影响（L-005）

**工作原理**:
1. 扫描全部270个基因组GFF3，提取所有IS元件注释
2. 从 `*_genomic.fna.gz` 中提取IS_unknown类型的蛋白质序列
3. 用PFAM PF01527（IS6家族HMM）通过pyhmmer `hmmsearch` 重新分类
4. 重新计算IS6负担统计，与注释级结果对比

**PFAM HMM自动下载**（如本地无缓存）:
```python
# 脚本自动从EBI下载，缓存至 data/db/pfam/PF01527.hmm
# EBI返回gzip压缩内容，脚本自动处理（无需手动解压）
```

**输出**:

| 文件 | 内容 |
|------|------|
| `data/processed/is_hmmer_results.tsv` | 21,760行IS特征HMMER分类结果 |
| `data/processed/is_burden_corrected.tsv` | 270行，每基因组含原始+校正IS6计数 |
| `data/processed/is_burden_corrected_stats.json` | 校正前后完整对比统计 |

**关键结果**:
```
IS_unknown查询总数:   4,148个蛋白质
PFAM IS6重分类:       1,002个 (826来自敏感株 / 176来自耐药株 → 4.7×偏差)
IS6 敏感株中位数:     0.0 → 5.0  (+5.0，注释偏差证实)
IS6 耐药株中位数:    13.0 → 16.0 (+3.0)
Cliff's δ:           0.903 → 0.899 (Δ=−0.004，可忽略)
AUC (可发表):        0.952 → 0.949 (Δ=−0.003，统计结论稳健)
```

**pyhmmer 0.12 API注意事项**（坑已踩过，不要改回旧写法）:
```python
# ✅ 正确写法
seq_name = hit.name           # str，不是bytes
query_name = hits.query.name  # 属性，不是hits.query_name

# ❌ 错误写法（AttributeError）
hit.name.decode()             # str没有.decode()
hits.query_name               # 不存在此属性
dom.alignment.identity        # 不存在此属性，用_domain_identity()函数代替
```

---

## Step 10 — AMR基因CARD序列验证（L-001解决）✅

```bash
python analysis/10_amr_hmmer_verify.py
```

**目的**: 验证所有碳青霉烯AMR基因检测的准确性（L-001）

**工作原理**:
1. 读取 `amr_hits.tsv` 中50个碳青霉烯AMR hits
2. 从对应 `*_genomic.fna.gz` 中提取蛋白质序列（流式读取，内存高效）
3. 用pyhmmer `phmmer` 搜索CARD v3.3蛋白质同源模型库
4. 按E-value最佳命中分类: CONFIRMED / NAME_MISMATCH / NO_HIT / EXTRACT_FAIL

**蛋白质提取逻辑**:
```python
# 1. 从GFF3坐标定位contig
# 2. 流式读取fna.gz，只加载目标contig（内存友好）
# 3. 切片 [start:stop]，如strand='-'则反向互补
# 4. 按标准密码子表翻译，跳过<50aa序列
```

**输出**:

| 文件 | 内容 |
|------|------|
| `data/processed/amr_hmmer_results.tsv` | 50行，含classification/card_gene/card_evalue |

**关键结果**:
```
总验证数:   217个碳青霉烯AMR hits（P1）
CONFIRMED:  217 (100%) — 所有命中均有CARD序列支持
NO_HIT:     0
NAME_MISMATCH: 0
EXTRACT_FAIL:  0

KPC-2确认: 293aa，N端MSLYRRLVLLSCLSW（经典KPC信号肽）
检测层级升级: GFF_KEYWORD → HMMER_VERIFIED
```

**常见运行问题**:
```
# 错误: KeyError: 'evalue' in best comparison
# 原因: 字典键是 'card_evalue' 不是 'evalue'
# 修正: if best is None or evalue < best['card_evalue']:

# 错误: AttributeError: 'str' object has no attribute 'decode'
# 原因: pyhmmer 0.12 hit.name已是str
# 修正: 删除所有 .decode() 调用

# 错误: card_meta lookup失败（bytes vs str键）
# 原因: card_meta字典键必须用str（header[1:].split()[0]），不能.encode()
```

---

## Step 11 — 时间分层与敏感性分析 ★P1新增

```bash
python analysis/11_subgroup_analysis.py
```

**依赖**: `data/processed/is_burden_all.tsv`（Step 08输出）, `is_context.tsv`（Step 04输出）, `amr_hits.tsv`（Step 03输出）

**分析模块**:

| 分析 | 内容 |
|------|------|
| A. 年份分层 | 各年份耐药率 + IS6负担 + AUC |
| B. 队列对比 | 全队列/Pre-2022/2022+/2025/2025-2026 的完整统计 |
| C. 复合转座子率 | 跨时代稳定性验证（主要科学发现不受偏差影响） |
| D. AMR基因分布 | 各时代KPC/NDM/IMP比例 |
| E. 2025队列详解 | 代表性队列主要分析（论文核心数字来源） |
| F. 提交偏差量化 | Fisher精确检验 OR=48.0, p=2×10⁻⁴⁶ |

**核心发现**:

```
NCBI提交偏差:
  Pre-2022:  115/123 = 93.5% resistant
  2022+:     85/369  = 23.0% resistant
  Fisher OR = 48.03,  p = 1.96×10⁻⁴⁶

2025代表性队列 (n=188):
  耐药率:     15.4% (29/188)
  IS6_res中位数: 13 (IQR 11–17)
  IS6_sus中位数: 0  (IQR 0–0)
  Cliff's δ:  0.952
  AUC:        0.976
  p:          2.2×10⁻²³

复合转座子率（跨时代稳健）:
  Pre-2022:   95.2%
  2022+:      94.2%
  2025:       88.0%
  全队列:     94.6%
```

**输出**:
- `data/processed/subgroup_stats.json` — 全部分层结果
- `data/processed/subgroup_table.tsv` — 年份统计表

---

## 一键重跑管道

```bash
cd e:/miniconda3/envs/llama-env/amr_project

# P1全流程（跳过下载，从审计开始）
python analysis/00_data_audit.py && \
python analysis/02_validate.py && \
python analysis/03_amr_scan.py && \
python analysis/04_is_context.py && \
python analysis/05_stats.py && \
python analysis/06_figures.py && \
python analysis/07_manuscript.py && \
python analysis/08_is_burden_all.py && \
python analysis/09_is_hmmer_verify.py && \
python analysis/10_amr_hmmer_verify.py && \
python analysis/11_subgroup_analysis.py

# 仅更新统计+图表+草稿（AMR扫描结果不变时最快）
python analysis/05_stats.py && python analysis/06_figures.py && python analysis/07_manuscript.py

# 仅重跑IS负担分析+时间分层（不涉及序列验证）
python analysis/08_is_burden_all.py && python analysis/11_subgroup_analysis.py

# 仅重跑序列验证（09/10独立，依赖08输出和CARD数据库）
python analysis/09_is_hmmer_verify.py && python analysis/10_amr_hmmer_verify.py
```

---

## 扩展到更大数据集

```python
# P1 已完成：MAX_GENOMES_CHINA = None，677个UIDs全量获取
# 未来P2建议：添加其他国家或物种

# config.py 修改示例
MAX_GENOMES_CHINA = None   # 当前P1设定

# 注意事项：
# - 01_download.py 已修复分页bug（None → retmax=10_000）
# - 00_data_audit.py 支持 --fast 模式（下载进行时的快速检查）
# - 181个GCF_前缀基因组仍有MISSING_GFF问题，可用datasets CLI补充下载
```

---

## 添加新物种（如 E. coli）

在 `config.py` 中 `TARGET_SPECIES` 添加：
```python
'Escherichia coli': {
    'genome_size_min_bp': 4_200_000,
    'genome_size_max_bp': 5_800_000,
    'min_cds': 3_000,
    'min_n50': 50_000,
}
```
在 `ASSEMBLY_SUMMARY_URLS` 中已有E. coli条目，`01_download.py` 中修改物种列表即可。

---

## 日志文件

| 文件 | 内容 |
|------|------|
| `logs/audit.log` / `audit_partial.log` | 数据完整性审计（Step 00，含每基因组检查结果） |
| `logs/download.log` | 每个基因组下载状态+MD5验证结果 |
| `logs/validate.log` | PASS/FAIL + 失败原因 |
| `logs/amr_scan.log` | 每基因组AMR hits数 + tier信息 |
| `logs/is_context.log` | 每基因组IS特征数 + composite数 |
| `logs/stats.log` | 全部统计结果文字报告 |
| `logs/subgroup.log` | 时间分层分析详细输出（Step 11，含bias量化） |

---

## 常见错误排查

| 错误信息 | 原因 | 解决 |
|---------|------|------|
| `KeyError: 'accession'` | download_status.tsv用`assembly_accession`列名 | 02_validate.py已自动重命名，检查是否旧版本 |
| `KeyError: 'species'` | download_status.tsv无species列 | 已修复：从organism_name派生 |
| `GFF parse error: not a gzip file` | 文件损坏或下载不完整 | 删除该目录，重新运行01_download.py |
| `ValidationError: N50 too low` | 组装质量差 | 正常，该基因组被排除 |
| `ModuleNotFoundError: scipy` | 环境未激活 | `conda activate llama-env` |
| 图表中文乱码 | matplotlib字体 | 使用英文标签（已实现） |
| `ModuleNotFoundError: pyhmmer` | pyhmmer未在当前环境安装 | `python -m pip install pyhmmer` |
| `AttributeError: 'str' has no attr decode` | pyhmmer 0.12 hit.name为str非bytes | 删除所有`.decode()`调用（见D-010） |
| `AttributeError: 'TopHits' has no attr query_name` | pyhmmer 0.12 API变更 | 改用`hits.query.name`（见D-010） |
| `AttributeError: 'Alignment' has no attr identity` | pyhmmer 0.12 Alignment无identity属性 | 用`_domain_identity()`手动计算（见D-010） |
| PFAM下载失败：非HMMER3格式 | EBI返回gzip压缩内容 | 检测`\x1f\x8b`魔术字节后先`gzip.decompress()`（见D-011） |
| CARD数据库不存在 | 需手动下载（非自动） | 从CARD官网下载broadstreet包，解压后放入`data/db/card/` |
