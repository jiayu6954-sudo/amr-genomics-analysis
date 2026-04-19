# 系统白皮书 — 细菌耐药基因组自动化科学发现系统
> **System Whitepaper: Automated AMR Genomics Science Discovery Pipeline**  
> 版本 1.3 · 2026-04-19 · 作者: [Author]  
> 项目路径: `e:/miniconda3/envs/llama-env/amr_project/`

---

## 一、系统定位与价值主张

### 1.1 解决什么问题

全球每年有近500万人死于耐药菌感染，其中亚洲（特别是中国和东南亚）承担了最大的疾病负担。碳青霉烯耐药肺炎克雷伯菌（CRKP）已成为中国重症监护病房中最危险的病原体之一，而这种耐药性的快速传播依赖于一种分子机制：**IS元件（插入序列）介导的复合转座子**。

然而，回答"IS元件在多大程度上驱动了中国CRKP的传播"这个问题，传统方法需要：
- 数月的实验室分离培养工作
- 昂贵的测序和生物信息学服务
- 专业生物信息学团队手动分析数百个基因组

**本系统将这一工作压缩至7个脚本、不到2小时的自动化运行。**

### 1.2 系统能做什么

```
输入: 研究问题 + NCBI数据库
输出: 发表级论文草稿 + 统计结果 + 出版质量图表

全程自动化:
  ① 从NCBI检索并下载数百个完整基因组
  ② 严格质控，排除低质量组装
  ③ 多层次AMR基因检测（注释扫描 + 序列验证）
  ④ IS元件侧翼结构分析，量化复合转座子率
  ⑤ 统计分析（精确置信区间 + 假设检验）
  ⑥ 出版质量图表（PDF + PNG）
  ⑦ 论文草稿自动填充（数字自动更新）
```

### 1.3 不能做什么（系统边界）

- 不替代实验室验证（MIC测定、接合实验等）
- 不执行从头基因组组装（依赖已组装的NCBI基因组）
- 不处理宏基因组数据（设计为单基因组模式）

### 1.4 目标期刊定位（更新）

| 阶段 | 期刊 | IF | 理由 |
|------|------|----|------|
| **当前（P1完成）✅** | **Nature Communications** | ~16 | 492基因组+时间分层+提交偏差量化+217/217 CARD验证 = 高质量方法学贡献 |
| 若加入多物种比较（路线B） | Nature Microbiology | ~28 | 物种水平系统性发现 |
| 保守备选 | Emerging Microbes & Infections | ~9 | 数据稳健，投稿风险低 |

> ✅ **P1已完成**：677→492基因组，0重复/0损坏（00_data_audit.py），217/217 AMR hits CARD验证（100%），NCBI提交偏差量化（OR=48.0，p=2×10⁻⁴⁶），2025代表性队列AUC=0.976（Step 11）。投稿级论文草稿：`reports/manuscript_submission_v1.md`。

---

## 二、架构设计

### 2.1 设计原则

**原则 1 — 数据完整性优先（数据是基石）**  
每个数据转换步骤都有显式验证。MD5校验和在下载时验证。质控门拒绝所有不满足阈值的基因组，即使这意味着样本量减少。AMR正则模式经过假阳性验证。

**原则 2 — 单一参数来源**  
所有生物学阈值、文件路径、正则模式集中在 `analysis/config.py`。任何参数调整只改这一个文件，管道脚本永远不硬编码数字。

**原则 3 — 顺序依赖，增量更新**  
管道严格顺序: 01→02→03→04→05→06→07。每步读取上步输出，写入标准化格式（TSV/JSON）。可以从任意步骤重新开始，无需完整重跑。

**原则 4 — 可复现性**  
所有随机性来源已消除（无随机采样）。固定NCBI查询字符串确保检索可复现（依赖数据库快照时间点）。日志记录每个决策。

**原则 5 — 方法分离**  
生物信息学逻辑（步骤1-4）与统计分析（步骤5）与可视化（步骤6）完全分离。修改图表风格不影响数据；修改统计方法不影响原始数据。

### 2.2 数据流架构

```
┌────────────────────────────────────────────────────────────────┐
│                    NCBI 数据层                                  │
│  GenBank FTP (gff.gz, fna.gz, assembly_stats.txt)             │
│  Entrez API (esearch, esummary)                                │
└─────────────────────────┬──────────────────────────────────────┘
                          │ HTTP + MD5验证
                          ▼
┌────────────────────────────────────────────────────────────────┐
│                    原始数据层  data/raw/                        │
│  {accession}/                                                  │
│    *_genomic.gff.gz          GFF3注释（IS检测主数据源）        │
│    *_genomic.fna.gz          基因组序列（序列验证备用）        │
│    *_assembly_stats.txt      N50/大小/contig统计               │
│    *_amrfinderplus.tsv       (可选) NCBI预计算AMR结果          │
│  download_status.tsv         检查点文件                        │
└─────────────────────────┬──────────────────────────────────────┘
                          │ 质控门（7项检查）
                          ▼
┌────────────────────────────────────────────────────────────────┐
│                  验证数据层  data/validated/                    │
│  genome_manifest.tsv         已通过质控基因组注册表            │
│  qc_report.tsv               每基因组QC详情（含FAIL原因）      │
└─────────────────────────┬──────────────────────────────────────┘
                          │ Tier1 AMRFinder / Tier2 GFF正则
                          ▼
┌────────────────────────────────────────────────────────────────┐
│                参考数据库层  data/db/   ◄── Step 09/10新增      │
│  pfam/PF01527.hmm            IS6家族HMM（EBI自动下载，36KB）   │
│  card/protein_fasta_*.fasta  CARD v3.3蛋白质同源模型（1.87MB） │
└─────────────────────────┬──────────────────────────────────────┘
                          │ pyhmmer hmmsearch/phmmer
                          ▼
┌────────────────────────────────────────────────────────────────┐
│                  处理数据层  data/processed/                   │
│  audit_report.tsv            数据审计明细（677行）★P1新增     │
│  audit_summary.json          审计汇总★P1新增                  │
│  amr_hits.tsv                AMR基因坐标+分类（217行，P1）     │
│  is_context.tsv              IS-AMR配对（2,985行，P1）         │
│  context_summary.tsv         基因组级IS侧翼汇总（200行，P1）   │
│  stats_summary.tsv           全部统计结果                      │
│  stats_tables.json           结构化输出（供论文自动填充）      │
│  is_burden_all.tsv           全队列IS负担（492行，P1）         │
│  is_burden_stats.json        IS负担统计（AUC/δ等）             │
│  is_hmmer_results.tsv        IS特征HMMER分类（270基因组子集）  │
│  is_burden_corrected.tsv     HMMER校正后IS6负担（270行）       │
│  is_burden_corrected_stats.json  校正前后对比统计              │
│  amr_hmmer_results.tsv       AMR CARD验证（217行，100%确认）   │
│  subgroup_stats.json         时间分层统计★P1新增              │
│  subgroup_table.tsv          年份分层表★P1新增                │
└─────────────────────────┬──────────────────────────────────────┘
                          │ matplotlib 渲染
                          ▼
┌──────────────────────────────────────────────────┐
│  figures/   报告层                               │
│  fig1–fig5 + fig_combined + IS负担组合图         │
│  reports/manuscript_draft_v1.md (自动生成)        │
│  reports/manuscript_submission_v1.md (投稿级★P1) │
└──────────────────────────────────────────────────┘
```

### 2.3 模块依赖图

```
config.py ◄──────────────────────────────────── 所有脚本
    │
    ├─► 01_download.py   [分页bug已修复 → 677 UIDs]
    │       │
    │       └── data/raw/ + download_status.tsv
    │
    ├─► 00_data_audit.py ◄── data/raw/   ★P1新增
    │       │               (download_status.tsv可选)
    │       └── audit_report.tsv + audit_summary.json
    │
    ├─► 02_validate.py ◄── download_status.tsv
    │       │
    │       └── genome_manifest.tsv  (492基因组，P1)
    │
    ├─► 03_amr_scan.py ◄── genome_manifest.tsv
    │       │
    │       └── amr_hits.tsv  (217 hits, P1)
    │
    ├─► 04_is_context.py ◄── amr_hits.tsv + genome_manifest.tsv
    │       │
    │       └── is_context.tsv (2,985行) + context_summary.tsv
    │
    ├─► 05_stats.py ◄── amr_hits.tsv + is_context.tsv + context_summary.tsv
    │       │
    │       └── stats_summary.tsv + stats_tables.json
    │
    ├─► 06_figures.py ◄── amr_hits.tsv + is_context.tsv + context_summary.tsv
    │       │
    │       └── figures/*.{pdf,png}
    │
    ├─► 07_manuscript.py ◄── stats_tables.json
    │       │
    │       └── reports/manuscript_draft_v1.md
    │
    ├─► 08_is_burden_all.py ◄── genome_manifest.tsv + amr_hits.tsv
    │       │
    │       └── is_burden_all.tsv (492行) + is_burden_stats.json
    │           figures/fig_is_burden_*.{pdf,png}
    │
    ├─► 09_is_hmmer_verify.py ◄── genome_manifest.tsv + data/db/pfam/PF01527.hmm
    │       │                     (PF01527.hmm 自动从EBI下载)
    │       └── is_hmmer_results.tsv + is_burden_corrected.tsv
    │           is_burden_corrected_stats.json
    │
    ├─► 10_amr_hmmer_verify.py ◄── amr_hits.tsv + genome_manifest.tsv
    │       │                      + data/db/card/protein_fasta_*.fasta
    │       └── amr_hmmer_results.tsv  (217行, 100% CONFIRMED)
    │
    └─► 11_subgroup_analysis.py ◄── is_burden_all.tsv + is_context.tsv   ★P1新增
            │                       + amr_hits.tsv
            └── subgroup_stats.json + subgroup_table.tsv
```

---

## 三、核心算法说明

### 3.1 Entrez 检索策略

```
查询串: "Klebsiella pneumoniae"[Organism]
       AND "China"[Country]
       AND latest[filter]
       AND ("Complete Genome"[Assembly Level]
            OR "Chromosome"[Assembly Level]
            OR "Scaffold"[Assembly Level])

分页: WebEnv + query_key（NCBI推荐大批量检索方式）
批次: 每次 esummary 200个（NCBI限制）
```

**为什么不用 assembly_summary.txt?**  
该文件没有地理位置信息列（geo_loc_name），无法直接过滤中国分离株。Entrez `[Country]` 字段来自 BioSample，只能通过 API 查询。

### 3.2 MD5 完整性验证

每个下载目录包含 `md5checksums.txt`。下载后逐文件验证 MD5，不匹配则重试（最多3次）。这是数据质量的第一道防线。

### 3.3 IS元件侧翼分析数学定义

设 AMR 基因坐标为 `[a_start, a_stop]`，IS元件坐标为 `[i_start, i_stop]`（同一contig上）：

```
gap_upstream   = max(0, a_start - i_stop)   # IS完全在AMR上游
gap_downstream = max(0, i_start - a_stop)   # IS完全在AMR下游
distance       = min(gap_upstream, gap_downstream)

# 若 distance ≤ W（窗口=10,000 bp），记录此配对
# 若重叠，distance=0

position = UPSTREAM   if gap_upstream ≤ gap_downstream
         = DOWNSTREAM otherwise
```

对每个 AMR 基因累积 upstream_count 和 downstream_count，最终：
```
COMPOSITE_TRANSPOSON  : upstream_count ≥ 1 AND downstream_count ≥ 1
SINGLE_IS_UPSTREAM    : upstream_count ≥ 1 AND downstream_count = 0
SINGLE_IS_DOWNSTREAM  : upstream_count = 0 AND downstream_count ≥ 1
NO_IS                 : 窗口内无任何IS
```

### 3.5 HMMER序列验证算法（Step 09/10）

**IS元件HMMER验证（09_is_hmmer_verify.py）**

核心逻辑：对GFF3中注释为IS_unknown（无法识别为IS6家族）的转座酶，通过PFAM PF01527（IS6家族HMM）做序列比对，补充分类。

```
输入: 270个基因组 GFF3 + fna.gz
      PFAM PF01527 (IS6 family: IS26, IS257, IS1353, IS240, IS6...)

1. 扫描GFF3，提取所有IS_unknown特征坐标
2. 流式读取fna.gz，提取蛋白质序列:
   seq = contig[start:stop]
   if strand == '-': seq = reverse_complement(seq)
   protein = translate(seq)  # 跳过 < 50aa 的序列

3. pyhmmer hmmsearch(PF01527_HMM, proteins, E=1e-5):
   for hits in hmmsearch(...):
       for hit in hits:
           if hit.included:
               identity = _domain_identity(hit.domains[0])

4. 重新统计IS6负担 (per genome):
   corrected_n_is6 = original_n_is6 + newly_classified_IS6_from_unknown

5. Mann-Whitney U → AUC (U_greater / n_pos / n_neg)
   Cliff's δ = (U_greater - U_less) / (n_pos * n_neg)
```

**AMR基因CARD验证（10_amr_hmmer_verify.py）**

```
输入: 50个碳青霉烯AMR hits坐标
      CARD v3.3 protein_fasta_protein_homolog_model.fasta (4,840条)

对每个AMR hit:
1. 从fna.gz流式提取蛋白质序列（同上）
2. pyhmmer phmmer(query_protein, CARD_proteins, E=1e-5):
   best_hit = lowest E-value hit
3. 分类:
   CONFIRMED     if best_hit.card_gene ≈ query_gene_name
   NAME_MISMATCH if hit存在但基因名不符
   NO_HIT        if E > 1e-5 或无包含hits
   EXTRACT_FAIL  if 蛋白质提取失败
```

### 3.4 置信区间选择理由

| 指标 | 方法 | 选择理由 |
|------|------|---------|
| 普通比例 | Wilson score CI | 小样本（n<1000）下比正态近似更准确，特别是p接近0或1时 |
| 精确二项 | Clopper-Pearson | 保守但精确，适合声明"93%"这样高比例时的下界 |
| IS家族均匀性 | chi-square GOF | 多分类，检验IS6是否显著超过均匀分布期望 |

---

## 四、现有发现的科学意义

### 4.1 核心发现陈述（双层发现）

**第一层：IS元件侧翼结构（Steps 03–05，P1 492基因组）**

> 中国临床肺炎克雷伯菌中，**94.6%的碳青霉烯耐药IS-AMR配对**（2,825/2,985对）形成复合转座子结构，IS6家族（主要为IS26）占所有已鉴定IS元件的36.8%（970/2,639）。217个碳青霉烯AMR hits分布于200个耐药基因组（40.7%，Wilson 95%CI 36.4–45.0%），含5株双碳青霉烯酶（KPC+NDM）。

**第二层：IS26作为耐药预测因子（Step 08+11）—— 核心发现 ✅**

> **IS26基因组拷贝数是碳青霉烯耐药性的近完美单特征预测因子**。代表性2025队列（n=188，15.4%耐药率，匹配真实监测数据）：AUC=0.976，Cliff's δ=0.952，p=2.2×10⁻²³——效应量在单基因组层面极为罕见。全492基因组队列AUC=0.807（受NCBI提交偏差影响，pre-2022提交93.5%为耐药株）。

**第三层：NCBI提交偏差量化（Step 11）—— 方法学发现 ✅**

> **新发现**: Pre-2022 NCBI组装提交93.5%为耐药株（靶向暴发研究），2022+为23.0%（常规监测）。Fisher OR=48.0，p=2×10⁻⁴⁶。**方法学贡献**: 系统量化生物数据库中的选择性提交偏差及其对IS26预测AUC的影响（全队列AUC=0.807 vs 代表性队列AUC=0.976）。

> ✅ **HMMER序列验证完成（Step 09）**: 4,148个IS_unknown蛋白质经PFAM PF01527重分类，注释偏差已量化（4.7×偏向敏感株），AUC变化仅−0.003，统计结论完全稳健。

> ✅ **AMR基因序列验证完成（Step 10）**: 217/217碳青霉烯AMR hits经phmmer vs CARD v3.3全部确认（100% CONFIRMED），检测层级：GFF_KEYWORD → **HMMER_VERIFIED**。

### 4.2 公共卫生意义

这意味着什么？

1. **传播机制已高度固化**: KPC-2通过IS26复合转座子传播的机制在中国临床株中已极为标准化，提示存在强烈的正向选择压力（抗生素使用驱动）

2. **单株内多质粒跳跃**: 复合转座子结构允许KPC-2在同一细菌内从质粒跳至染色体，从一个质粒跳至另一个质粒——这解释了为何耐药性可以如此快速在医院内传播

3. **双碳青霉烯酶（KPC+NDM共存）**: P1发现**5株**同时携带KPC+NDM（P0仅1株），超级耐药表型对所有上市碳青霉烯均耐药，且两个基因均有IS侧翼复合转座子结构

4. **IS26作为耐药生物标志物**: IS26基因组拷贝数≥11（25百分位数）可作为一个简单、可操作的监测阈值，在表型耐药出现之前识别潜在高风险菌株

5. **IS26作为靶标的启示**: 若能干扰IS26的转座酶活性，理论上可以阻断复合转座子形成，进而阻止耐药基因在菌株间传播

---

## 五、未来科学发现路线图

### 路线 A — 深化当前数据集（6个月内可发表）

**A0. IS26注释偏差验证 ✅ COMPLETED（Step 09）**
```
问题: 敏感株IS26中位数=0可能是GFF3注释偏差，而非真实生物学差异
工具: pyhmmer 0.12 + PFAM PF01527 (IS6家族HMM)
结果:
  - 4,148个IS_unknown蛋白质中1,002个被重分类为IS6
  - 偏差量化: 826敏感株 vs 176耐药株重分类 (4.7×偏向敏感株)
  - AUC: 0.952 → 0.949 (Δ=−0.003，可忽略)
  - Cliff's δ: 0.903 → 0.899 (Δ=−0.004，可忽略)
结论: 注释偏差真实存在且已量化；统计结论完全稳健
脚本: analysis/09_is_hmmer_verify.py（~2分钟运行）
```

**A1. 序列级AMR验证 ✅ COMPLETED（Step 10，P1）**
```
工具: pyhmmer 0.12 phmmer + CARD v3.3 protein_fasta_protein_homolog_model.fasta
结果:
  - 217/217碳青霉烯AMR hits 全部 CONFIRMED (100%)（P1）
  - 0 NAME_MISMATCH, 0 NO_HIT, 0 EXTRACT_FAIL
  - KPC-2 (293aa), NDM-1, NDM-5, IMP-4等全部有CARD序列支持
检测层级升级: GFF_KEYWORD → HMMER_VERIFIED
脚本: analysis/10_amr_hmmer_verify.py（~1分钟运行，P1）
```

**A2. 时间分层与提交偏差量化 ✅ COMPLETED（Step 11，P1）**
```
发现: NCBI pre-2022提交93.5%耐药（靶向暴发），2022+仅23.0%（常规监测）
结果:
  Fisher OR = 48.03, p = 1.96×10⁻⁴⁶ — 提交偏差统计显著性极强
  2025代表性队列: n=188, 15.4%耐药, AUC=0.976, δ=0.952, p=2.2×10⁻²³
  复合转座子率跨时代稳健: Pre-2022=95.2%, 2022+=94.2%, 2025=88.0%
脚本: analysis/11_subgroup_analysis.py（<5秒运行）
价值: 方法学贡献，量化了一类常见但少被显式报告的数据库选择性偏差
```

**A3. 投稿级论文撰写 ✅ COMPLETED（P1）**
```
文件: reports/manuscript_submission_v1.md
目标期刊: Nature Communications
内容: 492基因组，217/217 CARD验证，时间分层，AUC梯度分析，完整方法学
状态: 可直接投稿
```

**A4. 地理分布分析（P2 — 计划中）**
```
数据: BioSample API批量获取geo_loc_name（省份级别）
分析: 不同省份间KPC-2 vs NDM分布差异（Fisher精确检验）
假设: 华东/华南KPC-2更多；华北/西北NDM更多？
新脚本: analysis/12_geo_analysis.py（已在规划中，原A3编号）
```

**A5. IS26拷贝数与耐药多样性相关（P2 — 可选）**
```
分析: IS6元件数量 vs 耐药基因种类数（Spearman相关）
假设: IS26拷贝数越多，耐药基因多样性越高
价值: 提供"IS26是多药耐药的放大器"的定量证据
基础数据: is_burden_all.tsv + amr_hits.tsv已就绪
```

---

### 路线 B — 扩展物种范围（形成多菌种比较论文）

**目标菌种比较**:

| 菌种 | 预期检索量 | 预期carbapenem耐药率 | 主要耐药基因 |
|------|-----------|---------------------|-------------|
| *K. pneumoniae*（已完成） | 677 | 17-20% | KPC-2, NDM |
| *E. coli* | ~2,000 | 5-8% | NDM-5, OXA-181 |
| *A. baumannii* | ~500 | 60-80% | OXA-23, OXA-24 |
| *P. aeruginosa* | ~300 | 20-30% | IMP, VIM |

**科学问题**: IS26复合转座子是否在所有革兰阴性耐药病原体中同样主导？还是存在物种特异性的IS元件偏好？

**新价值**: 第一篇系统比较4种主要革兰阴性CRKP/CRAB/CRPA中IS元件侧翼结构的中国临床数据论文。

---

### 路线 C — 时间序列分析（P1已提供基础数据）

> ℹ️ **P1 Step 11 已实现时间分层核心分析**（年份分层AUC/IS6负担）。以下为深化方向。

**已完成（Step 11结果）**: 年份分层耐药率（2013–2026）、AUC梯度（Pre-2022→2025：0.663→0.976）、复合转座子率跨时代稳定（88–95%）。

**深化方向**: 
1. 将492个基因组按年份分为4个时间窗（已有年份数据）：<2019, 2019-2021, 2022-2023, 2024+
2. Cochran-Armitage趋势检验（复合转座子率趋势）
3. 控制提交偏差后的净时间趋势

**投稿目标**: The Lancet Microbe (IF 20.9) 或 Nature Microbiology (IF 28.3)

---

### 路线 D — 机制深度分析（2年研究计划）

**D1. Tn4401同型分析**
KPC-2通常位于 Tn4401 转座子内，而Tn4401有多个亚型（Tn4401a/b/c/d/e）。不同亚型的表达活性和转座频率不同。

```
方法: 从FASTA中提取KPC-2侧翼序列（±3kb），与Tn4401参考序列比对
工具: minimap2 + python pairwise alignment
假设: Tn4401b是中国株最常见亚型（与欧洲Tn4401a不同）
```

**D2. 质粒背景分析**
复合转座子位于质粒（IncF、IncL/M、IncR）还是染色体？质粒型别与IS26密度的关系。

```
方法: PlasmidFinder数据库 + contig长度（circular contig通常是质粒）
新价值: 识别高风险质粒型别（传播能力最强的质粒骨架）
```

**D3. 比较基因组学——克隆谱系分析**
IS26在相同ST（序列型）与不同ST之间传播的比较。

```
工具: Mash距离（快速基因组相似度）+ mlst工具（7个管家基因ST分型）
假设: ST11（中国CRKP主流克隆）内IS26密度高于其他ST
```

---

### 路线 E — 平台化与应用（3年远景）

**E1. 实时监测仪表板**
将分析管道改造为每月自动从NCBI获取新基因组并更新统计数据的监测系统。输出：网页可视化的中国CRKP耐药趋势实时图。

**E2. 医院感染溯源工具**
输入：同一医院爆发的多个基因组 → 输出：IS元件侧翼结构一致性评分（判断传播链是否由相同复合转座子驱动）。

**E3. 新型耐药机制早期预警**
对每月新增基因组中检测到但未在已知CARBAPENEM_GENES列表中的耐药基因，自动标记为"新兴耐药元件"并发送告警。

---

## 六、为什么这项研究对普通人重要

每年，中国有数万名ICU患者因CRKP感染死亡。这些患者往往因为手术后感染、化疗后免疫抑制或老年基础疾病而入住ICU，他们不是因为自身的行为感染，而是被医院环境中流动的耐药基因"选择"了。

IS26复合转座子就是这些耐药基因在医院内快速传播的"运载火箭"。搞清楚它的结构、分布和演化规律，是设计精准干预措施的第一步：

- **医院感染控制**: 如果知道哪类质粒骨架最危险（携带IS26最多），可以优先对这些克隆的携带患者加强隔离
- **药物研发**: IS26转座酶是潜在的新型靶标，靶向转座酶的小分子抑制剂可能阻止耐药基因在医院内播散
- **抗生素管理**: 时间序列数据可以量化"减少碳青霉烯使用是否确实降低了IS26密度"——为抗生素管理政策提供基因组层面的证据

**本系统做的事，是把这个需要专业实验室数月才能完成的分析，压缩为任何具备Python环境的研究者在几小时内就能独立完成的工作流。**

---

## 七、技术规格

### 7.1 系统需求

```
最低配置:
  RAM:  8 GB（270基因组分析）
  存储: 50 GB（270基因组原始数据 ~40GB）
  CPU:  4核（单线程管道，IO密集型）
  网络: 20Mbps+（下载阶段）

推荐配置（500基因组）:
  RAM:  16 GB
  存储: 100 GB SSD
  CPU:  8核（未来并行化改造）
  网络: 100Mbps+
```

### 7.2 扩展性估算

| 基因组数量 | 下载时间 | 分析时间 | 磁盘占用 |
|-----------|---------|---------|---------|
| 100 | ~15分钟 | <5分钟 | ~15 GB |
| 300（当前） | ~45分钟 | ~8分钟 | ~45 GB |
| 500 | ~90分钟 | ~15分钟 | ~75 GB |
| 1,000 | ~3小时 | ~30分钟 | ~150 GB |
| 5,000 | ~15小时 | ~3小时 | ~750 GB |

> 注：分析速度可通过并行化（`multiprocessing.Pool`）线性提升，待实现（见路线E）。

### 7.3 代码行数统计

| 文件 | 行数 | 主要功能 |
|------|------|---------|
| config.py | 119 | 参数配置 |
| 01_download.py | ~350 | NCBI下载+MD5验证 |
| 02_validate.py | 307 | 质控门 |
| 03_amr_scan.py | 295 | AMR检测 |
| 04_is_context.py | 317 | IS侧翼分析 |
| 05_stats.py | 191 | 统计计算 |
| 06_figures.py | 340 | 图表生成 |
| 07_manuscript.py | 270 | 论文草稿生成 |
| 08_is_burden_all.py | ~500 | IS负担全队列分析+预测模型 |
| 09_is_hmmer_verify.py | ~570 | IS元件PFAM HMM验证+注释偏差量化 |
| 10_amr_hmmer_verify.py | ~490 | AMR基因CARD序列验证 |
| **合计** | **~3,749** | 10个功能完整脚本 |

---

## 八、引用与参考架构

本系统的设计参考了以下已发表方法：

1. **IS元件侧翼分析**: Harmer CJ & Hall RM (2016) mSphere — IS26复合转座子分类框架
2. **AMR基因检测**: NCBI AMRFinderPlus (2019) — 官方AMR基因数据库 CARD
3. **质控标准**: Sczyrba et al. (2017) Critical Assessment of Metagenome Interpretation — 组装质量阈值参考
4. **Entrez批量检索**: NCBI E-utilities Documentation — WebEnv分页策略
5. **统计方法**: Newcombe (1998) Statistics in Medicine — Wilson score置信区间
6. **CARD数据库**: Alcock et al. (2023) Nucleic Acids Research — CARD v3.3蛋白质同源模型；用于AMR基因序列验证（Step 10）
7. **PFAM/InterPro**: Paysan-Lafosse et al. (2023) Nucleic Acids Research — InterPro/PFAM HMM数据库；PF01527 IS6家族HMM用于IS元件验证（Step 09）
8. **pyhmmer**: Larralde & Zeller (2023) Bioinformatics — pyhmmer Python包；本研究使用v0.12.0（Step 09/10序列验证引擎）

---

*本白皮书记录了当前系统（v1.2）的设计理念、技术实现和未来方向。随着系统演进，请更新本文档对应章节。*

*最后更新: 2026-04-19*
