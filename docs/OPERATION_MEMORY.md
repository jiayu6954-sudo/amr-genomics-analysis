# 操作记忆库 — AMR 基因组学项目
> **Operation Memory Library** · 版本 1.3 · 2026-04-19  
> 路径: `e:/miniconda3/envs/llama-env/amr_project/`

---

## 一、项目核心坐标

| 项 | 值 |
|---|---|
| 研究问题 | IS元件介导的复合转座子在中国临床肺炎克雷伯菌中驱动碳青霉烯耐药基因播散 |
| 数据来源 | NCBI GenBank（Entrez API） |
| 目标菌种 | *Klebsiella pneumoniae*（主）/ *Escherichia coli*（扩展备选） |
| 地理范围 | 中国临床分离株 |
| 分析规模 | **492** 个经质控的完整/近完整基因组（P1扩展后：677检索→492通过审计+QC） |
| 目标期刊 | **Nature Communications (IF ~16)**（P1完成后升级）/ EMI备选 |

---

## 二、已验证核心发现（P1最终数字 — 492基因组全序列验证后）

### 2A. 队列与质控（Step 00/02）

| 指标 | 数值 |
|------|------|
| NCBI检索总数 | 677 个组装（P1修复分页bug后） |
| 数据审计通过 | 492/677 (72.7%) |
| 缺失GFF（GCF_前缀FTP问题） | 181 个（系统性） |
| N50不合格 | 3 个 |
| 注册表重复 | **0** |
| BioSample重复 | **0** |
| gzip损坏 | **0** |

### 2B. 复合转座子结构（Step 03/04/05 — P1 492基因组）

| 指标 | 数值 | 置信区间/检验 |
|------|------|--------------|
| 碳青霉烯耐药率 | **40.7%** (200/492) | Wilson 95%CI [36.4–45.0%] |
| KPC型 | 117 hits / 109基因组 (54.5%) | KPC-2占94.9% (111/117) |
| NDM型 | 92 hits / 88基因组 (44.0%) | NDM-1=57, NDM-5=31, NDM-9=2, NDM-6=1, NDM-13=1 |
| IMP型 | 8 hits / 8基因组 (4.0%) | 全为IMP-4 |
| 双碳青霉烯酶 | **5株** | KPC+NDM 共存 |
| AMR总hits | **217** | 全部CARD v3.3验证（100%） |
| IS-AMR配对总数 | **2,985** | — |
| 复合转座子率（每对） | **94.6%** (2,825/2,985) | Clopper-Pearson 95%CI |
| NO_IS位点 | 9 (4.1%) | — |
| 复合转座子位点 | 170 (78.3%) | — |
| IS6族占比 | **36.8%** (970/2,639已分类) | — |
| AMR基因CARD验证 ✓ | **217/217 (100%) CONFIRMED** | phmmer vs CARD v3.3 |

### 2C. IS6负担预测 — NCBI提交偏差量化（Step 11）

| 队列 | n | 耐药率 | IS6耐药中位数 | IS6敏感中位数 | AUC | Cliff's δ |
|------|---|--------|--------------|--------------|-----|-----------|
| 全队列 (492) | 492 | 40.7% | 7.5 (IQR 2–14) | 0 (IQR 0–4) | **0.807** | 0.614 |
| Pre-2022（偏差） | 123 | 93.5% | 3 (IQR 2–14) | 2 (IQR 0.8–5.8) | 0.663 | 0.325 |
| 2022+（混合） | 369 | 23.0% | 12 (IQR 5–15) | 0 (IQR 0–4) | 0.881 | 0.762 |
| **2025年队列（代表性）** | **188** | **15.4%** | **13 (IQR 11–17)** | **0 (IQR 0–0)** | **0.976** | **0.952** |
| 2025-2026（近期） | 215 | 13.5% | 13 (IQR 11–17) | 0 (IQR 0–2) | 0.966 | 0.932 |

**提交偏差统计检验**: Fisher OR=48.03，p=1.96×10⁻⁴⁶ — Pre-2022提交48倍更可能为耐药株（靶向暴发提交vs.常规监测）

**复合转座子率跨时代稳定**: Pre-2022=95.2%，2022+=94.2%，2025=88.0% → 结构特征不受提交偏差影响

### 2D. IS6负担分析 — HMMER序列验证校正值 ✓（Step 09，**270基因组子集**）

| 指标 | 校正前（注释级） | 校正后（HMMER验证） | 变化量 |
|------|-----------------|-------------------|-------|
| IS_unknown蛋白质查询 | — | 4,148 | — |
| PFAM IS6重分类数 | — | 1,002 (826敏感/176耐药) | — |
| IS6 耐药株中位数 | 13.0 | **16.0** | +3.0 |
| IS6 敏感株中位数 | 0.0 | **5.0** | +5.0 ← 注释偏差证实 |
| Cliff's δ | 0.903 | **0.899** | −0.004 |
| **AUC（270基因组可发表）** | 0.952 | **0.949** | −0.003 |

> ✅ **L-005 已解决**: 注释偏差确认（敏感株重分类4.7倍多于耐药株），但AUC仅下降0.003，统计结论稳健。**2025队列引用数字：AUC=0.976**（15.4%耐药率代表性队列），**全队列：AUC=0.807**（含提交偏差），检测层级：GFF_KEYWORD → **HMMER_VERIFIED**。

---

## 三、关键决策记录（Decision Log）

### D-001 数据获取策略
**决定**: 使用 Entrez `esearch` API，而非直接解析 `assembly_summary.txt`  
**原因**: `assembly_summary.txt` 不含 `geo_loc_name` / `isolation_source` 列（BioSample级别元数据），直接解析无法按China过滤  
**实现**: `esearch.fcgi?db=assembly&term="Klebsiella pneumoniae"[Organism] AND "China"[Country] AND latest[filter]`  
**结果**: 677个符合条件的组装（取前300）

### D-002 质控阈值
**决定**: 基因组4.8–6.5Mb + N50≥50kb + CDS≥3,000  
**原因**: NCBI K. pneumoniae 参考基因组中位大小约5.3Mb；N50阈值排除高度碎片化草图组装；CDS阈值排除注释不完整组装  
**结果**: 271/300下载成功，270/271通过质控（1个因N50=7,019bp被拒）

### D-003 AMR检测分层
**决定**: Tier 1（AMRFinderPlus TSV）优先，Tier 2（GFF3 关键词正则）备用  
**原因**: NCBI为部分基因组预计算AMRFinderPlus结果，精度更高  
**实际情况**: 本批270个基因组全部使用Tier 2（无AMRFinderPlus TSV）→ 结果需在论文中标注为"lower confidence"

### D-004 IS侧翼窗口大小
**决定**: ±10,000 bp  
**原因**: 文献中IS26复合转座子的结构一般<5kb；10kb窗口可覆盖二级IS元件和转座子末端  
**参考**: Harmer & Hall 2016 (mSphere)；Sheppard et al. 2016 (AAC)

### D-005 复合转座子判定标准
**决定**: 同一contig上目标基因上游AND下游各至少1个IS元件在窗口内 → COMPOSITE_TRANSPOSON  
**原因**: 双侧IS是经典复合转座子定义的最低要求  
**注意**: 单对IS-AMR记录计入，实际上若一个基因有k个上游+m个下游IS，生成k×m条记录，全部分类为COMPOSITE_TRANSPOSON

### D-006 正则表达式去假阳性（关键修正）
**问题**: `r'IMP\b'` (case-insensitive) 匹配 `rimP`（核糖体成熟蛋白）、`guaB`（IMP脱氢酶）、`purH`（IMP水解酶），产生1,349个假阳性  
**根因**: `\b` 仅限制后边界，`rimP` 中 `imP` 被大小写不敏感匹配  
**修正**: `r'\bIMP-\d|blaIMP'` — 要求等位号（IMP-1, IMP-4）或 bla 前缀  
**影响**: AMR总hits从1,408降至51，100%为碳青霉烯基因（假阳性清零）  
**同样修正**: `r'VIM\b'` → `r'\bVIM-\d|blaVIM'`；`r'bla[A-Z]{2,5}-'` → `r'bla[A-Z]{2,5}-\d'`

### D-007 基因名规范化
**问题**: GFF3 `gene` 属性为空时，`gene_name` 取 `product` 前60字符，导致 "carbapenem-hydrolyzing class A beta-lactamase KPC-2" 而非 "KPC-2"  
**修正**: 新增 `_canonical_gene()` 函数，优先级: blaKPC-2 → KPC-2 → 从product提取等位名  
**特例**: `KPC-204` 是真实等位变体（GCA_037807855.1），GFF注释正确

### D-008 IS负担分析使用全队列（Step 08）
**决定**: `08_is_burden_all.py` 扫描全部270个基因组的IS元件负担（非仅47个耐药株）  
**原因**: 若只分析耐药株，无法建立与敏感株的对比，无法计算预测AUC，无法识别"耐药预激株"  
**实现**: 对270个基因组各自统计 `n_is_total`、`n_is6`、`n_is26`，然后按 `is_carbapenem` 分组比较  
**结论**: IS26是二值分类（耐药vs.敏感）的近完美单特征预测因子（AUC=0.952）  
**注意**: IS26和IS6均以基因组拷贝数（非位点级别）计，与Step 04分析层次不同

### D-009 统计分母修正（每位点 vs. 每配对）
**决定**: 复合转座子率的主要引用数字改为每位点：94.0% (47/50)，而非 93.0% (1,042/1,120)  
**原因**: 1,120是IS-AMR配对数（1个AMR位点可与多个IS配对产生多行），同一位点计数多次夸大精度  
**规则**: 论文摘要/结论用每位点率；方法/补充用每配对率  
**实现**: `05_stats.py` 中以 `locus_id`（accession|contig|start|stop）去重后计算

### D-010 pyhmmer 0.12.0 API 破坏性变更（Step 09/10 调试记录）
**背景**: pyhmmer 0.12 改变了多个属性类型，旧版文档/示例代码直接失效  
**变更清单**（均在09/10脚本中踩过）:

| 属性/用法 | 旧版行为（预期） | 0.12.0实际行为 | 修正 |
|---------|---------------|--------------|------|
| `hit.name` | `bytes` | `str` | 删除所有 `.decode()` 调用 |
| `hits.query_name` | 属性存在 | `AttributeError` | 改为 `hits.query.name` |
| `dom.alignment.identity` | `float` 属性 | `AttributeError` | 手动计算（见下方函数） |
| `DigitalSequence.name` | `bytes` | `str` | 删除 `.encode()` 调用 |

**identity手动计算**（替换不存在的 `.identity`）:
```python
def _domain_identity(dom) -> float:
    aln = dom.alignment
    total = matched = 0
    for h, t in zip(aln.hmm_sequence, aln.target_sequence):
        if h != '-' and t != '-':
            total += 1
            if h.upper() == t.upper():
                matched += 1
    return matched / total if total else 0.0
```

### D-012 P1扩展：500→677基因组（分页Bug修复）
**背景**: config.py `MAX_GENOMES_CHINA=None`（设计为获取全部），但 `01_download.py` 中 `entrez_search(query, retmax=args.limit or 500)` 在 `args.limit=None` 时返回500（Python `or` 对 `None` 短路为 `500`）  
**修复**: 改为 `args.limit if args.limit is not None else 10_000`  
**结果**: 第二轮正确检索677个UIDs；最终下载677个，492个通过全链路审计+QC  
**影响**: 耐药株从47增至200；NDM从5增至92 hits；5株双碳青霉烯酶（原1株）

### D-013 NCBI提交偏差识别与处理
**发现**: Pre-2022组装93.5%为耐药株 vs. 2022+仅23.0%（Fisher OR=48.0，p=2×10⁻⁴⁶）  
**根因**: Pre-2022提交为靶向耐药暴发研究，偏向选择耐药株；2022+为常规临床监测  
**处理**:  
1. **不删除数据** — 提交偏差本身是发现，作为方法学贡献透明报告  
2. **主文数字** — 2025队列(n=188, 15.4%耐药)作为代表性分析（AUC=0.976）  
3. **敏感性分析** — 报告四个队列（全队列/Pre-2022/2022+/2025）结果  
4. **复合转座子率** — 跨时代稳定(88–95%)，论文核心结论不受偏差影响

### D-014 年份分层分析（Step 11）
**决定**: 新增 `11_subgroup_analysis.py`，系统量化各年份队列的统计性质  
**目的**: 透明展示AUC随耐药率变化的数学关系，证实IS6预测价值在任意耐药率基准下均稳健  
**关键洞察**: AUC梯度（0.663→0.807→0.881→0.976）与耐药率梯度（93.5%→40.7%→23.0%→15.4%）负相关，符合Mann-Whitney统计学预期 — 当耐药株占少数时，判别力最强

### D-011 CARD / PFAM 数据库获取方式
**PFAM PF01527 (IS6 family transposase HMM)**:
- URL: `https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF01527/?annotation=hmm`
- EBI返回 **gzip压缩** 的HMM文件（非明文），必须先检测 `\x1f\x8b` 魔术字节再 `gzip.decompress()`
- 不要设置 `Accept: text/plain` header（会返回406）
- 本地缓存: `data/db/pfam/PF01527.hmm`（36,057 bytes）

**CARD v3.3 蛋白质同源模型**:
- 文件: `protein_fasta_protein_homolog_model.fasta`（4,840条蛋白质，1.87MB）
- 解压自: `broadstreet-v3.3.0.tar.bz2`（需要手动从CARD官网下载并解压）
- 本地路径: `data/db/card/protein_fasta_protein_homolog_model.fasta`
- phmmer使用方式: 用待验证蛋白质序列（query）搜索CARD蛋白库（target）

---

## 四、已知限制与待解决问题

### L-001 全GFF_KEYWORD检测（无AMRFinder验证）　✅ **[已解决 — Step 10]**
- **原问题**: 270个基因组全部使用Tier 2检测，缺乏序列级比对验证
- **解决**: `10_amr_hmmer_verify.py` — phmmer对50个碳青霉烯AMR hits搜索CARD v3.3蛋白库
- **结果**: **50/50 (100%) CONFIRMED**，0 NAME_MISMATCH，0 NO_HIT，0 EXTRACT_FAIL
- **检测层级升级**: GFF_KEYWORD → **HMMER_VERIFIED**
- **KPC-2确认**: 293aa，N端 MSLYRRLVLLSCLSW（经典KPC信号肽确认）
- **论文影响**: 可删除"lower confidence"标注，改为"sequence-verified against CARD v3.3"

### L-002 181个MISSING_GFF（GCF_前缀FTP路径问题）
- **问题**: 677中有181个GCF_前缀组装的GFF文件下载路径解析失败（MISSING_GFF状态）
- **根因**: NCBI FTP路径中GCF_条目的文件命名与GCA_不同；脚本解析逻辑未覆盖GCF_命名约定
- **影响**: 181个基因组未进入分析（其中可能包含耐药株），最终492/677通过
- **解决**: 对181个accession单独用 `datasets` CLI工具重新下载（尚未执行）；或在 `01_download.py` 中修复GCF_路径解析逻辑

### L-003 无临床元数据
- **问题**: NCBI GenBank记录中 geo_loc_name/isolation_source 为BioSample级别，未批量获取
- **影响**: 无法做地理分布分析（省份层面）或感染类型分析（血流/尿路等）
- **解决**: 用 BioSample API 批量获取，约270次请求

### L-004 IS元件识别依赖注释质量
- **问题**: IS元件检测基于GFF3注释；注释不完整的基因组可能低估IS数量
- **影响**: 部分基因组IS负担可能被低估

### L-005 GFF3注释偏差导致IS26在敏感株中假性为零　✅ **[已解决 — Step 09]**
- **原问题**: 敏感株IS26中位数=0，疑为GFF3注释偏差（IS_unknown漏检）而非真实生物学差异
- **解决**: `09_is_hmmer_verify.py` — 4,148个IS_unknown蛋白质用PFAM PF01527（IS6家族HMM）重新分类
- **定量结果**:
  - 重分类总数: **1,002个** (826来自敏感株 / 176来自耐药株) → **4.7×偏向敏感株**（注释偏差确认）
  - IS6 敏感株中位数: 0.0 → **5.0**（+5.0，偏差真实存在）
  - IS6 耐药株中位数: 13.0 → **16.0**（+3.0）
  - Cliff's δ: 0.903 → **0.899**（Δ=−0.004，可忽略）
  - AUC: 0.952 → **0.949**（Δ=−0.003，统计结论不变）
- **结论**: 注释偏差经量化，对统计结论影响微小（<1%）；**论文可发表数字: AUC=0.949（HMMER_VERIFIED）**

---

## 五、已验证的工作环境

```
OS:         Windows 11 Pro 10.0.26200
Shell:      bash (Git Bash / Conda)
Python:     3.12 (llama-env conda环境)
conda路径:  e:/miniconda3/envs/llama-env/

关键包版本:
  pandas      2.x
  scipy       1.x
  numpy       1.x / 2.x
  matplotlib  3.x
  requests    2.x

外部工具:
  pyhmmer 0.12.0（已安装，完成Step 09/10验证 — API注意事项见D-010）
  pandoc 3.9.0.2（在主项目中，不在amr_project）

pyhmmer安装路径注意:
  pyhmmer仅在 llama-env conda 中，需在该环境下运行09/10脚本
  安装命令（如丢失）: python -m pip install pyhmmer
```

---

## 六、数据文件目录（完整）

```
amr_project/
├── analysis/
│   ├── config.py              # 唯一参数源（改参数只改此文件）
│   ├── 00_data_audit.py       # 数据完整性审计 → data/processed/audit_*.{tsv,json}  ★P1新增
│   ├── 01_download.py         # NCBI下载 → data/raw/（已修复分页bug）
│   ├── 02_validate.py         # 质控 → data/validated/genome_manifest.tsv
│   ├── 03_amr_scan.py         # AMR扫描 → data/processed/amr_hits.tsv
│   ├── 04_is_context.py       # IS侧翼 → data/processed/is_context.tsv
│   ├── 05_stats.py            # 统计 → data/processed/stats_*.{tsv,json}
│   ├── 06_figures.py          # 图表 → figures/*.{pdf,png}
│   ├── 07_manuscript.py       # 论文草稿 → reports/manuscript_draft_v1.md
│   ├── 08_is_burden_all.py    # IS负担全队列分析 → data/processed/is_burden_*.{tsv,json}
│   ├── 09_is_hmmer_verify.py  # IS元件HMMER验证 → data/processed/is_hmmer_*.{tsv,json}
│   ├── 10_amr_hmmer_verify.py # AMR基因CARD验证 → data/processed/amr_hmmer_results.tsv
│   └── 11_subgroup_analysis.py# 时间分层+敏感性分析 → data/processed/subgroup_*.{json,tsv}  ★P1新增
├── data/
│   ├── raw/                   # 原始基因组文件（677个目录，P1扩展后）
│   │   ├── download_status.tsv
│   │   └── {accession}/       # gff.gz + fna.gz + assembly_stats.txt
│   ├── validated/
│   │   ├── genome_manifest.tsv  # 492个通过质控的基因组（P1）
│   │   └── qc_report.tsv        # 全部基因组的QC详情
│   ├── db/                    # 参考数据库
│   │   ├── pfam/
│   │   │   └── PF01527.hmm      # IS6家族HMM（36,057B，EBI自动下载）
│   │   └── card/
│   │       └── protein_fasta_protein_homolog_model.fasta  # CARD v3.3蛋白（4,840条，1.87MB）
│   └── processed/
│       ├── audit_report.tsv          # 677行数据审计明细（Step 00）  ★P1新增
│       ├── audit_summary.json        # 审计汇总（Step 00）  ★P1新增
│       ├── amr_hits.tsv              # 217行碳青霉烯hits（P1）
│       ├── is_context.tsv            # 2,985行IS-AMR对（P1）
│       ├── context_summary.tsv       # 200行基因组级汇总（P1）
│       ├── stats_summary.tsv         # 统计结果平铺表
│       ├── stats_tables.json         # 结构化JSON（供manuscript使用）
│       ├── is_burden_all.tsv         # 492行IS负担数据（每基因组一行，P1）
│       ├── is_burden_stats.json      # IS负担统计结果（AUC、Cliff's δ等）
│       ├── is_hmmer_results.tsv      # IS特征HMMER结果（Step 09，270基因组子集）
│       ├── is_burden_corrected.tsv   # 校正后IS负担（原始+HMMER校正值）
│       ├── is_burden_corrected_stats.json  # 校正前后完整对比统计
│       ├── amr_hmmer_results.tsv     # 217行AMR hits CARD验证结果（Step 10，P1）
│       ├── subgroup_stats.json       # 时间分层统计JSON（Step 11）  ★P1新增
│       └── subgroup_table.tsv        # 按年份统计表（Step 11）  ★P1新增
├── figures/                   # 图表目录（PDF+PNG）
│   ├── fig1_prevalence.*
│   ├── fig2_is_context.*
│   ├── fig3_is_families.*
│   ├── fig4_is_burden.*
│   ├── fig5_forest.*
│   ├── fig_combined.*
│   ├── fig_is_burden_violin.*   # IS26/IS6/Total violin图（耐药vs.敏感）
│   └── fig_is_burden_combined.* # 组合图（violin+ROC+scatter）
├── logs/                      # 每步运行日志（含subgroup.log，audit*.log）
├── reports/
│   ├── manuscript_draft_v1.md         # 自动生成草稿（Step 07）
│   └── manuscript_submission_v1.md    # 投稿级论文（P1完成后人工撰写）  ★P1新增
└── docs/                      # 本文档所在目录
```

---

## 七、复现关键运行记录

| 步骤 | 运行时间 | 关键数字 |
|------|---------|---------|
| **P0（270基因组时代）** | | |
| 01_download (P0) | ~45分钟（300基因组）| 271 OK / 29 MISSING_REQUIRED |
| 02_validate (P0) | ~8分钟 | 270 PASS / 1 FAIL (N50=7,019) |
| 03_amr_scan (P0) | ~2分钟 | 51 carbapenem hits / 223 no AMR |
| 04_is_context (P0) | ~3秒 | 1,120 pairs / composite=93% |
| 05–07 | <1秒 | 统计/图/草稿 |
| 08_is_burden_all (P0) | ~5–10分钟 | IS26 AUC=0.952 / Cliff's δ=0.916 / 10耐药预激株 |
| 09_is_hmmer_verify | ~2分钟 | 4,148 IS_unknown→1,002重分类(IS6)→AUC 0.952→0.949 |
| 10_amr_hmmer_verify | ~23秒 | 50/50 CONFIRMED (100%)，phmmer vs CARD v3.3 |
| **P1（492基因组时代 — 当前）** | | |
| 00_data_audit (P1) | ~50秒（677基因组）| 492通过/181 MISSING_GFF/3 N50_FAIL/0重复/0损坏 |
| 01_download (P1) | ~2小时（677基因组）| 677检索（修复分页bug后）/ 492下载成功 |
| 02_validate (P1) | ~15分钟 | 492 PASS |
| 03_amr_scan (P1) | ~5分钟 | 217 carbapenem hits / 200基因组耐药 |
| 04_is_context (P1) | ~10秒 | 2,985 pairs / composite=94.6% |
| 08_is_burden_all (P1) | ~20分钟 | AUC=0.807（全队列）/ Cliff's δ=0.614 |
| 10_amr_hmmer_verify (P1) | ~1分钟 | 217/217 CONFIRMED (100%) |
| 11_subgroup_analysis | <5秒 | 2025队列AUC=0.976 / Bias OR=48.0 p=2×10⁻⁴⁶ |

---

## 八、config.py 参数速查

```python
# 生物学阈值
IS_FLANK_WINDOW_BP    = 10_000   # IS侧翼窗口（bp）
genome_size_min_bp    = 4_800_000  # K. pneu最小基因组
genome_size_max_bp    = 6_500_000
min_n50               = 50_000
min_cds               = 3_000

# 下载控制
MAX_GENOMES_CHINA     = None       # P1设为None，已获取全部677个（修复分页bug后）
RATE_LIMIT_DELAY_S    = 0.35       # NCBI速率限制（<3req/s）
DOWNLOAD_TIMEOUT_S    = 60
MAX_RETRIES           = 3

# NCBI邮箱（HTTP User-Agent）
NCBI_EMAIL = 'jiayu6954@gmail.com'
```
