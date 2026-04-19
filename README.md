# AMR Genomics Analysis — IS26-Mediated Composite Transposons in Chinese Clinical *K. pneumoniae*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
[![CARD v3.3](https://img.shields.io/badge/CARD-v3.3-green.svg)](https://card.mcmaster.ca/)
[![Genome n=492](https://img.shields.io/badge/genomes-n%3D492-orange.svg)]()

> **Population genomic analysis of IS element architecture around carbapenem resistance genes in 492 Chinese clinical *Klebsiella pneumoniae* genomes.**

---

## Overview

This repository contains the complete automated analysis pipeline for:

> **ZJY. (2026).** *IS26-Mediated Composite Transposons Underlie the Near-Universal Flanking Architecture of Carbapenem Resistance Genes in Chinese Clinical Klebsiella pneumoniae: A Population Genomic Analysis.* (Submitted to *Nature Communications*)

### Key Findings

| Finding | Result |
|---------|--------|
| Cohort (after 7-criterion QC) | **492 / 677** genomes |
| Carbapenem resistance prevalence | **40.7%** (200/492; 95% CI 36.4–45.0%) |
| AMR loci CARD-verified (phmmer) | **217 / 217 (100%)** |
| Composite transposon rate (pair level) | **94.6%** (2,825/2,985) |
| IS6-family dominance | **36.8%** of all classified IS elements |
| 2025 cohort IS6 AUC | **0.976** (Cliff's δ = 0.952, p = 2.2×10⁻²³) |
| NCBI submission bias (Fisher OR) | **48.0** (p = 2.0×10⁻⁴⁶) |
| Dual carbapenemase genomes | **5** (KPC + NDM co-carriage) |

### Pipeline Architecture

```
NCBI GenBank (677 assemblies)
      │
      ▼
00_data_audit.py      ← 7-criterion integrity audit (0 duplicates, 0 corrupted)
      │
      ▼
01_download.py        ← Entrez API retrieval (pagination bug fixed)
02_validate.py        ← QC: size / N50 / CDS / format  →  492 PASS
      │
      ▼
03_amr_scan.py        ← Tier1: AMRFinder / Tier2: GFF3 regex → 217 hits
04_is_context.py      ← IS flanking analysis ±10 kb → 2,985 IS–AMR pairs
05_stats.py           ← Wilson CI / Clopper-Pearson / Mann-Whitney / Cliff's δ
06_figures.py         ← Publication-quality figures (PDF + PNG)
07_manuscript.py      ← Auto-populated draft
      │
      ├── 08_is_burden_all.py    ← Per-genome IS6 burden → AUC / ROC
      ├── 09_is_hmmer_verify.py  ← PFAM PF01527 HMM reclassification (pyhmmer)
      ├── 10_amr_hmmer_verify.py ← phmmer vs CARD v3.3 → 217/217 CONFIRMED
      └── 11_subgroup_analysis.py← Temporal stratification + submission bias
```

---

## Repository Structure

```
amr-genomics-analysis/
├── analysis/
│   ├── config.py               # All biological thresholds — edit only this file
│   ├── 00_data_audit.py        # Data integrity audit (7 checks)
│   ├── 01_download.py          # NCBI Entrez download (pagination fix applied)
│   ├── 02_validate.py          # Genome quality control
│   ├── 03_amr_scan.py          # AMR gene detection (Tier1 + Tier2)
│   ├── 04_is_context.py        # IS element flanking analysis
│   ├── 05_stats.py             # Statistical analysis
│   ├── 06_figures.py           # Figure generation
│   ├── 07_manuscript.py        # Auto manuscript generation
│   ├── 08_is_burden_all.py     # IS6 burden analysis + AUC
│   ├── 09_is_hmmer_verify.py   # IS6 HMMER reclassification
│   ├── 10_amr_hmmer_verify.py  # AMR CARD sequence verification
│   └── 11_subgroup_analysis.py # Temporal stratification + bias quantification
├── docs/
│   ├── WORKFLOW.md             # Pipeline workflow reference
│   ├── OPERATION_MEMORY.md     # Decision log and verified numbers
│   └── SYSTEM_WHITEPAPER.md    # System architecture and scientific rationale
├── reports/
│   ├── manuscript_submission_v1.md   # Submission-ready manuscript (Markdown)
│   └── manuscript_submission_v1.pdf  # Submission-ready manuscript (PDF)
├── data/                       # Not tracked in git (see .gitignore)
│   ├── raw/                    # Downloaded genome files (~677 directories)
│   ├── validated/              # QC-passed manifest (492 genomes)
│   ├── processed/              # Analysis outputs (TSV / JSON)
│   └── db/                     # Reference databases (CARD, PFAM)
├── figures/                    # Generated figures (PDF + PNG)
├── logs/                       # Per-step run logs
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Quickstart

### Prerequisites

```bash
conda create -n amr_env python=3.12
conda activate amr_env
pip install -r requirements.txt
```

### Reference Databases (manual step)

```bash
# CARD v3.3 Protein Homolog Model (required for Step 10)
# Download from: https://card.mcmaster.ca/download
# File: broadstreet-v3.3.0.tar.bz2
tar -xjf broadstreet-v3.3.0.tar.bz2
cp protein_fasta_protein_homolog_model.fasta data/db/card/

# PFAM PF01527 (IS6 family HMM — auto-downloaded by Step 09)
# Cached at: data/db/pfam/PF01527.hmm
```

### Configuration

Edit `analysis/config.py`:
```python
NCBI_EMAIL   = 'your.email@example.com'   # Required
NCBI_API_KEY = 'your_api_key'             # Optional (10× rate increase)
MAX_GENOMES_CHINA = None                  # None = retrieve all (~677)
```

### Run the Full Pipeline

```bash
# Step 0: Data audit (run after download)
python analysis/00_data_audit.py

# Steps 1–7: Core pipeline
python analysis/01_download.py    # ~2 hours for 677 genomes
python analysis/02_validate.py
python analysis/03_amr_scan.py
python analysis/04_is_context.py
python analysis/05_stats.py
python analysis/06_figures.py
python analysis/07_manuscript.py

# Steps 8–11: Extended analysis
python analysis/08_is_burden_all.py
python analysis/09_is_hmmer_verify.py
python analysis/10_amr_hmmer_verify.py
python analysis/11_subgroup_analysis.py
```

### One-Liner (after download)

```bash
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
```

---

## Methods Summary

| Step | Tool / Method | Key Output |
|------|---------------|-----------|
| Genome retrieval | NCBI Entrez API (`esearch` + `esummary`) | 677 UIDs |
| Integrity audit | 7-criterion gzip/GFF3/FASTA/size/N50/CDS/duplicate check | 492 PASS |
| AMR detection | GFF3 annotation scanning (Tier 2) with allele-anchored regex | 217 carbapenem hits |
| IS flanking | ±10-kb window, COMPOSITE_TRANSPOSON if upstream ≥1 AND downstream ≥1 | 2,985 pairs |
| AMR verification | phmmer vs CARD v3.3 protein homolog models | 217/217 CONFIRMED |
| IS reclassification | pyhmmer hmmsearch vs PFAM PF01527 (IS6 family HMM) | 1,002 reclassified |
| Temporal bias | Fisher's exact test, cohort stratification | OR = 48.0, p = 2×10⁻⁴⁶ |
| IS6 predictor | Mann-Whitney U → AUC, Cliff's δ | AUC = 0.976 (2025 cohort) |

---

## Known Limitations

- 181 GCF-prefixed assemblies excluded (FTP path naming issue — see `docs/OPERATION_MEMORY.md` L-002)
- IS element detection is annotation-dependent (GFF3 quality varies)
- No province-level geographic metadata (BioSample API not queried in this version)

---

## Citation

If you use this pipeline or data, please cite:

> ZJY. (2026). *IS26-Mediated Composite Transposons Underlie the Near-Universal Flanking Architecture of Carbapenem Resistance Genes in Chinese Clinical Klebsiella pneumoniae: A Population Genomic Analysis.* Preprint / Submitted to Nature Communications.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

**ZJY** — Independent Researcher  
Email: jiayu6954@gmail.com  
GitHub: [https://github.com/jiayu6954](https://github.com/jiayu6954)
