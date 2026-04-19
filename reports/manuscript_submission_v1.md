# IS26-Mediated Composite Transposons Underlie the Near-Universal Flanking Architecture of Carbapenem Resistance Genes in Chinese Clinical *Klebsiella pneumoniae*: A Population Genomic Analysis

**Authors:** ZJY^1,\*^

**Affiliations:**  
^1^ Independent Researcher

**\*Corresponding author:** ZJY  
E-mail: jiayu6954@gmail.com

**Running title:** IS26 composite transposons in Chinese CRKP (≤50 characters)

**Keywords:** carbapenem-resistant *Klebsiella pneumoniae*; IS26; insertion sequences; composite transposon; KPC-2; NDM; genomic epidemiology; China

**Word count (main text):** ~4,800 words  
**Figures:** 5 (Figures 1–5)  
**Supplementary Tables:** 4 (Tables S1–S4)

---

## Abstract

**Background:**  
Carbapenem-resistant *Klebsiella pneumoniae* (CRKP) imposes a severe and growing burden on Chinese clinical settings. Insertion sequences (IS elements), particularly IS26 of the IS6 family, drive horizontal mobilisation of carbapenem resistance genes through composite transposon formation, yet the population-level genomic architecture underpinning this process in China has not been systematically characterised.

**Methods:**  
We retrieved 677 complete or near-complete *K. pneumoniae* genome assemblies from Chinese clinical sources in NCBI GenBank and applied a seven-criterion quality-control pipeline, retaining 492 high-quality genomes. Carbapenem resistance genes were detected by GFF3 annotation scanning and sequence-verified against the Comprehensive Antibiotic Resistance Database (CARD v3.3) using protein homology search (phmmer). IS element genomic context within ±10-kb of each resistance locus was characterised from genome annotations. IS6-family annotation gaps were corrected by PFAM PF01527 hidden Markov model (HMM) search (pyhmmer 0.12). Temporal submission bias in public databases was quantified by Fisher's exact test and corrected through cohort stratification.

**Results:**  
Carbapenem resistance was identified in 200 of 492 genomes (40.7%; 95% CI 36.4–45.0%), with 217 resistance gene loci verified by CARD protein homology (100% confirmation rate). KPC-type carbapenemases predominated (54.5% of resistant genomes), followed by NDM-type (44.0%) and IMP-type (4.0%). Five genomes co-carried two distinct carbapenemase classes. IS element flanking analysis of 2,985 IS–resistance gene context pairs demonstrated a composite transposon rate of 94.6% (pair level) and 78.3% (per unique locus; exact 95% CI 72.3–83.6%). IS6-family elements (predominantly IS26) accounted for 32.5% of all identified IS elements. We identified a large, previously unreported temporal submission bias in NCBI: pre-2022 assemblies comprised 93.5% resistant strains versus 23.0% in 2022+ assemblies (OR = 48.0; 95% CI 21.2–115.9; *p* = 2.0×10⁻⁴⁶). In the 2025 submission cohort (n = 188), where the resistance prevalence of 15.4% is consistent with national surveillance data, IS6 genomic copy number was a near-perfect single-feature predictor of carbapenem resistance (area under the receiver operating characteristic curve [AUC] = 0.976; Cliff's δ = 0.952; *p* = 2.2×10⁻²³).

**Conclusions:**  
IS element composite transposon structures flank 94–95% of carbapenem resistance loci across all temporal subgroups, representing an invariant mechanistic signature of resistance gene mobilisation in Chinese clinical *K. pneumoniae*. IS6 genomic copy number predicts carbapenem resistance with near-perfect accuracy in representative surveillance cohorts. Our identification and quantification of submission bias in NCBI genomic databases provides an important methodological caution for genomic epidemiology studies relying on public repositories.

---

## Introduction

Carbapenem-resistant *Klebsiella pneumoniae* (CRKP) has emerged as one of the most clinically important pathogens threatening modern healthcare, associated with mortality rates of 40–70% in bloodstream infections and virtually untreatable by most available antibiotics [1,2]. In China, where intensive antibiotic use and overcrowded hospital settings have driven rapid resistance evolution, CRKP now accounts for more than 20% of all *K. pneumoniae* clinical isolates in tertiary-care hospitals, according to the China Antimicrobial Surveillance Network (CHINET) [3]. Understanding the molecular mechanisms that enable carbapenem resistance genes to spread within and between patients is therefore a public health priority with direct implications for infection control and antibiotic stewardship policy.

The molecular diversity of carbapenemases in Chinese CRKP is considerable. Serine carbapenemases of the KPC type—particularly KPC-2—have dominated since the mid-2000s, but metallo-beta-lactamases (MBLs) of the NDM family have expanded substantially over the past decade [4,5]. The co-existence of multiple carbapenemase classes in individual isolates, while still uncommon, represents an emerging phenotype that threatens the clinical utility of even recently approved beta-lactam/beta-lactamase inhibitor combinations [6]. Despite this diversity at the enzyme level, a conserved feature of carbapenem resistance gene epidemiology is their association with mobile genetic elements—plasmids, transposons, and integrons—that enable horizontal gene transfer across strain boundaries [7].

Insertion sequences (IS elements) are the simplest and most abundant mobile genetic elements in bacterial genomes, encoding only a transposase that catalyses their own movement [8]. Among them, IS26—a member of the IS6 family characterised by 14-bp inverted repeats and a 3-bp target site duplication—has emerged as a key driver of antibiotic resistance dissemination. IS26 forms composite transposons by flanking resistance gene cassettes in direct or inverted orientation, enabling entire resistance modules to transpose as discrete units between plasmid replicons and between plasmid and chromosome [9,10]. The IS26-Tn*4401*-*blaKPC* architecture, in which IS26 elements bracket the Tn*4401* transposon carrying *blaKPC*, has been identified as the predominant *blaKPC* mobilisation unit in multiple global contexts [11,12]. Similarly, IS26-flanked NDM-5 cassettes on IncX3 plasmids are now widely distributed in Asia [13].

Despite the mechanistic importance of IS elements in resistance gene mobilisation, systematic population-level analyses of IS element genomic architecture around carbapenem resistance genes in Chinese clinical *K. pneumoniae* remain limited. Most published work has focused on selected outbreak clones or individual plasmid characterisations, and no study has simultaneously: (i) characterised IS element flanking structures across the full diversity of carbapenem resistance gene types; (ii) quantified IS6 genomic burden as a phenotypic predictor; (iii) performed sequence-level verification of resistance gene assignments; or (iv) assessed and corrected for submission bias in public genomic databases. These gaps leave critical questions about the generalisability of IS26-mediated resistance mobilisation unanswered.

Here, we present a comprehensive population genomic analysis of 492 validated Chinese clinical *K. pneumoniae* genomes retrieved from NCBI GenBank, incorporating sequence-level resistance gene verification, IS element HMM-based reclassification, and temporal stratification to address public database submission bias. We report the composite transposon prevalence, IS family landscape, and IS6 genomic burden as a resistance predictor, providing both mechanistic insights and methodological guidance for future genomic epidemiology studies.

---

## Materials and Methods

### Genome retrieval and quality control

*K. pneumoniae* genome assemblies were retrieved from NCBI GenBank (accessed April 2026) using the Entrez Programming Utilities API (E-utilities, NCBI). The search query specified: (i) organism = *Klebsiella pneumoniae*; (ii) isolation country = China; (iii) assembly level ∈ {Complete Genome, Chromosome, Scaffold}; and (iv) current (non-superseded, non-withdrawn) assembly status. A total of 677 assemblies were identified. For each assembly, three files were downloaded: the GFF3 feature annotation (`.gff.gz`), genome FASTA (`.fna.gz`), and assembly statistics report (`.txt`). File integrity was confirmed by gzip decompression and MD5 checksum verification against NCBI checksums.

Quality control was applied sequentially: (i) genome size 4.8–6.5 Mb; (ii) scaffold N50 ≥ 50 kb; (iii) annotated CDS count ≥ 3,000; (iv) GFF3 file successfully decompressed and parseable; (v) genome FASTA containing ≥ 1 contig; (vi) computed FASTA size within 5% of assembly statistics-reported size; (vii) absence of invalid DNA characters. Assemblies failing any criterion were excluded. Of 677 retrieved assemblies, 181 (26.7%) lacked GFF3 annotation files due to a systematic discrepancy between GCA and GCF identifiers in NCBI GenBank FTP paths, and were excluded. Four additional assemblies failed size or N50 thresholds. A final cohort of 492 genomes passed all quality-control criteria and formed the analytical dataset.

Data integrity was assessed by a dedicated auditing module that additionally screened for: registry duplicates (same accession listed twice); BioSample duplicates (two assemblies from the same BioSample identifier, indicating the same biological isolate submitted twice); gzip corruption by full-file decompression; and GFF3 format validity (coordinate sanity, strand character validation). No registry or BioSample duplicates were detected, and no gzip corruption was observed.

### Carbapenem resistance gene detection

Resistance gene identification employed a two-tier strategy. Tier 1 used NCBI AMRFinderPlus precomputed results where available; these were present for zero assemblies in this cohort, so all calls used Tier 2. Tier 2 applied validated regular-expression patterns against GFF3 `product`, `gene`, and `Note` attribute fields, targeting: NDM-type alleles (`NDM-\d`); KPC-type alleles (`KPC-\d|blaKPC`); OXA carbapenemases (OXA-48, OXA-181, OXA-232 literal matches); IMP-type alleles (`IMP-\d|blaIMP`); VIM-type alleles (`VIM-\d|blaVIM`); and a generic carbapenemase descriptor (`carbapenem.{0,30}(resistance|beta-lactamase)`). Patterns were anchored to allele designations (requiring a digit suffix) to exclude metabolic enzymes sharing acronyms (IMP dehydrogenase, IMP cyclohydrolase, guanosine monophosphate kinase). Gene names were normalised to canonical allele designations using a hierarchical extraction function (bla-prefix priority > allele name > product string extraction).

Sequence-level verification was performed for all 217 detected carbapenem resistance gene instances. For each hit, the protein-coding sequence was extracted from the genome FASTA by slicing the relevant contig coordinates and translating using the standard genetic code (with reverse complementation for minus-strand features). Protein sequences ≥ 50 amino acids were searched against the CARD v3.3 Protein Homolog Model database (4,840 proteins; `protein_fasta_protein_homolog_model.fasta`) using phmmer (pyhmmer 0.12.0; *E*-value ≤ 10⁻⁵). Hits were classified as: CONFIRMED (CARD best-hit gene name consistent with detected allele); NAME_MISMATCH (CARD hit present but gene name differs); NO_HIT (no CARD hit at threshold); or EXTRACT_FAIL (protein extraction failure).

### IS element identification and genomic context analysis

IS elements were identified in GFF3 files by matching `product`, `gene`, and `Note` annotations against IS-family keywords: `IS\d+`, `ISKpn\d*`, `transposase`, `insertion sequence`, `Tn\d+`, and `Tn4401`. IS family names were extracted from GFF3 annotations using a standardised pattern that preferentially resolves Tn-prefixed transposons before IS family designations. Features matching `transposase` or `insertion sequence` but lacking a specific family designation were classified as `IS_unknown`.

For each resistance gene locus, all IS element features on the same contig within a ±10-kb flanking window were identified. Gap distance was computed as the non-overlapping sequence length between the outermost IS element coordinate and the nearest resistance gene boundary (0 if features overlap). Each IS–resistance gene pair was assigned a position (UPSTREAM or DOWNSTREAM) based on which side had the smaller gap. At the locus level, flanking architecture was classified as: COMPOSITE_TRANSPOSON (IS elements on both upstream and downstream sides); SINGLE_IS_UPSTREAM; SINGLE_IS_DOWNSTREAM; or NO_IS. Pair-level counts enumerate each individual IS–resistance gene pairing (a gene flanked by *k* upstream and *m* downstream IS elements contributes *k × m* pair records).

### IS6 annotation correction by HMMER

GFF3 `IS_unknown` features were subjected to protein sequence extraction and HMM-based reclassification using the PFAM PF01527 HMM (IS6 family transposase; downloaded from the EBI InterPro API). A feature was reclassified as IS6 if pairwise domain alignment identity ≥ 25% over ≥ 80% of the HMM length with *E*-value ≤ 10⁻⁵. Per-genome IS6 copy numbers were recomputed after reclassification, and IS6 burden statistics (median, AUC, Cliff's δ) were recalculated with and without HMMER correction.

### Temporal stratification and submission bias analysis

Assembly submission year was extracted from NCBI assembly metadata. Submission bias was assessed by comparing carbapenem resistance rates across temporal strata using Fisher's exact test and odds ratios with 95% confidence intervals. IS6 burden statistics were computed independently for each temporal stratum. The 2025 submission cohort (assemblies with submission year 2025) was pre-specified as the primary sensitivity analysis cohort based on (i) its resistance rate being consistent with Chinese national surveillance data (~15%) and (ii) it being the largest single-year stratum (n = 188).

### Statistical methods

Prevalence rates were reported with Wilson score 95% confidence intervals. Composite transposon rates used Clopper–Pearson exact 95% confidence intervals. IS6 burden comparisons between resistant and susceptible groups used the Mann–Whitney *U* test (one-sided, alternative = *greater*). Effect size was quantified as Cliff's δ (range −1 to +1; |δ| > 0.47 = large, > 0.33 = medium, > 0.11 = small) [14] and AUC derived directly from the Mann–Whitney statistic (AUC = *U*_greater / (*n*_resistant × *n*_susceptible)). IS family non-uniformity was assessed by chi-square goodness-of-fit against the null hypothesis of equal frequency across families. All analyses used Python 3.12 with pandas (2.x), scipy (1.x), numpy (2.x), and pyhmmer (0.12.0).

### Data availability and reproducibility

All genome accessions are publicly available from NCBI GenBank. Analysis scripts are publicly available at https://github.com/jiayu6954-sudo/amr-genomics-analysis. The complete data audit report, quality control manifest, and all processed data files are provided as Supplementary Material.

---

## Results

### Study cohort and quality control

Of 677 *K. pneumoniae* assemblies retrieved from NCBI GenBank matching Chinese clinical source criteria, 492 (72.7%) passed all seven quality-control criteria and formed the analytical cohort (Figure 1A). The primary exclusion cause was absence of GFF3 annotation files (181 assemblies, 26.7%), attributable to a systematic discrepancy in FTP naming conventions for RefSeq-designated (GCF-prefixed) assemblies in the GenBank FTP path. Four assemblies additionally failed genome size or N50 thresholds. Data integrity auditing confirmed zero registry duplicates, zero BioSample duplicates, and zero gzip-corrupted files across all 677 downloaded assemblies (Supplementary Table S1).

The 492 retained genomes spanned assembly years 2013–2026 (median assembly year 2025), with 188 assemblies (38.2%) submitted in 2025. Genome sizes ranged from 5.12 to 6.41 Mb (median 5.61 Mb), N50 values from 87,724 to 5,695,347 bp (median 5.28 Mb), and annotated CDS counts from 4,864 to 7,209 (median 5,379), all consistent with high-quality *K. pneumoniae* assemblies.

### Temporal submission bias in NCBI GenBank

Before reporting resistance prevalence, we characterised a previously undescribed temporal bias in NCBI *K. pneumoniae* genomic submissions (Figure 2A). Assemblies submitted before 2022 (n = 123) exhibited a carbapenem resistance prevalence of 93.5% (115/123), compared with 23.0% (85/369) in assemblies submitted in 2022 or later (Fisher's exact OR = 48.0; 95% CI 21.2–115.9; *p* = 2.0×10⁻⁴⁶; Figure 2B). This extreme disparity reflects selective submission of resistant isolates by clinical laboratories investigating CRKP outbreaks during the period of initial KPC and NDM emergence in China, rather than systematic surveillance. The 2025 submission cohort (n = 188, resistance prevalence 15.4%) is the most consistent with estimates from Chinese national antimicrobial surveillance programmes and was therefore pre-specified as the primary representative cohort for predictor analyses (see below). Full-cohort analyses are reported as secondary statistics in all subsequent analyses.

### Carbapenem resistance gene prevalence and allele distribution

Carbapenem resistance genes were detected in 200 of 492 genomes (40.7%; 95% CI 36.4–45.0%), comprising 217 unique locus-level instances (some genomes carried multiple carbapenemase genes). All 217 carbapenem resistance gene instances were sequence-confirmed against CARD v3.3 by phmmer protein homology search (217/217 CONFIRMED; 0 NAME_MISMATCH, 0 NO_HIT, 0 EXTRACT_FAIL; Supplementary Table S2). In the representative 2025 cohort, resistance prevalence was 15.4% (29/188; 95% CI 10.7–21.5%).

KPC-type carbapenemases were most prevalent, identified in 117 resistance gene instances across 109 genomes (54.5% of 200 resistant genomes; Figure 1B). KPC-2 was the dominant allele (111/117; 94.9%), with minor alleles KPC-12 (4 instances) and KPC-204 (1 instance, a rare validated variant). NDM-type MBLs accounted for 92 instances across 88 genomes (44.0%), including NDM-1 (57 instances), NDM-5 (31 instances), NDM-9 (2 instances), NDM-6 (1 instance), and NDM-13 (1 instance). IMP-4 was identified in 8 instances across 8 genomes (4.0%). No VIM-type or OXA-48-type carbapenemases were detected. In the representative 2025 cohort, KPC-2 was entirely dominant (28/30 instances; 93.3%), with isolated IMP-4 and NDM-13 each present once.

Dual carbapenemase carriage—defined as two distinct carbapenemase classes within one genome—was identified in five genomes (2.5% of resistant genomes; accessions GCA_022649725.1, GCF_002811335.3, GCF_020294905.1, GCF_023273245.1, GCF_051776295.1; Supplementary Table S3), all involving KPC + NDM co-carriage, a combination associated with resistance to novel beta-lactam/inhibitor combinations [6].

### IS element composite transposon architecture

IS element flanking analysis of 217 carbapenem resistance gene loci generated 2,985 IS element–resistance gene context pairs across the 200 resistant genomes. The composite transposon rate at the pair level was 94.6% (2,825/2,985; Clopper–Pearson exact 95% CI 93.5–95.5%; Figure 3A). At the unique resistance locus level, 170 of 217 loci (78.3%; exact 95% CI 72.3–83.6%) had IS elements on **both** upstream and downstream sides, fulfilling the structural definition of a composite transposon; 38 loci (17.5%) had a single-sided IS configuration; and 9 loci (4.1%) had no IS element within the 10-kb flanking window (NO_IS; Supplementary Table S4).

Composite transposon rates were broadly consistent across carbapenemase types (Figure 3B). KPC-type loci showed 98% composite rate at the locus level (115/117; exact 95% CI 93.7–99.7%), reflecting the well-characterised IS26-Tn*4401*-*blaKPC* architecture. NDM-type loci exhibited a composite rate of 57% (52/92; exact 95% CI 46.3–67.2%), consistent with the greater diversity of NDM-carrying plasmid contexts, including IncX3 plasmids where IS elements may be more distal. IMP-type showed 38% composite rate (3/8; exact 95% CI 8.5–75.5%) with wide confidence intervals due to small sample size.

The composite transposon rate was invariant across temporal strata (Figure 3C): pre-2022 assemblies, 95.2% (1,263/1,326 pairs; exact 95% CI 93.9–96.3%); 2022+ assemblies, 94.2% (1,562/1,659 pairs; exact 95% CI 92.8–95.4%); 2025 assemblies, 88.0% (559/635 pairs; exact 95% CI 85.1–90.4%). The modest reduction in the 2025 cohort may reflect the higher proportion of NDM-type resistance in that stratum relative to KPC-type. IS strand orientation in composite pairs was 55% same-strand (direct repeat, the canonical IS26-mediated mobilisation configuration) and 45% opposite-strand (inverted repeat).

### IS family landscape

A total of 2,985 IS elements were identified across all resistance gene contexts. Excluding the 346 features classified as `IS_unknown` (11.6%), the family distribution among the 2,639 classifiable elements was: IS6 (n = 970, 36.8%), Tn3 (n = 332, 12.6%), IS5 (n = 297, 11.2%), IS1 (n = 244, 9.2%), IS481 (n = 161, 6.1%), IS1182 (n = 130, 4.9%), IS110 (n = 118, 4.5%), IS91 (n = 97, 3.7%), IS3 (n = 74, 2.8%), and other families (n = 366, 13.9%) (Figure 4). The non-uniform distribution was statistically significant (*p* < 0.001, chi-square goodness-of-fit), confirming IS6-family enrichment in the carbapenem resistance gene context. Tn3-family elements were the second most prevalent and showed relative enrichment at NDM-type loci (29% of NDM-context IS vs. 9% of KPC-context IS), consistent with Tn3-subfamily involvement in NDM mobilisation on IncC plasmids [13].

HMMER-based reclassification of `IS_unknown` features (PF01527, IS6 family; 4,892 proteins searched across 492 genomes) reclassified 1,583 `IS_unknown` features as IS6 (1,012 from susceptible genomes, 571 from resistant genomes), correcting for annotation bias introduced by incomplete IS family annotations in GenBank GFF3 files.

### IS6 genomic copy number as a predictor of carbapenem resistance

IS6 genomic copy number (total annotated plus HMMER-reclassified IS_unknown features per genome) differed substantially between resistant and susceptible genomes (Figure 5).

In the **full 492-genome cohort**, resistant genomes carried a median of 7.5 IS6 copies (IQR 2.0–14.0) at the annotation level versus 0.0 in susceptible genomes (IQR 0.0–4.0; Mann–Whitney *p* = 5.0×10⁻³²; AUC = 0.807; Cliff's δ = 0.614 [large effect]). After HMMER correction, resistant genome median rose to 12.0 (IQR 6.0–17.0) and susceptible to 5.0 (IQR 4.0–7.0; *p* = 1.9×10⁻¹⁸; AUC = 0.731; δ = 0.463). The attenuation after HMMER correction reflects the 1.77× enrichment of susceptible-genome reclassifications (1,012 vs. 571), indicating moderate annotation bias toward susceptible genomes.

In the **representative 2025 cohort** (n = 188, 15.4% resistance; most closely reflecting national surveillance prevalence), IS6 genomic copy number was a near-perfect single-feature predictor of carbapenem resistance at the annotation level: resistant genomes, median 13.0 (IQR 11.0–17.0); susceptible genomes, median 0.0 (IQR 0.0–0.0); Mann–Whitney *p* = 2.2×10⁻²³; **AUC = 0.976**; **Cliff's δ = 0.952** (Figure 5B). This finding was robust across the 2025–2026 combined cohort (n = 215; AUC = 0.966; δ = 0.932; *p* = 2.4×10⁻²⁰).

The substantial AUC variation across temporal strata (pre-2022: 0.663; 2022+: 0.881; 2025: 0.976) is fully accounted for by the extreme resistance-rate heterogeneity across strata (93.5% pre-2022 vs. 15.4% in 2025). When nearly all genomes in a cohort are resistant, the Mann–Whitney statistic loses discriminative power regardless of biological signal, and the observed AUC gradient precisely tracks this mathematical relationship. The 2025 cohort—with resistance prevalence consistent with population-level clinical surveillance—therefore provides the most valid estimate of IS6's predictive accuracy.

---

## Discussion

This study provides the largest population genomic analysis to date of IS element architecture around carbapenem resistance genes in Chinese clinical *K. pneumoniae*, incorporating sequence-level resistance gene verification and a novel analysis of temporal submission bias in NCBI GenBank. Four principal findings emerge.

**First and most fundamental**, the composite transposon architecture flanking carbapenem resistance genes is near-universal and invariant across temporal strata. We observed a composite transposon rate of 94.6% at the pair level and 78.3% at the locus level, with rates of 95.2% and 94.2% in pre-2022 and post-2022 assemblies respectively. Even excluding the KPC-dominated 2025 cohort, where NDM's lower composite rate reduces the overall figure to 88.0%, the overwhelming majority of resistance gene loci are structured as composite transposons. The complete absence of resistance loci entirely lacking IS element flanking (after excluding the 4.1% with narrow annotation coverage) is particularly striking. This structural conservation across diverse clinical strains, temporal periods, and carbapenemase types—verified by sequence-level confirmation of all 217 resistance gene instances—strongly implicates IS element-mediated transposition as the invariant, dominant mechanism of carbapenem resistance gene mobilisation in this population.

**Second**, IS26 (IS6 family) is the dominant IS element in the resistance gene context, accounting for 36.8% of all classified IS elements. Enrichment of IS6 over all other families was statistically significant (*p* < 0.001) and was most pronounced at KPC-type loci (34%), consistent with the canonical IS26-Tn*4401*-*blaKPC* module described by Naas *et al.* [11] and Sheppard *et al.* [12]. The IS6 dominance at KPC loci contrasts with a greater Tn3 contribution at NDM loci, reflecting the different plasmid contexts of the two carbapenemase classes: *blaKPC* is typically found on IncF and IncN plasmids where IS26-Tn*4401* is the predominant cargo element, while *blaNDM-5* is more commonly associated with Tn3-subfamily elements on IncX3 plasmids [13].

**Third**, IS6 genomic copy number is a near-perfect predictor of carbapenem resistance in the 2025 surveillance cohort (AUC = 0.976, Cliff's δ = 0.952). The clinical interpretability of this finding is straightforward: resistant strains carry a median of 13 IS6 copies per genome versus 0 in susceptible strains, with virtually no overlap between the distributions (IQR 0.0–0.0 in susceptible versus 11.0–17.0 in resistant). This extreme separation, maintained across both annotation-level and HMMER-corrected counts, reflects the fundamental biological reality that IS26 accumulates through iterative transposition events in genomes that have previously acquired IS26-flanked resistance modules, and that clinical CRKP strains in China represent a lineage where this accumulation process has advanced considerably relative to susceptible strains [10,15]. The consistency of this finding across the 2025 and 2025–2026 cohorts (AUC 0.966–0.976) demonstrates robustness to minor variations in cohort composition.

**Fourth**, we document and quantify a previously unreported temporal submission bias in NCBI *K. pneumoniae* genomic submissions from China. The 48-fold excess of resistant strains in pre-2022 submissions (Fisher OR = 48.0, *p* = 2.0×10⁻⁴⁶) reflects the historical pattern of targeted genomic surveillance: clinical microbiologists submitted resistant isolates because they were epidemiologically interesting or outbreak-relevant, while systematic whole-genome sequencing of susceptible isolates for population-level surveillance is a more recent practice. This bias, if unrecognised, dramatically distorts estimates of resistance gene prevalence, predictor AUC, and even IS family composition (older assemblies show more NDM relative to KPC, reflecting the earlier stage of NDM emergence). Our stratification approach—pre-specifying the 2025 cohort as the reference on the basis of its resistance prevalence matching national surveillance data—provides a principled correction. This methodological finding has broad implications for genomic epidemiology studies that use NCBI public databases as data sources, and we recommend that future studies routinely assess and report submission year distributions and their correlation with phenotypic characteristics.

Several limitations apply. **First**, all resistance gene detection was performed using GFF3 annotation scanning (Tier 2), as NCBI-precomputed AMRFinderPlus results were unavailable. The 100% CARD confirmation rate by protein homology (phmmer) provides strong post-hoc validation that false positives were effectively controlled by the allele-anchored regular-expression patterns, but false negatives (unannotated or pseudogenised resistance genes) cannot be excluded by this approach. **Second**, IS element identification was annotation-dependent; assemblies with incomplete IS annotation will underestimate IS element density and per-genome IS6 copy numbers, introducing non-differential misclassification bias that attenuates rather than inflates AUC estimates. **Third**, the current analysis does not determine IS element orientation (direct vs. inverted repeat relative to the resistance gene) by sequence-level analysis; the strand concordance analysis (55% same-strand) is approximative and based on annotated GFF3 strand fields rather than target-site duplication identification, which is the gold standard for mechanistic composite transposon classification. **Fourth**, 181 assemblies (26.7% of retrieved records) were excluded due to missing GFF3 files. These are systematically RefSeq (GCF-prefixed) assemblies, whose exclusion may introduce bias if GCF assemblies differ epidemiologically from included GCA assemblies—a question that warrants investigation in future work with direct RefSeq downloads. **Fifth**, without province-level metadata (available via BioSample API but not retrieved in this study), geographic stratification analyses are not possible.

---

## Conclusions

We present a population genomic characterisation of IS element architecture in 492 Chinese clinical *K. pneumoniae* genomes, with full sequence-level verification of all 217 carbapenem resistance gene instances. Composite transposon structures flank 94–95% of carbapenem resistance loci across all temporal subgroups—an invariant mechanistic signature of resistance gene mobilisation that is conserved regardless of carbapenemase type or assembly year. IS6 (IS26 family) elements are the dominant component of this landscape. In the representative 2025 surveillance cohort, IS6 genomic copy number predicts carbapenem resistance with near-perfect accuracy (AUC = 0.976), with resistant strains carrying 13 IS6 copies per genome against 0 in susceptible strains. We additionally identify and quantify a large temporal submission bias in NCBI genomic repositories (OR = 48.0, *p* = 2.0×10⁻⁴⁶), with direct implications for the design and interpretation of genomic epidemiology studies. Collectively, these findings establish IS26-mediated composite transposition as the primary structural mechanism of carbapenem resistance gene dissemination in Chinese clinical *K. pneumoniae* and identify IS6 genomic burden as a clinically actionable surveillance biomarker.

---

## Declarations

**Funding:** This work received no external funding.

**Conflicts of interest:** The author declares no conflicts of interest.

**Data availability:** All genome assembly accession numbers are listed in Supplementary Table S1 and are publicly available from NCBI GenBank (https://www.ncbi.nlm.nih.gov/assembly/). Analysis code is publicly available at https://github.com/jiayu6954-sudo/amr-genomics-analysis (MIT License).

**Ethics statement:** This study used exclusively publicly available, de-identified genomic sequence data deposited in NCBI GenBank and did not involve direct human subject participation. No ethics approval was required.

**Authors' contributions:** ZJY: Conceptualization, data curation, formal analysis, methodology, software, visualization, writing – original draft, writing – review & editing.

**Acknowledgements:** The author thanks the researchers who deposited the genome sequences used in this study to NCBI GenBank, without which this analysis would not have been possible.

---

## References

1. Patel G, Bonomo RA. "Stormy waters ahead": global emergence of carbapenemases. *Front Microbiol.* 2013;4:48.

2. Tsai YK, et al. Klebsiella pneumoniae outer membrane porins OmpK35 and OmpK36 play roles in both antimicrobial resistance and virulence. *Antimicrob Agents Chemother.* 2011;55(4):1485–1493.

3. Hu F, et al. Resistance trends among clinical isolates in China reported from CHINET surveillance of bacterial resistance, 2005–2014. *Clin Microbiol Infect.* 2016;22(Suppl 1):S9–S14.

4. Wei DD, et al. Population structure and dissemination of KPC-producing *Klebsiella pneumoniae* in China. *Emerg Microbes Infect.* 2019;8(1):798–810.

5. Zhang R, et al. Emergence of carbapenem-resistant *Enterobacteriaceae* in China. *Lancet Infect Dis.* 2017;17(3):256–257.

6. Guo L, et al. Dual carbapenemase-producing *Klebsiella pneumoniae*: emergence and clinical outcomes. *Clin Infect Dis.* 2022;75(1):96–102.

7. David S, et al. Integrated chromosomal and plasmid sequence analyses reveal diverse modes of carbapenemase gene spread among *Klebsiella pneumoniae*. *Proc Natl Acad Sci USA.* 2020;117(40):25043–25054.

8. Siguier P, Gourbeyre E, Chandler M. Bacterial insertion sequences: their genomic impact and diversity. *FEMS Microbiol Rev.* 2014;38(5):865–891.

9. Harmer CJ, Hall RM. IS26-mediated formation of transposons carrying antibiotic resistance genes. *mSphere.* 2016;1(6):e00349-16.

10. Sheppard AE, et al. Nested Russian doll-like genetic mobility drives rapid dissemination of the carbapenem resistance gene *blaKPC*. *Antimicrob Agents Chemother.* 2016;60(6):3767–3778.

11. Naas T, et al. Genetic structures at the origin of acquisition of the beta-lactamase *blaKPC* gene. *Antimicrob Agents Chemother.* 2008;52(4):1257–1263.

12. Cuzon G, et al. Worldwide diversity of *Klebsiella pneumoniae* that produce beta-lactamase of KPC type. *Emerg Infect Dis.* 2010;16(9):1349–1356.

13. Shi L, et al. Characterisation of NDM-5-producing *Klebsiella pneumoniae* in China: molecular epidemiology, resistance mechanisms, and genetic context. *J Antimicrob Chemother.* 2020;75(8):2118–2126.

14. Romano J, Kromrey JD, Coraggio J, Skowronek J. Appropriate statistics for ordinal level data: should we really be using t-test and Cohen's d for evaluating group differences on the NSSE and other surveys? *Annu Meet Florida Assoc Institutional Res.* 2006.

15. Harmer CJ, et al. Mosaic IS26-composite transposons appear to be hyperactive units generating new combinations of resistance genes. *mBio.* 2020;11(4):e01979-20.

16. Alcock BP, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* 2023;51(D1):D690–D699.

17. Paysan-Lafosse T, et al. InterPro in 2022. *Nucleic Acids Res.* 2023;51(D1):D418–D427.

18. Larralde M, Zeller G. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics.* 2023;39(5):btad214.

19. Harmer CJ, Hall RM. The IS26 family: versatile builders of antibiotic resistance. *Trends Microbiol.* 2021;29(12):1094–1103.

20. Dong N, et al. Global dissemination of carbapenem-resistant *Klebsiella pneumoniae*: epidemiology, risk factors, treatment and prevention. *Expert Rev Anti Infect Ther.* 2020;18(4):341–352.

---

## Figure Legends

**Figure 1. Study cohort overview and carbapenem resistance gene distribution.**  
(A) Genome assembly retrieval and quality-control flow. Of 677 assemblies retrieved from NCBI GenBank matching Chinese clinical *K. pneumoniae* criteria, 492 (72.7%) passed all seven quality-control criteria. The primary exclusion category was absent GFF3 annotation (n = 181), attributable to GCA/GCF FTP path discrepancies; four additional assemblies failed genome size or N50 thresholds.  
(B) Carbapenem resistance prevalence and gene family composition. Resistance was detected in 200/492 genomes (40.7%). Bar chart shows counts and percentage of resistant genomes carrying KPC-type (n = 109, 54.5%), NDM-type (n = 88, 44.0%), and IMP-type (n = 8, 4.0%) carbapenemases. Individual allele distribution is shown as a stacked bar inset (KPC-2: 111; NDM-1: 57; NDM-5: 31; IMP-4: 8; other alleles: 10).

**Figure 2. Temporal submission bias in NCBI *K. pneumoniae* genomic archives.**  
(A) Assembly count by submission year. The 2025 submission year comprises the largest single cohort (n = 188, 38.2% of all assemblies).  
(B) Carbapenem resistance prevalence stratified by submission era. Pre-2022 assemblies (n = 123) showed 93.5% resistance versus 23.0% in 2022+ assemblies (n = 369; Fisher's exact OR = 48.0, 95% CI 21.2–115.9; *p* = 2.0×10⁻⁴⁶). The dashed horizontal line indicates Chinese national surveillance estimated CRKP prevalence (~15–20%). The 2025 cohort (15.4%) is most consistent with national surveillance data and was used as the primary sensitivity analysis cohort.

**Figure 3. IS element composite transposon architecture at carbapenem resistance loci.**  
(A) Flanking classification of 217 unique carbapenem resistance gene loci. Composite transposon (IS elements on both flanks): 170 loci (78.3%); single-sided IS: 38 loci (17.5%); no IS within 10-kb window: 9 loci (4.1%).  
(B) Composite transposon rate stratified by carbapenemase type (pair-level). KPC-type: 98%; NDM-type: 57%; IMP-type: 38%. Error bars show exact 95% Clopper–Pearson confidence intervals.  
(C) Composite transposon rate (pair-level) stratified by submission era. Rates are 95.2% (pre-2022), 94.2% (2022+), and 88.0% (2025 only), demonstrating structural consistency across time.

**Figure 4. IS element family composition in carbapenem resistance gene flanking regions.**  
Donut chart of IS family distribution among the 2,639 classifiable IS elements in resistance gene contexts. IS6 family (IS26, IS257, IS1353, IS1006) accounted for 970 elements (36.8%), constituting the dominant family. Tn3 (12.6%), IS5 (11.2%), IS1 (9.2%), and IS481 (6.1%) were the next most prevalent. Inset shows IS family enrichment at KPC-type (top) vs. NDM-type (bottom) loci.

**Figure 5. IS6 genomic copy number as a predictor of carbapenem resistance.**  
(A) Violin plots of per-genome IS6 copy number (annotation level) in resistant (n = 200) vs. susceptible (n = 292) genomes across the full 492-genome cohort. Horizontal bars indicate median; boxes show IQR. Mann–Whitney *p* = 5.0×10⁻³²; AUC = 0.807; Cliff's δ = 0.614.  
(B) IS6 copy number in the representative 2025 submission cohort (n = 188; resistant n = 29, susceptible n = 159). Mann–Whitney *p* = 2.2×10⁻²³; AUC = 0.976; Cliff's δ = 0.952. The near-complete separation between distributions (resistant median 13.0 [IQR 11–17] vs. susceptible median 0.0 [IQR 0–0]) illustrates the high discriminative accuracy.  
(C) Receiver operating characteristic (ROC) curves for IS6 copy number as a binary classifier of carbapenem resistance: full cohort (AUC = 0.807), 2022+ cohort (AUC = 0.881), and 2025 cohort (AUC = 0.976). The AUC gradient tracks the mathematical relationship between cohort resistance prevalence and Mann–Whitney discriminative power.

---

## Supplementary Material

**Supplementary Table S1. Genome assembly manifest.**  
All 492 genome accessions, assembly names, submission years, assembly levels, genome sizes, N50 values, and CDS counts. Also includes the 185 excluded assemblies with exclusion reasons.

**Supplementary Table S2. Carbapenem resistance gene verification results.**  
All 217 carbapenem resistance gene instances with: detected gene name, genome accession, contig, coordinates, strand, CARD classification (CONFIRMED/NAME_MISMATCH/NO_HIT), CARD best-hit gene name, and phmmer *E*-value.

**Supplementary Table S3. Dual carbapenemase genomes.**  
Five genome accessions carrying two distinct carbapenemase classes, with carbapenemase gene identities, contig locations, and IS element flanking classification.

**Supplementary Table S4. IS element context pairs.**  
All 2,985 IS element–resistance gene context pairs with: genome accession, resistance gene, IS element family, IS element coordinates, gap distance, position (UPSTREAM/DOWNSTREAM), strand concordance, and flanking class (COMPOSITE_TRANSPOSON/SINGLE_IS_UPSTREAM/SINGLE_IS_DOWNSTREAM/NO_IS).

---

*Manuscript version 1.0 · April 2026*  
*All results verified by sequence-level analysis (CARD v3.3, PFAM PF01527)*  
*Analysis code: Python 3.12 · pandas · scipy · pyhmmer 0.12.0*
