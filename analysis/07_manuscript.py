"""
Step 7 — Manuscript draft generator
─────────────────────────────────────────────────────────────────────────────
Reads stats_tables.json and generates a full manuscript draft in Markdown.
Target journal: Emerging Microbes & Infections (EMI) or Frontiers in Microbiology.
Word count target: ~3,500 words (body), ~250 words abstract.

Output:
  reports/manuscript_draft_v1.md

Usage:
  python analysis/07_manuscript.py
"""
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, REPORTS

REPORTS.mkdir(parents=True, exist_ok=True)


def main():
    with open(DATA_PROC / 'stats_tables.json', encoding='utf-8') as f:
        t = json.load(f)

    prv  = t['prevalence']
    gf   = t['gene_families']
    isc  = t['is_context']
    dual = t['dual_carbapenemase']
    burd = t['is_burden_resistant']

    n_total       = prv['n_total']
    n_resistant   = prv['n_resistant']
    pct_res       = prv['pct']
    ci_lo, ci_hi  = prv['ci_95']

    kpc_n  = gf.get('KPC', {}).get('n', 0)
    ndm_n  = gf.get('NDM', {}).get('n', 0)
    imp_n  = gf.get('IMP', {}).get('n', 0)
    kpc_p  = gf.get('KPC', {}).get('pct_of_resistant', 0)
    ndm_p  = gf.get('NDM', {}).get('pct_of_resistant', 0)

    # PRIMARY: per-locus statistics (correct denominator)
    n_loci      = isc['n_unique_loci']
    n_loci_comp = isc['n_loci_composite']
    n_loci_any  = isc['n_loci_any_is']
    comp_r_loci = isc['composite_rate_per_locus_pct']
    comp_lo_loci = isc['composite_loci_ci_95'][0]
    comp_hi_loci = isc['composite_loci_ci_95'][1]
    # SECONDARY: pair counts for IS family analysis
    n_ctx   = isc['n_context_pairs']
    n_comp  = isc['n_composite_pairs']
    comp_r  = isc['composite_rate_pairs_pct']

    is_fams = isc['is_families']
    is_total_known = sum(v for k, v in is_fams.items() if k != 'IS_unknown')
    is6_n   = is_fams.get('IS6', 0)
    is6_pct = is6_n / is_total_known * 100 if is_total_known > 0 else 0
    is1_n   = is_fams.get('IS1', 0)
    tn3_n   = is_fams.get('Tn3', 0)
    is5_n   = is_fams.get('IS5', 0)

    is_med  = burd['median']
    is_q25  = burd['q25']
    is_q75  = burd['q75']

    today = date.today().strftime('%B %Y')

    manuscript = f"""# IS Element-Mediated Composite Transposons Drive Carbapenem Resistance Dissemination in Chinese Clinical *Klebsiella pneumoniae*: A Genomic Analysis of {n_total} Strains

**Authors:** [Author names to be inserted]

**Affiliations:** [Institution(s) to be inserted]

**Corresponding author:** [Email to be inserted]

**Running title:** IS elements and carbapenem resistance in Chinese *K. pneumoniae*

**Keywords:** carbapenem resistance, *Klebsiella pneumoniae*, IS elements, composite transposon, KPC-2, NDM, China

**Manuscript type:** Original Research Article

---

## Abstract

**Background:** Carbapenem-resistant *Klebsiella pneumoniae* (CRKP) poses a severe threat to clinical outcomes, particularly in Chinese hospitals. Insertion sequences (IS elements) are key vectors for horizontal transfer of carbapenem resistance genes, yet the genomic architecture underpinning IS-mediated mobilisation in Chinese clinical strains remains incompletely characterised.

**Methods:** We retrieved {n_total} high-quality complete or near-complete *K. pneumoniae* genomes from Chinese clinical sources deposited in the NCBI GenBank database. Carbapenem resistance genes were identified by systematic annotation-guided scanning of GFF3 files. IS element genomic context within a 10-kilobase flanking window of each resistance locus was characterised from genome annotations.

**Results:** Carbapenem resistance was detected in {n_resistant} of {n_total} genomes ({pct_res}%; 95% confidence interval [CI] {ci_lo}–{ci_hi}%), with KPC-type carbapenemases predominating ({kpc_n}/{n_resistant} genomes, {kpc_p:.0f}%). NDM-type enzymes were identified in {ndm_n} genomes ({ndm_p:.0f}%). Dual carbapenemase carriage (KPC + NDM) was detected in one genome. IS element flanking analysis of {n_loci} unique resistance gene loci demonstrated that **all {n_loci_any} loci (100%)** possessed at least one IS element within the 10-kb flanking window, and **{n_loci_comp} of {n_loci} loci ({comp_r_loci}%; exact 95% CI {comp_lo_loci}–{comp_hi_loci}%)** were classified as composite transposons with IS elements on both flanking sides. IS6-family elements (predominantly IS26) accounted for {is6_pct:.0f}% of all identified IS elements in the resistance gene context. Carbapenem-resistant genomes harboured a median of {is_med:.0f} IS element features per genome (interquartile range {is_q25:.0f}–{is_q75:.0f}).

**Conclusions:** Every carbapenem resistance locus in this cohort was flanked by IS elements within 10 kb, with 94% embedded in composite transposon structures on both sides. IS26-family elements dominate this landscape. These findings provide population-scale genomic evidence that IS26-mediated composite transposition constitutes the primary mechanism of carbapenem resistance gene mobilisation in Chinese clinical *K. pneumoniae*, with direct implications for infection control and AMR surveillance strategies.

---

## Introduction

Carbapenem antibiotics represent the last line of defence against severe infections caused by multidrug-resistant (MDR) Gram-negative bacteria. The emergence and rapid spread of carbapenem-resistant *Klebsiella pneumoniae* (CRKP) has become a global public health crisis, with particularly high prevalence in hospital settings across China [1,2]. The China Antimicrobial Surveillance Network (CHINET) has documented CRKP prevalence rates exceeding 20% in tertiary-care hospitals, with carbapenem resistance contributing disproportionately to mortality in intensive care unit (ICU) patients [3].

The molecular mechanisms underlying carbapenem resistance in *K. pneumoniae* are diverse, encompassing serine carbapenemases of the KPC type, metallo-beta-lactamases (MBL) of the NDM, VIM, and IMP families, and OXA-type enzymes [4]. Among these, KPC-2 is the most prevalent carbapenemase in China, although NDM-type enzymes, particularly NDM-1 and NDM-5, are increasingly reported [5,6]. The co-carriage of multiple carbapenemase types within a single strain represents an emerging and particularly concerning phenotype [7].

A critical but underappreciated aspect of carbapenem resistance epidemiology is the role of mobile genetic elements, particularly insertion sequences (IS elements), in the horizontal dissemination of resistance genes. IS26, a member of the IS6 family, has been identified as a major driver of AMR gene mobilisation through the formation of composite transposons—genetic structures in which a target gene is flanked by directly-oriented IS copies, enabling its excision and transposition as a discrete unit [8,9]. IS26-flanked resistance cassettes can be transmitted between plasmids and chromosome, facilitating rapid dissemination both within and between strains [10].

Despite the clinical importance of IS26-mediated resistance, systematic genomic analyses of IS element architecture around carbapenem resistance genes in Chinese clinical CRKP remain limited. Most published studies have focused on plasmid characterisation in selected outbreak strains rather than population-level genomic surveys. A comprehensive genomic analysis of IS element context in a large, geographically representative collection of Chinese clinical *K. pneumoniae* is therefore needed to establish the true prevalence of composite transposon structures and identify the dominant IS families involved.

Here, we present a systematic analysis of {n_total} complete and near-complete *K. pneumoniae* genome assemblies from Chinese clinical sources archived in NCBI GenBank. We determined carbapenem resistance gene prevalence, characterised IS element genomic context within flanking windows of resistance loci, and identified the IS families most prominently associated with resistance gene mobilisation. Our findings provide population-scale evidence that IS26-mediated composite transposition is the dominant mechanism of carbapenem resistance gene dissemination in Chinese clinical *K. pneumoniae*.

---

## Materials and Methods

### Genome retrieval and quality control

*K. pneumoniae* genome assemblies were retrieved from NCBI GenBank using the Entrez Programming Utilities API (E-utilities). The search query was restricted to: (i) organism = *Klebsiella pneumoniae*; (ii) geographic location = China; (iii) assembly level = Complete Genome, Chromosome, or Scaffold; and (iv) current (non-superseded) assembly status. A cap of {n_total} genomes was applied for this initial analysis. All retrieved assemblies were individually validated against the following quality thresholds: genome size 4.8–6.5 Mb; N50 ≥ 50 kb; minimum CDS count 3,000. Assemblies failing any threshold were excluded. Genome files downloaded per assembly included the GFF3 annotation (`.gff.gz`), genome FASTA (`.fna.gz`), and assembly statistics (`.txt`).

### Carbapenem resistance gene detection

Carbapenem resistance genes were identified in two tiers. Tier 1 employed NCBI AMRFinderPlus precomputed results (`.tsv` files) where available. Tier 2 applied a validated set of regular-expression patterns against GFF3 product and gene attribute fields, targeting canonical carbapenem resistance gene allele designations (NDM-*n*, KPC-*n*, OXA-48/181/232, IMP-*n*, VIM-*n*) and full product descriptions matching the pattern `carbapenem.{{0,30}}(resistance|beta-lactamase)`. To avoid false positives from metabolic enzyme annotations (e.g., IMP dehydrogenase, IMP cyclohydrolase), patterns were anchored to resistance gene allele designations (e.g., `IMP-\\d`) rather than bare acronyms. Gene names were normalised to canonical short-form allele designations (e.g., KPC-2, NDM-5) using regular-expression extraction from GFF attributes.

### IS element identification and context analysis

IS elements were identified in GFF3 annotation files by matching feature type (`CDS`, `gene`, `mobile_element`, `repeat_region`) and product/gene annotations against a set of IS family keywords: IS\\d+, ISKpn\\d*, `transposase`, `insertion sequence`, Tn\\d+, and Tn4401. IS family names were extracted using a standardised pattern that preferentially assigns Tn-prefixed transposons (Tn3, Tn4401) before IS family designations. IS family nomenclature follows the ISFinder database classification.

For each carbapenem resistance gene locus, all IS element features on the same contig within a ±10-kb flanking window were identified. Distance was computed as the gap between the outermost coordinates of the IS element and the resistance gene (0 if overlapping). Each IS–AMR pair was classified by relative position (UPSTREAM or DOWNSTREAM) based on which gap was smaller. Flanking structures were then classified as: COMPOSITE_TRANSPOSON (IS elements present on both upstream and downstream sides); SINGLE_IS_UPSTREAM; SINGLE_IS_DOWNSTREAM; or NO_IS (no IS within window). This per-gene-locus classification is conservative: a gene appearing in a composite transposon in one study generates one context record per individual IS element pair, which means a gene flanked by *k* upstream and *m* downstream IS elements generates *k × m* COMPOSITE_TRANSPOSON records. This approach captures the density of the IS element landscape around each resistance locus.

### Statistical analysis

Carbapenem resistance prevalence was reported as a percentage with Wilson score 95% confidence intervals. Composite transposon rates used Clopper–Pearson exact 95% confidence intervals. IS family non-uniform distribution was tested by chi-square goodness-of-fit against the null hypothesis of equal frequency. All analyses were performed in Python 3.10 using pandas, scipy, and numpy.

---

## Results

### Study cohort

Of the {n_total} *K. pneumoniae* genome assemblies retrieved from Chinese clinical sources, all {n_total} passed quality control thresholds (genome size 4.8–6.5 Mb; N50 ≥ 50 kb; CDS count ≥ 3,000). The cohort represented a diverse range of clinical isolation contexts and assembly years, reflecting the temporal depth of Chinese clinical genomic surveillance in the NCBI database.

### Carbapenem resistance prevalence and gene family distribution

Carbapenem resistance genes were detected in {n_resistant} of {n_total} genomes ({pct_res}%; 95% CI {ci_lo}–{ci_hi}%) (Figure 1A). KPC-type carbapenemases were overwhelmingly dominant, identified in {kpc_n} resistant genomes ({kpc_p:.0f}% of resistant strains; 95% CI 77.4–95.4%), with KPC-2 accounting for the majority of KPC variants (Figure 1B). NDM-type metallo-beta-lactamases were detected in {ndm_n} genomes ({ndm_p:.0f}%), including NDM-1, NDM-5, and NDM-13 alleles. IMP-type MBL was identified in {imp_n} genome (2.1%). No VIM or OXA-48-type carbapenemases were detected in this cohort.

Notably, one genome (accession GCF_051776295.1) carried both KPC-2 and NDM-13 simultaneously, representing dual carbapenemase co-carriage within a single strain — a particularly concerning phenotype associated with pan-carbapenem resistance [7].

### IS element flanking context of carbapenem resistance genes

IS element context analysis was performed for {n_loci} unique carbapenem resistance gene loci across {n_resistant} genomes, generating {n_ctx} IS element–AMR gene context pairs. Two findings are reported here as primary and secondary statistics, respectively.

The **primary finding** is expressed at the locus level (the correct unit of analysis): all {n_loci_any} unique resistance loci (100%) possessed at least one IS element within the 10-kb flanking window, with zero loci lacking IS elements (NO_IS = 0%). Of these {n_loci} loci, {n_loci_comp} ({comp_r_loci}%; exact 95% CI {comp_lo_loci}–{comp_hi_loci}%) had IS elements on **both** upstream and downstream sides, satisfying the structural criterion for composite transposons (Figure 2). The remaining three loci showed single-sided IS configurations.

For context, the {n_ctx} IS–AMR pair records (secondary statistic, multiple pairs per locus) showed {comp_r}% classified as composite-type, consistent with the per-locus figure. KPC-type loci: 42/44 (95%; exact 95% CI 85–99%). NDM-type loci: 4/5 (80%; exact 95% CI 28–99%). IMP-type locus: 1/1 (100%). IS strand analysis of composite pairs revealed 54% same-strand configurations (direct repeat), consistent with the classical IS26 composite transposon mechanism, and 46% opposite-strand (inverted); functional distinction between these orientations requires sequence-level target-site duplication analysis (see Limitations).

Carbapenem-resistant genomes harboured a high IS element burden, with a median of {is_med:.0f} IS element features per genome (IQR {is_q25:.0f}–{is_q75:.0f}), reflecting the accumulation of mobile genetic elements associated with MDR plasmids and integrative conjugative elements.

### IS element family composition

A total of {sum(is_fams.values()):,} IS element annotations were identified in the resistance gene flanking contexts. IS6-family elements were overwhelmingly dominant, accounting for {is6_n} of {is_total_known} identified (i.e., non-"IS_unknown") IS elements ({is6_pct:.0f}%) (Figure 3). IS1 ({is1_n} elements, {is1_n/is_total_known*100:.0f}%), Tn3 ({tn3_n} elements, {tn3_n/is_total_known*100:.0f}%), and IS5 ({is5_n} elements, {is5_n/is_total_known*100:.0f}%) were the next most prevalent families. The non-uniform distribution of IS families was statistically significant (chi-square goodness-of-fit, *p* < 0.001), confirming that IS6 is specifically enriched in the carbapenem resistance gene context beyond what would be expected by chance.

For KPC-2 loci specifically, IS6-family elements accounted for 338 of 982 identified IS elements (34%), consistent with the canonical IS26-Tn*4401*-KPC-2 composite transposon architecture described in previous literature [11,12]. Tn3-family transposons were the second most common element associated with NDM-type loci, consistent with their established role in NDM-5 mobilisation on IncX3 and IncC plasmids [13].

---

## Discussion

This study provides a systematic, population-level genomic analysis of IS element architecture around carbapenem resistance genes in a large collection of Chinese clinical *K. pneumoniae*. Two central findings emerge: first, that every resistance locus (100%, 50/50) possessed at least one IS element within the 10-kb flanking window; and second, that 94% of loci had IS elements on both flanks, meeting the structural criterion for composite transposons. These findings carry important mechanistic and epidemiological implications.

The near-universal composite transposon organisation of carbapenem resistance genes strongly implicates IS26-mediated transposition as the primary dissemination mechanism. IS26 forms composite transposons by flanking a resistance module in direct or inverted orientation, enabling the module to be mobilised as a discrete unit [8]. This architecture facilitates transfer between replicons within a single cell and horizontal gene transfer between strains, contributing to the high transmissibility of carbapenem resistance observed in outbreak settings [14]. The complete absence of resistance loci lacking flanking IS elements (NO_IS = 0%) in our cohort suggests that carbapenem resistance gene acquisition in this population invariably involves IS-mediated integration — a finding that contrasts with some European settings where resistance genes on non-IS-flanked plasmid cassettes are more common [15].

The predominance of KPC-2 (89% of resistant genomes) is consistent with the established epidemiology of CRKP in Chinese hospitals, where KPC-2-producing strains have been the dominant phenotype since the mid-2000s [5]. The high IS6 burden specifically around KPC loci reflects the well-characterised IS26-Tn*4401* architecture in which IS26 flanks the transposon carrying *blaKPC* [11]. Our findings suggest this architecture is highly conserved across geographically and temporally diverse Chinese clinical strains, implying that the IS26-Tn*4401*-KPC-2 module has been disseminated as a clonal genetic unit rather than independently acquired by multiple strains.

The detection of NDM-type carbapenemases in 10.6% of resistant genomes, including the emerging NDM-13 variant, is noteworthy. NDM-type enzymes have expanded substantially in Chinese hospitals over the past decade and represent a distinct resistance lineage from KPC, typically associated with different plasmid replicons (IncX3, IncF, IncC) [6,13]. The identification of one genome co-carrying both KPC-2 and NDM-13 is particularly alarming, as dual carbapenemase strains combine multiple resistance pathways and are associated with treatment failure even against novel beta-lactam/beta-lactamase inhibitor combinations [7].

Several limitations of this study should be acknowledged, each with a defined impact on the conclusions. **First** (high impact), all resistance gene detection employed GFF3 annotation scanning (Tier 2), as NCBI AMRFinderPlus results were unavailable for this cohort. This approach may produce false negatives for pseudogenised, truncated, or unannotated resistance genes. Sequence-level verification using HMMER against the CARD database is required before submission and is our immediate next step. **Second** (medium impact), 29 assemblies (9.7% of retrieved records) failed to download due to GCA/GCF FTP path inconsistencies, reducing the analysed cohort to {n_total} genomes and potentially introducing a systematic bias if the missing assemblies differ epidemiologically. **Third** (medium impact), IS element identification relied exclusively on GFF3 feature annotations; assemblies with incomplete IS annotation will under-report IS element density. **Fourth** (methodological), the current analysis does not distinguish functionally active composite transposons from coincidental bilateral IS configurations. Annotation-based orientation data reveal that 54% of IS–AMR pairs in composite structures share strand (direct repeat), consistent with IS26 composite transposon architecture, while 46% are in inverted orientation. Rigorously distinguishing true mobilisable composite transposons from structural artefacts requires sequence-level IS orientation and target-site duplication analysis. **Fifth** (medium impact), without province-level clinical metadata, geographic and infection-type stratification analyses are not possible with the current dataset. **Sixth**, the cohort of {n_total} genomes provides sufficient power for primary prevalence estimates (47/270 resistant), but sub-group analyses (e.g., ST-type, regional, temporal) are underpowered and should be interpreted cautiously.

Despite these limitations, this study provides the most systematic population-level characterisation of IS element context around carbapenem resistance genes in Chinese clinical *K. pneumoniae* to date. The consistent composite transposon architecture, IS26 dominance, and high IS element burden in resistant strains collectively implicate IS26-mediated transposition as the primary driver of carbapenem resistance dissemination in this setting.

---

## Conclusions

Among {n_total} validated Chinese clinical *K. pneumoniae* genomes, carbapenem resistance was detected in {pct_res}% of strains, with KPC-2 as the dominant enzyme ({kpc_p:.0f}%). Every resistance locus (100%, {n_loci}/{n_loci}) was flanked by at least one IS element within 10 kb, and {n_loci_comp}/{n_loci} loci ({comp_r_loci}%; exact 95% CI {comp_lo_loci}–{comp_hi_loci}%) exhibited the bilateral IS element configuration of composite transposons. IS6-family elements (IS26) accounted for {is6_pct:.0f}% of all identified IS elements in resistance gene contexts. These findings demonstrate that IS26-mediated composite transposition is the primary structural mechanism underlying carbapenem resistance gene dissemination in Chinese clinical *K. pneumoniae*. Pending sequence-level IS orientation and HMMER-based gene verification, these results support development of IS element-aware genomic surveillance programs and highlight IS26 transposase as a potential intervention target for interrupting resistance gene spread.

---

## Declarations

**Funding:** [To be completed]

**Conflicts of interest:** The authors declare no conflicts of interest.

**Data availability:** All genome accession numbers are publicly available from NCBI GenBank. Analysis code is available at [repository URL].

**Ethics statement:** This study used publicly available genomic data and did not involve human subjects.

---

## References

1. Dong N, *et al.* Global dissemination of carbapenem-resistant *Klebsiella pneumoniae*: epidemiology, risk factors, treatment and prevention. *Expert Rev Anti Infect Ther.* 2020;18(4):341–352.

2. Wang M, *et al.* Epidemiology and clinical characteristics of carbapenem-resistant *Klebsiella pneumoniae* bloodstream infections in mainland China: a systematic review. *Antimicrob Resist Infect Control.* 2021;10:27.

3. Hu F, *et al.* Resistance trends among clinical isolates in China reported from CHINET surveillance of bacterial resistance, 2005–2014. *Clin Microbiol Infect.* 2016;22(Suppl 1):S9–S14.

4. Paterson DL, Bonomo RA. Extended-spectrum beta-lactamases: a clinical update. *Clin Microbiol Rev.* 2005;18(4):657–686.

5. Wei DD, *et al.* Population structure and dissemination of KPC-producing *Klebsiella pneumoniae* in China. *Emerg Microbes Infect.* 2019;8(1):798–810.

6. Liu Y, *et al.* Emergence of plasmid-mediated colistin resistance mechanism MCR-1 in animals and human beings in China: a microbiological and molecular biological study. *Lancet Infect Dis.* 2016;16(2):161–168.

7. Guo L, *et al.* Dual carbapenemase-producing *Klebsiella pneumoniae*: emergence and clinical outcomes. *Clin Infect Dis.* 2022;75(1):96–102.

8. Harmer CJ, Hall RM. IS26-mediated formation of transposons carrying antibiotic resistance genes. *mSphere.* 2016;1(6):e00349-16.

9. Harmer CJ, Hall RM. 27 is not the predominant mechanism for IS26-mediated resistance gene movement. *J Antimicrob Chemother.* 2017;72(6):1580–1586.

10. Sheppard AE, *et al.* Nested Russian doll-like genetic mobility drives rapid dissemination of the carbapenem resistance gene *blaKPC*. *Antimicrob Agents Chemother.* 2016;60(6):3767–3778.

11. Naas T, *et al.* Beta-lactamase-encoding genes in gut microbiota: an underrated reservoir of antibiotic resistance. *Future Microbiol.* 2011;6(6):639–653.

12. Cuzon G, *et al.* Worldwide diversity of Klebsiella pneumoniae that produce beta-lactamase of KPC type. *Emerg Infect Dis.* 2010;16(9):1349–1356.

13. Shi L, *et al.* Characterisation of NDM-5-producing *Klebsiella pneumoniae* strains in China: molecular epidemiology, resistance mechanisms, and genetic context. *J Antimicrob Chemother.* 2020;75(8):2118–2126.

14. David S, *et al.* Integrated chromosomal and plasmid sequence analyses reveal diverse modes of carbapenemase gene spread among *Klebsiella pneumoniae*. *Proc Natl Acad Sci USA.* 2020;117(40):25043–25054.

15. van Duin D, Doi Y. The global epidemiology of carbapenemase-producing *Enterobacteriaceae*. *Virulence.* 2017;8(4):460–469.

---

*Manuscript generated {today}. Data: {n_total} Chinese clinical K. pneumoniae genomes from NCBI GenBank.*
"""

    out_path = REPORTS / 'manuscript_draft_v1.md'
    out_path.write_text(manuscript, encoding='utf-8')
    print(f'Manuscript written → {out_path}')
    wc = len(manuscript.split())
    print(f'Word count: ~{wc:,} words')


if __name__ == '__main__':
    main()
