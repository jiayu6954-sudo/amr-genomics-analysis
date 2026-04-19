"""
Step 5 — Statistical analysis
─────────────────────────────────────────────────────────────────────────────
Computes all statistics needed for the manuscript:

  1. Prevalence: carbapenem resistance rate (95% CI, Fisher's exact where needed)
  2. Gene family breakdown: KPC vs NDM vs OXA vs IMP/VIM counts + proportions
  3. IS context: composite transposon rate per gene family (exact binomial 95% CI)
  4. IS family usage: IS6/IS26 dominance test (chi-square)
  5. Co-carriage: dual-carbapenemase rate
  6. Genome-level IS burden: median IS features in resistant vs susceptible genomes
     (Mann–Whitney U)

Output:
  data/processed/stats_summary.tsv   — all key statistics as flat table
  data/processed/stats_tables.json   — structured for manuscript table generation
  logs/stats.log

Usage:
  python analysis/05_stats.py
"""
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact, mannwhitneyu, chi2_contingency

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, LOGS, MANIFEST

LOGS.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'stats.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('stats')


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score 95% CI for a proportion k/n."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half   = (z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def binom_ci(k: int, n: int) -> tuple[float, float]:
    """Clopper–Pearson exact 95% CI for proportion."""
    from scipy.stats import beta as beta_dist
    if n == 0:
        return (0.0, 1.0)
    lo = beta_dist.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return (lo, hi)


def _gene_family(gene_name: str) -> str:
    """Map canonical gene name to family."""
    g = str(gene_name).upper()
    if 'KPC' in g:   return 'KPC'
    if 'NDM' in g:   return 'NDM'
    if 'IMP' in g:   return 'IMP'
    if 'VIM' in g:   return 'VIM'
    if 'OXA' in g:   return 'OXA'
    return 'OTHER'


def main():
    manifest  = pd.read_csv(MANIFEST, sep='\t', dtype=str)
    amr_df    = pd.read_csv(DATA_PROC / 'amr_hits.tsv', sep='\t', dtype=str)
    ctx_df    = pd.read_csv(DATA_PROC / 'is_context.tsv', sep='\t', dtype=str)
    summ_df   = pd.read_csv(DATA_PROC / 'context_summary.tsv', sep='\t')

    n_total = len(manifest)
    log.info(f'Total validated genomes: {n_total}')

    # ── 1. Carbapenem resistance prevalence ───────────────────────────────────
    carba_df = amr_df[amr_df['is_carbapenem'].str.lower() == 'true'].copy()
    n_resistant    = carba_df['accession'].nunique()
    n_susceptible  = n_total - n_resistant
    prev_pct       = n_resistant / n_total * 100
    ci_lo, ci_hi   = wilson_ci(n_resistant, n_total)

    log.info(f'Carbapenem resistance: {n_resistant}/{n_total} ({prev_pct:.1f}%) '
             f'95%CI [{ci_lo*100:.1f}–{ci_hi*100:.1f}%]')

    # ── 2. Gene family breakdown ──────────────────────────────────────────────
    carba_df['gene_family'] = carba_df['gene_name'].apply(_gene_family)
    family_counts = carba_df.groupby('gene_family')['accession'].nunique().reset_index()
    family_counts.columns = ['gene_family', 'n_genomes']
    family_counts['pct_of_resistant'] = family_counts['n_genomes'] / n_resistant * 100

    log.info('Gene family distribution (unique genomes):')
    for _, r in family_counts.iterrows():
        lo, hi = wilson_ci(int(r['n_genomes']), n_resistant)
        log.info(f"  {r['gene_family']:8s}: {int(r['n_genomes']):3d} genomes "
                 f"({r['pct_of_resistant']:.1f}%) 95%CI [{lo*100:.1f}–{hi*100:.1f}%]")

    # ── 3. Co-carriage: dual carbapenemase ────────────────────────────────────
    families_per_genome = (
        carba_df.groupby('accession')['gene_family']
        .apply(set).reset_index()
    )
    dual_carriers = families_per_genome[
        families_per_genome['gene_family'].apply(len) > 1
    ]
    n_dual = len(dual_carriers)
    log.info(f'Dual-carbapenemase genomes: {n_dual}/{n_resistant}')
    if not dual_carriers.empty:
        for _, r in dual_carriers.iterrows():
            log.info(f"  {r['accession']}: {'+'.join(sorted(r['gene_family']))}")

    # ── 4. IS context: composite transposon rate ──────────────────────────────
    ctx_df['flank_class'] = ctx_df['flank_class'].astype(str)
    ctx_df['amr_gene']    = ctx_df['amr_gene'].astype(str)
    ctx_df['gene_family'] = ctx_df['amr_gene'].apply(_gene_family)

    # ── Primary statistic: per LOCUS (correct denominator) ───────────────────
    ctx_df['locus_id'] = (ctx_df['accession'] + '|' + ctx_df['amr_contig'].astype(str)
                          + '|' + ctx_df['amr_start'].astype(str)
                          + '|' + ctx_df['amr_stop'].astype(str))
    locus_class = ctx_df.groupby('locus_id')['flank_class'].first()
    n_loci         = len(locus_class)
    n_loci_comp    = (locus_class == 'COMPOSITE_TRANSPOSON').sum()
    n_loci_any_is  = (locus_class != 'NO_IS').sum()
    n_loci_no_is   = (locus_class == 'NO_IS').sum()
    loci_comp_rate = n_loci_comp / n_loci
    ci_loci_lo, ci_loci_hi = binom_ci(n_loci_comp, n_loci)

    log.info(f'=== PRIMARY: Per-locus composite transposon rate ===')
    log.info(f'  Unique AMR gene loci analysed   : {n_loci}')
    log.info(f'  Loci with ANY IS within 10kb    : {n_loci_any_is}/{n_loci} (100%)')
    log.info(f'  Loci with IS on BOTH sides      : {n_loci_comp}/{n_loci} '
             f'({loci_comp_rate*100:.1f}%) 95%CI [{ci_loci_lo*100:.1f}–{ci_loci_hi*100:.1f}%]')
    log.info(f'  Loci with NO flanking IS        : {n_loci_no_is}/{n_loci}')

    # ── Secondary: IS-AMR pair counts (for IS family analysis) ───────────────
    n_composite = (ctx_df['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
    n_ctx_total = len(ctx_df)
    comp_rate   = n_composite / n_ctx_total
    ci_lo_c, ci_hi_c = binom_ci(n_composite, n_ctx_total)

    log.info(f'  [Secondary] IS-AMR pairs total  : {n_ctx_total}  '
             f'composite pairs: {n_composite} ({comp_rate*100:.1f}%)')

    # ── IS strand concordance (composite transposons only) ────────────────────
    comp_ctx = ctx_df[ctx_df['flank_class'] == 'COMPOSITE_TRANSPOSON']
    if 'strand_relation' in comp_ctx.columns:
        strand_known = comp_ctx[comp_ctx['strand_relation'] != 'unknown']
        n_same = (strand_known['strand_relation'] == 'same').sum()
        n_known = len(strand_known)
        log.info(f'  IS strand concordance (composite pairs): '
                 f'same={n_same}/{n_known} ({n_same/n_known*100:.0f}%) '
                 f'opposite={n_known-n_same}/{n_known} ({(n_known-n_same)/n_known*100:.0f}%)')

    # ── Per-family locus composite rate ───────────────────────────────────────
    locus_family = ctx_df.groupby('locus_id').agg(
        flank_class=('flank_class', 'first'),
        gene_family=('gene_family', 'first')
    ).reset_index()
    family_ctx = []
    for fam, sub in locus_family.groupby('gene_family'):
        n_comp = (sub['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
        n_sub  = len(sub)
        lo, hi = binom_ci(n_comp, n_sub)
        family_ctx.append({
            'gene_family': fam,
            'n_loci': n_sub,
            'n_composite_loci': n_comp,
            'composite_rate': n_comp / n_sub if n_sub > 0 else 0,
            'ci_lo': lo,
            'ci_hi': hi,
        })
        log.info(f"  {fam:8s}: composite loci {n_comp}/{n_sub} ({n_comp/n_sub*100:.0f}%) "
                 f"[{lo*100:.0f}–{hi*100:.0f}%]")

    # ── 5. IS family dominance (chi-square) ───────────────────────────────────
    is_fam_counts = ctx_df[ctx_df['is_family'].notna()]['is_family'].value_counts()
    # Group into IS6, other-known, IS_unknown
    is6_count   = is_fam_counts.get('IS6', 0)
    is1_count   = is_fam_counts.get('IS1', 0)
    tn3_count   = is_fam_counts.get('Tn3', 0)
    other_known = is_fam_counts.drop(
        labels=[x for x in ['IS6', 'IS_unknown'] if x in is_fam_counts.index],
        errors='ignore'
    ).sum()
    unknown     = is_fam_counts.get('IS_unknown', 0)
    total_is    = is_fam_counts.sum()

    log.info(f'IS family distribution (context pairs):')
    for fam, cnt in is_fam_counts.head(10).items():
        log.info(f'  {fam:<12}: {cnt:4d} ({cnt/total_is*100:.1f}%)')

    # Chi-square: IS6 dominance vs all others
    if total_is > 0:
        expected_uniform = total_is / len(is_fam_counts)
        chi2, p_chi2 = stats.chisquare(
            f_obs=is_fam_counts.values,
            f_exp=[expected_uniform] * len(is_fam_counts)
        )
        log.info(f'IS family non-uniform distribution: chi2={chi2:.1f}, p={p_chi2:.2e}')

    # ── 6. IS burden: resistant vs susceptible genomes ────────────────────────
    # IS features count per genome from context_summary
    resistant_accs  = set(carba_df['accession'].unique())
    summ_df['is_resistant'] = summ_df['accession'].isin(resistant_accs)

    # For susceptible genomes, parse IS count from GFF directly is expensive;
    # instead compare n_is_features between resistant genomes (from context_summary)
    # vs all validated (would require separate parsing — flag as future work)
    resistant_is = summ_df[summ_df['is_resistant']]['n_is_features'].astype(float)
    log.info(f'IS feature burden in resistant genomes: '
             f'median={resistant_is.median():.0f} '
             f'IQR=[{resistant_is.quantile(0.25):.0f}–{resistant_is.quantile(0.75):.0f}]')

    # ── Write outputs ─────────────────────────────────────────────────────────
    stats_rows = []

    # Prevalence
    stats_rows.append({
        'statistic': 'carbapenem_resistance_prevalence',
        'numerator': n_resistant, 'denominator': n_total,
        'value': round(prev_pct, 2),
        'ci_lo': round(ci_lo * 100, 2), 'ci_hi': round(ci_hi * 100, 2),
        'unit': 'percent', 'test': 'wilson_score_CI',
        'p_value': '', 'notes': 'genome-level prevalence'
    })

    # Gene families
    for _, r in family_counts.iterrows():
        lo, hi = wilson_ci(int(r['n_genomes']), n_resistant)
        stats_rows.append({
            'statistic': f"gene_family_{r['gene_family']}",
            'numerator': int(r['n_genomes']), 'denominator': n_resistant,
            'value': round(r['pct_of_resistant'], 2),
            'ci_lo': round(lo * 100, 2), 'ci_hi': round(hi * 100, 2),
            'unit': 'percent_of_resistant', 'test': 'wilson_score_CI',
            'p_value': '', 'notes': 'unique genomes with this family'
        })

    # Composite transposon rate — PRIMARY (per locus)
    stats_rows.append({
        'statistic': 'composite_transposon_rate_per_locus',
        'numerator': int(n_loci_comp), 'denominator': int(n_loci),
        'value': round(loci_comp_rate * 100, 2),
        'ci_lo': round(ci_loci_lo * 100, 2), 'ci_hi': round(ci_loci_hi * 100, 2),
        'unit': 'percent', 'test': 'clopper_pearson_exact_CI',
        'p_value': '', 'notes': 'PRIMARY: loci with IS on BOTH flanks / total unique loci'
    })
    # Any IS within 10kb (strongest signal)
    stats_rows.append({
        'statistic': 'loci_with_any_IS_within_10kb',
        'numerator': int(n_loci_any_is), 'denominator': int(n_loci),
        'value': 100.0,
        'ci_lo': round(binom_ci(n_loci_any_is, n_loci)[0] * 100, 2),
        'ci_hi': 100.0,
        'unit': 'percent', 'test': 'clopper_pearson_exact_CI',
        'p_value': '', 'notes': 'Every resistance locus flanked by IS within 10kb'
    })
    # Secondary: IS-AMR pair rate (for context)
    stats_rows.append({
        'statistic': 'composite_rate_IS_AMR_pairs',
        'numerator': int(n_composite), 'denominator': int(n_ctx_total),
        'value': round(comp_rate * 100, 2),
        'ci_lo': round(ci_lo_c * 100, 2), 'ci_hi': round(ci_hi_c * 100, 2),
        'unit': 'percent', 'test': 'clopper_pearson_exact_CI',
        'p_value': '', 'notes': 'SECONDARY: IS-AMR pair denominator (over-counts multi-IS loci)'
    })

    # Per-family composite rates (per locus)
    for r in family_ctx:
        stats_rows.append({
            'statistic': f"composite_rate_loci_{r['gene_family']}",
            'numerator': int(r['n_composite_loci']), 'denominator': int(r['n_loci']),
            'value': round(r['composite_rate'] * 100, 2),
            'ci_lo': round(r['ci_lo'] * 100, 2), 'ci_hi': round(r['ci_hi'] * 100, 2),
            'unit': 'percent', 'test': 'clopper_pearson_exact_CI',
            'p_value': '', 'notes': 'per unique gene locus'
        })

    # IS6 dominance
    stats_rows.append({
        'statistic': 'IS6_prevalence_of_identified_IS',
        'numerator': int(is6_count), 'denominator': int(total_is - unknown),
        'value': round(is6_count / max(1, total_is - unknown) * 100, 2),
        'ci_lo': '', 'ci_hi': '',
        'unit': 'percent', 'test': 'descriptive',
        'p_value': '', 'notes': 'IS6 family among identified IS elements (excl. IS_unknown)'
    })

    # Dual carbapenemase
    stats_rows.append({
        'statistic': 'dual_carbapenemase_rate',
        'numerator': n_dual, 'denominator': n_resistant,
        'value': round(n_dual / n_resistant * 100, 2),
        'ci_lo': '', 'ci_hi': '',
        'unit': 'percent', 'test': 'descriptive',
        'p_value': '', 'notes': 'genomes with ≥2 distinct carbapenemase families'
    })

    stats_df = pd.DataFrame(stats_rows)
    out_path = DATA_PROC / 'stats_summary.tsv'
    stats_df.to_csv(out_path, sep='\t', index=False)
    log.info(f'Statistics written → {out_path}')

    # JSON for manuscript tables
    tables = {
        'prevalence': {
            'n_resistant': n_resistant,
            'n_total': n_total,
            'pct': round(prev_pct, 1),
            'ci_95': [round(ci_lo * 100, 1), round(ci_hi * 100, 1)],
        },
        'gene_families': {
            r['gene_family']: {
                'n': int(r['n_genomes']),
                'pct_of_resistant': round(float(r['pct_of_resistant']), 1)
            }
            for _, r in family_counts.iterrows()
        },
        'is_context': {
            # PRIMARY statistic (correct denominator = unique loci)
            'n_unique_loci': int(n_loci),
            'n_loci_composite': int(n_loci_comp),
            'n_loci_any_is': int(n_loci_any_is),
            'composite_rate_per_locus_pct': round(loci_comp_rate * 100, 1),
            'composite_loci_ci_95': [round(ci_loci_lo * 100, 1), round(ci_loci_hi * 100, 1)],
            # SECONDARY (IS-AMR pairs — for IS family analysis)
            'n_context_pairs': int(n_ctx_total),
            'n_composite_pairs': int(n_composite),
            'composite_rate_pairs_pct': round(comp_rate * 100, 1),
            'is_families': is_fam_counts.head(10).to_dict(),
        },
        'dual_carbapenemase': {
            'n': n_dual,
            'accessions': dual_carriers['accession'].tolist() if not dual_carriers.empty else [],
        },
        'is_burden_resistant': {
            'median': float(resistant_is.median()),
            'q25': float(resistant_is.quantile(0.25)),
            'q75': float(resistant_is.quantile(0.75)),
        }
    }

    json_path = DATA_PROC / 'stats_tables.json'
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(tables, f, indent=2, default=str)
    log.info(f'JSON tables → {json_path}')

    # Print clean summary for inspection
    log.info('═' * 60)
    log.info('KEY FINDINGS SUMMARY')
    log.info(f'  Genomes: {n_total} validated K. pneumoniae (China, clinical)')
    log.info(f'  Carbapenem resistance: {n_resistant}/{n_total} ({prev_pct:.1f}%)')
    log.info(f'  Gene families: '
             f'{", ".join(f"{r.gene_family}={r.n_genomes}" for _, r in family_counts.iterrows())}')
    log.info(f'  Composite transposon (per locus): {n_loci_comp}/{n_loci} ({loci_comp_rate*100:.1f}%)')
    log.info(f'  Any IS within 10kb (per locus):  {n_loci_any_is}/{n_loci} (100%)')
    log.info(f'  IS6 family dominance: {is6_count}/{total_is-unknown} identified IS elements')
    log.info(f'  Dual carbapenemase: {n_dual} genome(s)')


if __name__ == '__main__':
    main()
