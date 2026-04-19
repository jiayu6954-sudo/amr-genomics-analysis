"""
Step 11 — Temporal Subgroup & Sensitivity Analysis
─────────────────────────────────────────────────────────────────────────────
Addresses NCBI submission bias identified in P1 (492-genome full cohort).

Analyses:
  A. Year-stratified resistance rates and IS6 burden
  B. 2025 cohort sensitivity analysis (n=188, presumed most representative)
  C. Pre-2022 vs 2022+ cohort comparison (biased vs less-biased)
  D. Composite transposon rate stability across year groups

Output:
  data/processed/subgroup_stats.json
  data/processed/subgroup_table.tsv
"""

import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, LOGS

LOGS.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'subgroup.log', encoding='utf-8'),
        logging.StreamHandler(
            open(sys.stdout.fileno(), mode='w', encoding='utf-8', closefd=False)
        ),
    ]
)
log = logging.getLogger('subgroup')


# ── stats helpers ─────────────────────────────────────────────────────────────

def cliffs_delta(a, b):
    if len(a) == 0 or len(b) == 0:
        return float('nan')
    mat = np.sign(np.array(a)[:, None] - np.array(b)[None, :])
    return float(mat.mean())


def auc_from_mw(a, b):
    """AUC = P(resistant > susceptible). U_greater / (n_r * n_s)."""
    if len(a) == 0 or len(b) == 0:
        return float('nan')
    U, _ = mannwhitneyu(a, b, alternative='greater')
    return float(U / (len(a) * len(b)))


def cohort_stats(df_sub: pd.DataFrame, label: str) -> dict:
    """Compute key stats for a sub-cohort."""
    n = len(df_sub)
    if n == 0:
        return {'label': label, 'n': 0}

    res = df_sub[df_sub['is_resistant']]
    sus = df_sub[~df_sub['is_resistant']]
    n_res, n_sus = len(res), len(sus)

    r_is6 = res['n_is6'].values
    s_is6 = sus['n_is6'].values

    result = {
        'label':          label,
        'n_total':        n,
        'n_resistant':    n_res,
        'n_susceptible':  n_sus,
        'resistance_pct': round(100 * n_res / n, 1) if n > 0 else None,
        'is6_resistant_median':   float(np.median(r_is6)) if n_res > 0 else None,
        'is6_resistant_iqr':      [float(np.percentile(r_is6, 25)),
                                   float(np.percentile(r_is6, 75))] if n_res > 0 else None,
        'is6_susceptible_median': float(np.median(s_is6)) if n_sus > 0 else None,
        'is6_susceptible_iqr':    [float(np.percentile(s_is6, 25)),
                                   float(np.percentile(s_is6, 75))] if n_sus > 0 else None,
        'cliffs_delta':   round(cliffs_delta(r_is6, s_is6), 4) if (n_res > 0 and n_sus > 0) else None,
        'auc':            round(auc_from_mw(r_is6, s_is6), 4) if (n_res > 0 and n_sus > 0) else None,
    }

    if n_res > 0 and n_sus > 0:
        _, p = mannwhitneyu(r_is6, s_is6, alternative='greater')
        result['mannwhitney_p'] = float(p)
    else:
        result['mannwhitney_p'] = None

    return result


def main():
    log.info('=' * 70)
    log.info('TEMPORAL SUBGROUP & SENSITIVITY ANALYSIS')
    log.info('=' * 70)

    # ── Load data ─────────────────────────────────────────────────────────────
    burden_path = DATA_PROC / 'is_burden_all.tsv'
    is_ctx_path = DATA_PROC / 'is_context.tsv'
    amr_path    = DATA_PROC / 'amr_hits.tsv'

    df = pd.read_csv(burden_path, sep='\t')
    ctx = pd.read_csv(is_ctx_path, sep='\t')
    amr = pd.read_csv(amr_path, sep='\t')

    log.info(f'Loaded {len(df)} genomes from is_burden_all.tsv')
    log.info(f'Year range: {df["assembly_year"].min()} – {df["assembly_year"].max()}')

    # ── A. Year-stratified table ──────────────────────────────────────────────
    log.info('')
    log.info('A. YEAR-STRATIFIED ANALYSIS')
    log.info('-' * 50)

    year_rows = []
    for yr in sorted(df['assembly_year'].unique()):
        sub = df[df['assembly_year'] == yr]
        st  = cohort_stats(sub, str(yr))
        year_rows.append(st)
        log.info(
            f'  {yr}: n={st["n_total"]:>3}  '
            f'resistant={st["n_resistant"]:>3} ({st["resistance_pct"]:>5.1f}%)  '
            f'IS6_res_med={st["is6_resistant_median"]}  '
            f'IS6_sus_med={st["is6_susceptible_median"]}  '
            f'AUC={st["auc"]}'
        )

    # ── B. Four cohort comparison ─────────────────────────────────────────────
    log.info('')
    log.info('B. COHORT COMPARISON')
    log.info('-' * 50)

    cohorts = {
        'Full (492)':         df,
        'Pre-2022 (biased)':  df[df['assembly_year'] < 2022],
        '2022+ (mixed)':      df[df['assembly_year'] >= 2022],
        '2025 only':          df[df['assembly_year'] == 2025],
        '2025-2026 (recent)': df[df['assembly_year'] >= 2025],
    }

    cohort_results = {}
    for label, sub in cohorts.items():
        st = cohort_stats(sub, label)
        cohort_results[label] = st
        log.info(
            f'  {label:<25}  n={st["n_total"]:>3}  '
            f'resist={st["resistance_pct"]:>5.1f}%  '
            f'IS6_res_med={st["is6_resistant_median"]}  '
            f'IS6_sus_med={st["is6_susceptible_median"]}  '
            f'delta={st["cliffs_delta"]}  AUC={st["auc"]}'
        )

    # ── C. Composite transposon rate by cohort ────────────────────────────────
    log.info('')
    log.info('C. COMPOSITE TRANSPOSON RATE BY COHORT')
    log.info('-' * 50)

    # Map accession → year
    acc_year = df.set_index('accession')['assembly_year'].to_dict()
    ctx['year'] = ctx['accession'].map(acc_year)

    comp_results = {}
    for label, yrs in [
        ('Pre-2022', lambda y: y < 2022),
        ('2022+',    lambda y: y >= 2022),
        ('2025 only',lambda y: y == 2025),
        ('All',      lambda y: True),
    ]:
        sub = ctx[ctx['year'].apply(yrs)]
        n_pairs = len(sub)
        if n_pairs == 0:
            continue
        n_comp = (sub['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
        pct = 100 * n_comp / n_pairs
        comp_results[label] = {'n_pairs': n_pairs, 'n_composite': int(n_comp), 'pct': round(pct, 1)}
        log.info(f'  {label:<15}  pairs={n_pairs:>4}  composite={n_comp} ({pct:.1f}%)')

    # ── D. Gene family breakdown by era ──────────────────────────────────────
    log.info('')
    log.info('D. AMR GENE FAMILY BREAKDOWN BY ERA')
    log.info('-' * 50)

    amr['year'] = amr['accession'].map(acc_year)
    for label, mask in [
        ('Pre-2022', amr['year'] < 2022),
        ('2022+',    amr['year'] >= 2022),
        ('2025',     amr['year'] == 2025),
    ]:
        sub = amr[mask & amr['is_carbapenem']]
        n = len(sub)
        if n == 0:
            continue
        dist = sub['gene_name'].value_counts()
        top = ', '.join(f'{g}={c}({100*c/n:.0f}%)' for g, c in dist.head(5).items())
        log.info(f'  {label:<10}  n={n}  {top}')

    # ── E. 2025 cohort — detailed sensitivity analysis ────────────────────────
    log.info('')
    log.info('E. 2025 COHORT — KEY SENSITIVITY ANALYSIS')
    log.info('=' * 50)

    df25 = df[df['assembly_year'] == 2025]
    st25 = cohort_stats(df25, '2025')

    log.info(f'  n_total:        {st25["n_total"]}')
    log.info(f'  Resistant:      {st25["n_resistant"]} ({st25["resistance_pct"]}%)')
    log.info(f'  Susceptible:    {st25["n_susceptible"]}')
    log.info(f'  IS6 res median: {st25["is6_resistant_median"]} IQR {st25["is6_resistant_iqr"]}')
    log.info(f'  IS6 sus median: {st25["is6_susceptible_median"]} IQR {st25["is6_susceptible_iqr"]}')
    log.info(f'  Cliff\'s delta:  {st25["cliffs_delta"]}')
    log.info(f'  AUC:            {st25["auc"]}')
    log.info(f'  p-value:        {st25["mannwhitney_p"]:.2e}' if st25['mannwhitney_p'] else '  p: N/A')

    # ── F. Bias test: Fisher exact on submission era vs resistance ────────────
    log.info('')
    log.info('F. SUBMISSION BIAS QUANTIFICATION')
    log.info('-' * 50)

    pre22  = df[df['assembly_year'] < 2022]
    post22 = df[df['assembly_year'] >= 2022]
    table = [
        [pre22['is_resistant'].sum(),  (~pre22['is_resistant']).sum()],
        [post22['is_resistant'].sum(), (~post22['is_resistant']).sum()],
    ]
    OR, p_fisher = fisher_exact(table)
    log.info(f'  Pre-2022:  {table[0][0]} resistant / {table[0][1]} susceptible')
    log.info(f'  2022+:     {table[1][0]} resistant / {table[1][1]} susceptible')
    log.info(f'  Fisher OR={OR:.2f}  p={p_fisher:.2e}')
    log.info(f'  Interpretation: Pre-2022 submissions are {OR:.1f}x more likely to')
    log.info(f'  be resistant than 2022+ submissions (p={p_fisher:.2e})')

    # ── Save outputs ──────────────────────────────────────────────────────────
    out_json = DATA_PROC / 'subgroup_stats.json'
    out_tsv  = DATA_PROC / 'subgroup_table.tsv'

    summary = {
        'year_stratified':    year_rows,
        'cohort_comparison':  cohort_results,
        'composite_by_era':   comp_results,
        'submission_bias': {
            'pre2022_n': int(len(pre22)),
            'pre2022_resistance_pct': round(100 * pre22['is_resistant'].mean(), 1),
            'post2022_n': int(len(post22)),
            'post2022_resistance_pct': round(100 * post22['is_resistant'].mean(), 1),
            'fisher_OR': round(OR, 2),
            'fisher_p': float(p_fisher),
        },
    }

    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    pd.DataFrame(year_rows).to_csv(out_tsv, sep='\t', index=False)

    # ── Final summary ─────────────────────────────────────────────────────────
    log.info('')
    log.info('=' * 70)
    log.info('SENSITIVITY ANALYSIS SUMMARY')
    log.info('=' * 70)
    log.info('')
    log.info('  Composite transposon rate — ROBUST across all eras:')
    for lbl, v in comp_results.items():
        log.info(f'    {lbl:<15}  {v["pct"]}%')
    log.info('')
    log.info('  IS6 AUC — COHORT-SENSITIVE (evidence of submission bias):')
    for lbl, v in cohort_results.items():
        log.info(f'    {lbl:<25}  resist={v["resistance_pct"]}%  AUC={v["auc"]}  delta={v["cliffs_delta"]}')
    log.info('')
    log.info(f'  Outputs: {out_json}')
    log.info(f'           {out_tsv}')
    log.info('=' * 70)


if __name__ == '__main__':
    main()
