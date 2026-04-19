"""
Step 8 — Full-cohort IS element burden analysis
─────────────────────────────────────────────────────────────────────────────
Extends IS element analysis from 47 resistant genomes to ALL 270 validated
genomes. Core questions:

  Q1. Is IS6-family element density significantly higher in carbapenem-
      resistant genomes than in susceptible ones?
      → Mann-Whitney U test + Cliff's delta effect size

  Q2. Does IS6 burden correlate with assembly year?
      → Spearman correlation (temporal accumulation hypothesis)

  Q3. Can IS6 burden predict carbapenem resistance?
      → Logistic regression + ROC/AUC analysis

  Q4. Do high-IS6 susceptible genomes represent "resistance-primed" strains?
      → Density comparison; identify high-IS6 susceptible outliers

Outputs:
  data/processed/is_burden_all.tsv       per-genome IS counts + metadata
  data/processed/is_burden_stats.json    all statistical results
  figures/fig_is_burden_violin.pdf/png
  figures/fig_is_burden_temporal.pdf/png
  figures/fig_is_burden_roc.pdf/png
  figures/fig_is_burden_combined.pdf/png
  logs/is_burden_all.log

Usage:
  python analysis/08_is_burden_all.py
"""
import gzip
import json
import logging
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, FIGURES, LOGS, MANIFEST, IS_GFF_KEYWORDS

LOGS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'is_burden_all.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('is_burden_all')

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.8,
})

# IS6 family identifiers (IS26, IS257, IS1353, etc.)
_IS6_NAMES = re.compile(r'\bIS(26|257|1353|240|1006|6)\b', re.IGNORECASE)
_IS_PATTERNS = [re.compile(p, re.IGNORECASE) for p in IS_GFF_KEYWORDS]
_ATTRS_RE = re.compile(r'(\w+)=([^;]+)')


def _parse_attrs(s: str) -> dict:
    return dict(m.groups() for m in _ATTRS_RE.finditer(s))


# ── Per-genome IS element counting ───────────────────────────────────────────

def count_is_elements(gff_gz: Path) -> dict:
    """
    Parse GFF3 and count IS elements by family.
    Returns dict: n_is_total, n_is6, n_is26, n_transposase, n_mobile_element
    """
    n_total = n_is6 = n_is26 = 0
    try:
        with gzip.open(gff_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                feat_type = cols[2]
                if feat_type not in ('CDS', 'gene', 'mobile_element', 'repeat_region'):
                    continue
                attrs   = _parse_attrs(cols[8])
                product = attrs.get('product', '')
                gene    = attrs.get('gene', '')
                note    = attrs.get('Note', '')
                combined = f'{product} {gene} {note}'

                if not any(p.search(combined) for p in _IS_PATTERNS):
                    continue

                n_total += 1
                if _IS6_NAMES.search(combined):
                    n_is6 += 1
                    if re.search(r'\bIS26\b', combined, re.IGNORECASE):
                        n_is26 += 1

    except Exception:
        pass

    return {
        'n_is_total': n_total,
        'n_is6':      n_is6,
        'n_is26':     n_is26,
    }


def extract_assembly_date(stats_txt: Path) -> str | None:
    """Extract # Date: YYYY-MM-DD from assembly_stats.txt."""
    try:
        for line in stats_txt.read_text(errors='replace').splitlines():
            m = re.match(r'#\s*Date:\s*(\d{4}-\d{2}-\d{2})', line)
            if m:
                return m.group(1)
    except Exception:
        pass
    return None


# ── Statistics ────────────────────────────────────────────────────────────────

def cliffs_delta(a, b) -> float:
    """Cliff's delta effect size for two independent samples."""
    a, b = np.array(a), np.array(b)
    if len(a) == 0 or len(b) == 0:
        return float('nan')
    mat = np.sign(a[:, None] - b[None, :])
    return float(mat.mean())


def logistic_regression_roc(X, y):
    """Simple logistic regression without sklearn; returns AUC via trapezoidal rule."""
    from scipy.special import expit
    from scipy.optimize import minimize

    X = np.array(X, dtype=float)
    y = np.array(y, dtype=float)
    X_std = (X - X.mean()) / (X.std() + 1e-10)
    X_design = np.column_stack([np.ones(len(X_std)), X_std])

    def neg_log_likelihood(params):
        p = expit(X_design @ params)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))

    res = minimize(neg_log_likelihood, [0.0, 0.0], method='BFGS')
    params = res.x

    scores = expit(X_design @ params)
    thresholds = np.sort(np.unique(scores))[::-1]
    tprs, fprs = [], []
    pos, neg = y.sum(), (1 - y).sum()
    for t in thresholds:
        pred = (scores >= t).astype(float)
        tp = (pred * y).sum()
        fp = (pred * (1 - y)).sum()
        tprs.append(tp / pos if pos > 0 else 0)
        fprs.append(fp / neg if neg > 0 else 0)
    tprs, fprs = np.array(tprs), np.array(fprs)
    # Sort by fprs ascending for correct AUC
    order = np.argsort(fprs)
    fprs_s, tprs_s = fprs[order], tprs[order]
    trapz_fn = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
    auc = float(trapz_fn(tprs_s, fprs_s))
    fprs, tprs = fprs_s, tprs_s
    return scores, fprs, tprs, auc, params, X.mean(), X.std()


def _save(fig, stem):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGURES / f'{stem}.{ext}', bbox_inches='tight', dpi=300)
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str)
    amr_df   = pd.read_csv(DATA_PROC / 'amr_hits.tsv', sep='\t', dtype=str)

    resistant_accs = set(
        amr_df[amr_df['is_carbapenem'].str.lower() == 'true']['accession']
    )
    n_total = len(manifest)
    log.info(f'Full-cohort IS burden analysis: {n_total} genomes')
    log.info(f'Resistant: {len(resistant_accs)}  Susceptible: {n_total - len(resistant_accs)}')

    # ── Scan all 270 GFF3 files ───────────────────────────────────────────────
    records = []
    for i, (_, row) in enumerate(manifest.iterrows(), 1):
        accession = row['accession']
        local_dir = Path(row['local_dir'])

        gff_files   = list(local_dir.glob('*_genomic.gff.gz'))
        stats_files = list(local_dir.glob('*_assembly_stats.txt'))

        counts = count_is_elements(gff_files[0]) if gff_files else \
                 {'n_is_total': 0, 'n_is6': 0, 'n_is26': 0}
        asm_date = extract_assembly_date(stats_files[0]) if stats_files else None
        asm_year = int(asm_date[:4]) if asm_date else None

        rec = {
            'accession':    accession,
            'is_resistant': accession in resistant_accs,
            'assembly_date': asm_date,
            'assembly_year': asm_year,
            **counts,
        }
        records.append(rec)

        if i % 50 == 0:
            log.info(f'  Scanned {i}/{n_total}')

    df = pd.DataFrame(records)
    df.to_csv(DATA_PROC / 'is_burden_all.tsv', sep='\t', index=False)
    log.info(f'Per-genome IS counts saved → {DATA_PROC / "is_burden_all.tsv"}')

    # ── Split resistant / susceptible ─────────────────────────────────────────
    res_df = df[df['is_resistant']].copy()
    sus_df = df[~df['is_resistant']].copy()

    log.info(f'\n=== Q1: IS6 BURDEN — RESISTANT vs SUSCEPTIBLE ===')
    for col, label in [('n_is_total', 'Total IS'), ('n_is6', 'IS6 family'), ('n_is26', 'IS26 specific')]:
        r_vals = res_df[col].dropna().values
        s_vals = sus_df[col].dropna().values
        stat, p = mannwhitneyu(r_vals, s_vals, alternative='two-sided')
        delta   = cliffs_delta(r_vals, s_vals)
        r_med   = np.median(r_vals)
        s_med   = np.median(s_vals)
        fold    = r_med / s_med if s_med > 0 else float('inf')
        log.info(f'  {label}:')
        log.info(f'    Resistant:   median={r_med:.0f}  IQR=[{np.percentile(r_vals,25):.0f}–{np.percentile(r_vals,75):.0f}]')
        log.info(f'    Susceptible: median={s_med:.0f}  IQR=[{np.percentile(s_vals,25):.0f}–{np.percentile(s_vals,75):.0f}]')
        log.info(f'    Fold change: {fold:.2f}×   Mann-Whitney U p={p:.4e}   Cliff\'s δ={delta:.3f}')

    # ── Q2: Temporal correlation ──────────────────────────────────────────────
    log.info(f'\n=== Q2: IS6 BURDEN TEMPORAL TREND ===')
    year_df = df[df['assembly_year'].notna()].copy()
    log.info(f'  Genomes with date: {len(year_df)}/{n_total}')
    log.info(f'  Year range: {year_df["assembly_year"].min()} – {year_df["assembly_year"].max()}')

    if len(year_df) > 10:
        for col, label in [('n_is6', 'IS6'), ('n_is26', 'IS26')]:
            rho, p_rho = stats.spearmanr(year_df['assembly_year'], year_df[col])
            log.info(f'  {label} vs year: Spearman ρ={rho:.3f}  p={p_rho:.4e}')
        # Year × Resistance interaction
        for yr_grp, sub in year_df.groupby('assembly_year'):
            n_res = sub['is_resistant'].sum()
            log.info(f'  {yr_grp}: n={len(sub)}  resistant={n_res} ({n_res/len(sub)*100:.0f}%)  '
                     f'median IS6={sub["n_is6"].median():.0f}')

    # ── Q3: Logistic regression / ROC ────────────────────────────────────────
    log.info(f'\n=== Q3: IS6 AS RESISTANCE PREDICTOR ===')
    complete = df[df['n_is6'].notna()].copy()
    X = complete['n_is6'].values
    y = complete['is_resistant'].astype(int).values
    scores, fprs, tprs, auc, params, x_mean, x_std = logistic_regression_roc(X, y)
    log.info(f'  Logistic regression AUC (IS6 → resistance): {auc:.3f}')
    log.info(f'  Model intercept={params[0]:.3f}  IS6_coef={params[1]:.3f}')
    log.info(f'  Interpretation: each 1-SD increase in IS6 count → '
             f'OR={np.exp(params[1]):.2f}×')

    # ── Q4: Resistance-primed susceptible strains ─────────────────────────────
    log.info(f'\n=== Q4: HIGH-IS6 SUSCEPTIBLE STRAINS (resistance-primed?) ===')
    res_is6_75pct = np.percentile(res_df['n_is6'], 25)   # lowest quartile of resistant
    primed = sus_df[sus_df['n_is6'] >= res_is6_75pct]
    log.info(f'  IS6 threshold (resistant 25th pct): {res_is6_75pct:.0f}')
    log.info(f'  Susceptible strains above threshold: {len(primed)}/{len(sus_df)} '
             f'({len(primed)/len(sus_df)*100:.1f}%)')
    log.info(f'  → These are "resistance-primed" candidates for surveillance')
    if not primed.empty:
        log.info(f'  Primed strain accessions:')
        for acc in primed['accession'].tolist()[:10]:
            log.info(f'    {acc}')

    # ── Collect all stats for JSON ────────────────────────────────────────────
    r_is6 = res_df['n_is6'].values
    s_is6 = sus_df['n_is6'].values
    r_tot = res_df['n_is_total'].values
    s_tot = sus_df['n_is_total'].values
    stat_is6, p_is6 = mannwhitneyu(r_is6, s_is6, alternative='two-sided')
    stat_tot, p_tot = mannwhitneyu(r_tot, s_tot, alternative='two-sided')
    delta_is6 = cliffs_delta(r_is6, s_is6)
    delta_tot  = cliffs_delta(r_tot, s_tot)
    rho_is6 = rho_p_is6 = None
    if len(year_df) > 10:
        rho_is6, rho_p_is6 = stats.spearmanr(year_df['assembly_year'], year_df['n_is6'])

    stats_out = {
        'cohort': {'n_total': n_total, 'n_resistant': int(len(res_df)),
                   'n_susceptible': int(len(sus_df))},
        'IS6_resistant_vs_susceptible': {
            'resistant_median': float(np.median(r_is6)),
            'resistant_iqr': [float(np.percentile(r_is6,25)), float(np.percentile(r_is6,75))],
            'susceptible_median': float(np.median(s_is6)),
            'susceptible_iqr': [float(np.percentile(s_is6,25)), float(np.percentile(s_is6,75))],
            'fold_change': float(np.median(r_is6)/np.median(s_is6)) if np.median(s_is6)>0 else None,
            'mann_whitney_p': float(p_is6),
            'cliffs_delta': float(delta_is6),
        },
        'total_IS_resistant_vs_susceptible': {
            'resistant_median': float(np.median(r_tot)),
            'susceptible_median': float(np.median(s_tot)),
            'mann_whitney_p': float(p_tot),
            'cliffs_delta': float(delta_tot),
        },
        'temporal': {
            'n_with_date': int(len(year_df)),
            'year_range': [int(year_df['assembly_year'].min()), int(year_df['assembly_year'].max())],
            'spearman_rho_IS6_vs_year': float(rho_is6) if rho_is6 is not None else None,
            'spearman_p': float(rho_p_is6) if rho_p_is6 is not None else None,
        },
        'predictive': {
            'auc_IS6_vs_resistance': float(auc),
            'logistic_IS6_coef': float(params[1]),
            'OR_per_SD_increase': float(np.exp(params[1])),
        },
        'primed_susceptible': {
            'threshold_IS6': float(res_is6_75pct),
            'n_primed': int(len(primed)),
            'pct_of_susceptible': float(len(primed)/len(sus_df)*100),
            'accessions': primed['accession'].tolist(),
        }
    }

    with open(DATA_PROC / 'is_burden_stats.json', 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2)
    log.info(f'\nStats JSON → {DATA_PROC / "is_burden_stats.json"}')

    # ── Figures ───────────────────────────────────────────────────────────────
    _make_figures(df, res_df, sus_df, primed, year_df, fprs, tprs, auc, stats_out)

    log.info('═' * 60)
    log.info('NOVEL FINDINGS SUMMARY')
    log.info(f'  IS6 fold-change (resistant/susceptible): '
             f'{np.median(r_is6)/np.median(s_is6):.2f}×  p={p_is6:.2e}')
    log.info(f'  Cliff\'s δ = {delta_is6:.3f}  '
             f'(|δ|>0.47=large, >0.33=medium, >0.11=small)')
    log.info(f'  AUC (IS6 as resistance predictor): {auc:.3f}')
    log.info(f'  Resistance-primed susceptible strains: {len(primed)} '
             f'({len(primed)/len(sus_df)*100:.1f}%)')
    if rho_is6 is not None:
        log.info(f'  IS6 temporal trend: ρ={rho_is6:.3f}  p={rho_p_is6:.2e}')


# ── Figure generation ────────────────────────────────────────────────────────

def _make_figures(df, res_df, sus_df, primed, year_df, fprs, tprs, auc, stats_out):
    C_RES = '#E63946'
    C_SUS = '#457B9D'
    C_PRM = '#F4A261'  # primed (high IS6 susceptible)

    # ── Fig A: violin / boxplot IS6 burden ───────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(10, 4.5))

    ax = axes[0]
    data_v = [sus_df['n_is6'].values, res_df['n_is6'].values]
    parts  = ax.violinplot(data_v, positions=[0, 1], showmedians=True)
    for i, (body, col) in enumerate(zip(parts['bodies'], [C_SUS, C_RES])):
        body.set_facecolor(col); body.set_alpha(0.7)
    parts['cmedians'].set_colors(['white', 'white']); parts['cmedians'].set_linewidth(2)
    parts['cmins'].set_colors([C_SUS, C_RES]); parts['cmaxes'].set_colors([C_SUS, C_RES])
    parts['cbars'].set_colors([C_SUS, C_RES])
    p_val = stats_out['IS6_resistant_vs_susceptible']['mann_whitney_p']
    delta = stats_out['IS6_resistant_vs_susceptible']['cliffs_delta']
    fold  = stats_out['IS6_resistant_vs_susceptible']['fold_change']
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Susceptible\n(n={len(sus_df)})', f'Resistant\n(n={len(res_df)})'])
    ax.set_ylabel('IS6-family element count per genome')
    fold_str = f'{fold:.1f}×' if fold is not None and np.isfinite(fold) else '∞×'
    ax.set_title(f'A  IS6 burden by resistance status\n'
                 f'p={p_val:.2e}  δ={delta:.2f}  {fold_str} fold',
                 loc='left', fontweight='bold', fontsize=9)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

    # ── Fig B: scatter IS6 vs year ───────────────────────────────────────────
    ax2 = axes[1]
    if not year_df.empty:
        jitter = np.random.default_rng(42).uniform(-0.3, 0.3, len(year_df))
        colors_scatter = [C_RES if r else C_SUS for r in year_df['is_resistant']]
        ax2.scatter(year_df['assembly_year'] + jitter, year_df['n_is6'],
                    c=colors_scatter, alpha=0.5, s=18, linewidths=0)
        # Trend line
        yr = year_df['assembly_year'].values
        is6 = year_df['n_is6'].values
        z = np.polyfit(yr, is6, 1)
        xs = np.linspace(yr.min(), yr.max(), 100)
        ax2.plot(xs, np.polyval(z, xs), 'k--', lw=1.5, alpha=0.7)
        rho  = stats_out['temporal']['spearman_rho_IS6_vs_year']
        p_r  = stats_out['temporal']['spearman_p']
        ax2.set_xlabel('Assembly year')
        ax2.set_ylabel('IS6-family element count')
        ax2.set_title(f'B  IS6 burden over time\nSpearman ρ={rho:.3f}  p={p_r:.2e}',
                      loc='left', fontweight='bold', fontsize=9)
        patches = [mpatches.Patch(color=C_RES, label='Resistant'),
                   mpatches.Patch(color=C_SUS, label='Susceptible')]
        ax2.legend(handles=patches, frameon=False, fontsize=7.5)
        for sp in ('top','right'): ax2.spines[sp].set_visible(False)

    # ── Fig C: ROC curve ─────────────────────────────────────────────────────
    ax3 = axes[2]
    ax3.plot(fprs, tprs, color=C_RES, lw=2)
    ax3.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
    ax3.fill_between(fprs, tprs, alpha=0.15, color=C_RES)
    ax3.set_xlabel('False positive rate')
    ax3.set_ylabel('True positive rate')
    ax3.set_title(f'C  ROC: IS6 count → resistance\nAUC = {auc:.3f}',
                  loc='left', fontweight='bold', fontsize=9)
    ax3.set_xlim(0, 1); ax3.set_ylim(0, 1.05)
    for sp in ('top','right'): ax3.spines[sp].set_visible(False)

    fig.tight_layout()
    _save(fig, 'fig_is_burden_violin')

    # ── Fig D: combined 4-panel ───────────────────────────────────────────────
    fig2, axes2 = plt.subplots(2, 2, figsize=(11, 9))

    # D1: violin (same as A)
    ax = axes2[0, 0]
    parts2 = ax.violinplot(data_v, positions=[0, 1], showmedians=True)
    for body, col in zip(parts2['bodies'], [C_SUS, C_RES]):
        body.set_facecolor(col); body.set_alpha(0.7)
    parts2['cmedians'].set_colors(['white','white']); parts2['cmedians'].set_linewidth(2)
    parts2['cmins'].set_colors([C_SUS, C_RES]); parts2['cmaxes'].set_colors([C_SUS, C_RES])
    parts2['cbars'].set_colors([C_SUS, C_RES])
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Susceptible (n={len(sus_df)})', f'Resistant (n={len(res_df)})'])
    ax.set_ylabel('IS6-family count')
    ax.set_title(f'A  IS6 burden: p={p_val:.2e}  δ={delta:.2f}  {fold_str}',
                 loc='left', fontweight='bold', fontsize=9)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

    # D2: temporal scatter
    ax = axes2[0, 1]
    if not year_df.empty:
        ax.scatter(year_df['assembly_year'] + jitter, year_df['n_is6'],
                   c=colors_scatter, alpha=0.5, s=18, linewidths=0)
        ax.plot(xs, np.polyval(z, xs), 'k--', lw=1.5, alpha=0.7)
        ax.set_xlabel('Assembly year'); ax.set_ylabel('IS6-family count')
        ax.set_title(f'B  Temporal trend: ρ={rho:.3f}  p={p_r:.2e}',
                     loc='left', fontweight='bold', fontsize=9)
        ax.legend(handles=patches, frameon=False, fontsize=7.5)
        for sp in ('top','right'): ax.spines[sp].set_visible(False)

    # D3: ROC
    ax = axes2[1, 0]
    ax.plot(fprs, tprs, color=C_RES, lw=2)
    ax.plot([0,1],[0,1],'k--',lw=1,alpha=0.5)
    ax.fill_between(fprs, tprs, alpha=0.15, color=C_RES)
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title(f'C  Predictive AUC = {auc:.3f}', loc='left', fontweight='bold', fontsize=9)
    ax.set_xlim(0,1); ax.set_ylim(0,1.05)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

    # D4: "primed susceptible" IS6 distribution with threshold
    ax = axes2[1, 1]
    n_bins = 30
    threshold = stats_out['primed_susceptible']['threshold_IS6']
    s_is6_arr = sus_df['n_is6'].values
    r_is6_arr = res_df['n_is6'].values
    all_max = max(s_is6_arr.max(), r_is6_arr.max())
    bins = np.linspace(0, all_max, n_bins + 1)
    ax.hist(s_is6_arr, bins=bins, color=C_SUS, alpha=0.6, label='Susceptible', density=True)
    ax.hist(r_is6_arr, bins=bins, color=C_RES, alpha=0.6, label='Resistant', density=True)
    ax.axvline(threshold, color=C_PRM, lw=2, ls='--',
               label=f'Primed threshold\n(IS6≥{threshold:.0f}, n={len(primed)})')
    ax.set_xlabel('IS6-family element count')
    ax.set_ylabel('Density')
    ax.set_title(f'D  Resistance-primed susceptible strains\n'
                 f'{len(primed)}/{len(sus_df)} ({len(primed)/len(sus_df)*100:.0f}%) above threshold',
                 loc='left', fontweight='bold', fontsize=9)
    ax.legend(frameon=False, fontsize=7.5)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)

    fig2.suptitle('IS6-family element burden across 270 Chinese clinical K. pneumoniae',
                  fontweight='bold', fontsize=11)
    fig2.tight_layout()
    _save(fig2, 'fig_is_burden_combined')

    log.info(f'Figures saved → {FIGURES}')


if __name__ == '__main__':
    main()
