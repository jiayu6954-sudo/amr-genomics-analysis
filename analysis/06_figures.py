"""
Step 6 — Publication-quality figures
─────────────────────────────────────────────────────────────────────────────
Generates all manuscript figures:

  Fig 1 — Carbapenem resistance prevalence and gene family distribution
  Fig 2 — IS element flanking context (stacked bar per gene family)
  Fig 3 — IS family composition (horizontal bar chart)
  Fig 4 — IS element burden in carbapenem-resistant genomes (boxplot)
  Fig 5 — Composite transposon rates across gene families (forest plot)

Output:
  figures/fig1_prevalence.pdf + .png
  figures/fig2_is_context.pdf + .png
  figures/fig3_is_families.pdf + .png
  figures/fig4_is_burden.pdf + .png
  figures/fig5_forest.pdf + .png
  figures/fig_combined.pdf + .png   (all panels 2×3 grid)

Usage:
  python analysis/06_figures.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, FIGURES, MANIFEST

FIGURES.mkdir(parents=True, exist_ok=True)

# ── Publication style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'DejaVu Sans'],
    'font.size':         9,
    'axes.titlesize':    10,
    'axes.labelsize':    9,
    'xtick.labelsize':   8,
    'ytick.labelsize':   8,
    'legend.fontsize':   8,
    'figure.dpi':        150,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
})

PALETTE = {
    'KPC':  '#E63946',
    'NDM':  '#457B9D',
    'IMP':  '#2A9D8F',
    'VIM':  '#E9C46A',
    'OXA':  '#F4A261',
    'OTHER':'#ADB5BD',
}

IS_PALETTE = {
    'IS6':      '#E63946',
    'IS1':      '#457B9D',
    'Tn3':      '#2A9D8F',
    'IS5':      '#E9C46A',
    'IS481':    '#F4A261',
    'IS1182':   '#A8DADC',
    'IS110':    '#6A994E',
    'IS1380':   '#C77DFF',
    'IS3':      '#F7B2BD',
    'IS_unknown':'#ADB5BD',
}

FLANK_PALETTE = {
    'COMPOSITE_TRANSPOSON': '#E63946',
    'SINGLE_IS_UPSTREAM':   '#457B9D',
    'SINGLE_IS_DOWNSTREAM': '#2A9D8F',
    'NO_IS':                '#CED4DA',
}
FLANK_LABELS = {
    'COMPOSITE_TRANSPOSON': 'Composite transposon',
    'SINGLE_IS_UPSTREAM':   'Single IS (upstream)',
    'SINGLE_IS_DOWNSTREAM': 'Single IS (downstream)',
    'NO_IS':                'No IS element',
}


def _gene_family(g: str) -> str:
    g = str(g).upper()
    for k in ('KPC', 'NDM', 'IMP', 'VIM', 'OXA'):
        if k in g:
            return k
    return 'OTHER'


def _save(fig: plt.Figure, stem: str) -> None:
    for ext in ('pdf', 'png'):
        fig.savefig(FIGURES / f'{stem}.{ext}', bbox_inches='tight', dpi=300)
    plt.close(fig)


# ── Figure 1: Prevalence + gene family bar ────────────────────────────────────
def fig1_prevalence(manifest, amr_df):
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.5))

    # Panel A — prevalence bar with CI
    n_total    = len(manifest)
    carba_df   = amr_df[amr_df['is_carbapenem'].str.lower() == 'true']
    n_res      = carba_df['accession'].nunique()
    n_sus      = n_total - n_res
    pct_res    = n_res / n_total * 100
    pct_sus    = n_sus / n_total * 100

    from scipy.stats import beta as beta_dist
    k, n = n_res, n_total
    ci_lo = beta_dist.ppf(0.025, k, n - k + 1) * 100 if k > 0 else 0
    ci_hi = beta_dist.ppf(0.975, k + 1, n - k) * 100 if k < n else 100

    ax = axes[0]
    ax.bar(['Susceptible', 'Resistant'],
           [pct_sus, pct_res],
           color=['#ADB5BD', '#E63946'],
           width=0.5, zorder=2)
    ax.errorbar(['Resistant'], [pct_res],
                yerr=[[pct_res - ci_lo], [ci_hi - pct_res]],
                fmt='none', color='black', capsize=4, linewidth=1.2)
    ax.set_ylabel('Percentage of genomes (%)')
    ax.set_title(f'A  Carbapenem resistance prevalence\n'
                 f'(n = {n_total} genomes)', loc='left', fontweight='bold')
    ax.set_ylim(0, 105)
    ax.text(1, pct_res + (ci_hi - pct_res) + 3,
            f'{pct_res:.1f}%\n[{ci_lo:.1f}–{ci_hi:.1f}%]',
            ha='center', va='bottom', fontsize=7.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    # Panel B — gene family breakdown
    carba_df = carba_df.copy()
    carba_df['gene_family'] = carba_df['gene_name'].apply(_gene_family)
    fam_counts = carba_df.groupby('gene_family')['accession'].nunique()
    fam_order  = fam_counts.sort_values(ascending=False).index.tolist()

    ax2 = axes[1]
    colors = [PALETTE.get(f, '#ADB5BD') for f in fam_order]
    bars = ax2.bar(fam_order, [fam_counts[f] for f in fam_order],
                   color=colors, width=0.5, zorder=2)
    for bar, fam in zip(bars, fam_order):
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2, h + 0.4,
                 f'{int(h)}\n({h/n_res*100:.0f}%)',
                 ha='center', va='bottom', fontsize=7.5)
    ax2.set_ylabel('Number of genomes')
    ax2.set_title('B  Carbapenemase gene families\n(resistant genomes)',
                  loc='left', fontweight='bold')
    ax2.set_ylim(0, fam_counts.max() * 1.25)
    for sp in ('top', 'right'):
        ax2.spines[sp].set_visible(False)

    fig.tight_layout()
    _save(fig, 'fig1_prevalence')
    return fig


# ── Figure 2: IS flanking context stacked bar per gene family ─────────────────
def fig2_is_context(ctx_df):
    fig, ax = plt.subplots(figsize=(6, 4))

    ctx_df = ctx_df.copy()
    ctx_df['gene_family'] = ctx_df['amr_gene'].apply(_gene_family)
    classes = list(FLANK_PALETTE.keys())
    fam_order = ctx_df.groupby('gene_family')['accession'].nunique().sort_values(ascending=False).index.tolist()

    pivot = pd.crosstab(ctx_df['gene_family'], ctx_df['flank_class'])
    for cls in classes:
        if cls not in pivot.columns:
            pivot[cls] = 0
    pivot = pivot[classes].loc[fam_order]
    pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100

    bottom = np.zeros(len(fam_order))
    for cls in classes:
        vals = pivot_pct[cls].values
        ax.bar(fam_order, vals, bottom=bottom,
               color=FLANK_PALETTE[cls], label=FLANK_LABELS[cls],
               width=0.55, zorder=2)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 4:
                ax.text(i, b + v / 2, f'{v:.0f}%',
                        ha='center', va='center', fontsize=7.5,
                        color='white', fontweight='bold')
        bottom += vals

    ax.set_ylabel('Percentage of IS–AMR context pairs (%)')
    ax.set_title('IS element flanking context\nby carbapenemase gene family',
                 loc='left', fontweight='bold')
    ax.set_ylim(0, 110)
    ax.legend(loc='upper right', frameon=False)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    _save(fig, 'fig2_is_context')
    return fig


# ── Figure 3: IS family composition (horizontal bar) ─────────────────────────
def fig3_is_families(ctx_df):
    fig, ax = plt.subplots(figsize=(6, 4.5))

    fam_counts = (
        ctx_df[ctx_df['is_family'].notna() & (ctx_df['is_family'] != '')]
        ['is_family'].value_counts()
    )
    total = fam_counts.sum()
    fam_top = fam_counts.head(10)
    labels  = fam_top.index.tolist()
    values  = fam_top.values
    colors  = [IS_PALETTE.get(l, '#ADB5BD') for l in labels]

    ax.barh(labels[::-1], values[::-1], color=colors[::-1], zorder=2)
    for v, (i, l) in zip(values[::-1], enumerate(labels[::-1])):
        ax.text(v + 2, i, f'{v}  ({v/total*100:.1f}%)', va='center', fontsize=8)

    ax.set_xlabel('Number of IS–AMR context pairs')
    ax.set_title('IS element family distribution\n(carbapenem resistance gene context)',
                 loc='left', fontweight='bold')
    ax.set_xlim(0, values.max() * 1.4)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    _save(fig, 'fig3_is_families')
    return fig


# ── Figure 4: IS burden boxplot (resistant genomes by gene family) ────────────
def fig4_is_burden(ctx_df, summ_df):
    fig, ax = plt.subplots(figsize=(5.5, 4))

    summ_df = summ_df.copy()
    summ_df['n_is_features'] = pd.to_numeric(summ_df['n_is_features'], errors='coerce')

    # Get dominant gene family per genome
    gene_fam_map = (
        ctx_df.groupby('accession')['amr_gene']
        .apply(lambda x: _gene_family(x.mode()[0]))
        .reset_index()
    )
    gene_fam_map.columns = ['accession', 'dominant_family']
    merged = summ_df.merge(gene_fam_map, on='accession', how='left')
    merged['dominant_family'] = merged['dominant_family'].fillna('OTHER')

    fam_order = merged.groupby('dominant_family')['n_is_features'].median().sort_values(ascending=False).index.tolist()
    data   = [merged[merged['dominant_family'] == f]['n_is_features'].dropna().values for f in fam_order]
    colors = [PALETTE.get(f, '#ADB5BD') for f in fam_order]

    bp = ax.boxplot(data, positions=range(len(fam_order)), patch_artist=True,
                    widths=0.4, showfliers=True,
                    medianprops={'color': 'black', 'linewidth': 1.5})
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)

    ax.set_xticks(range(len(fam_order)))
    ax.set_xticklabels(fam_order)
    ax.set_ylabel('Number of IS element features in genome')
    ax.set_title('IS element burden\nin carbapenem-resistant genomes',
                 loc='left', fontweight='bold')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    _save(fig, 'fig4_is_burden')
    return fig


# ── Figure 5: Forest plot — composite rate by gene family ─────────────────────
def fig5_forest(ctx_df):
    from scipy.stats import beta as beta_dist

    ctx_df = ctx_df.copy()
    ctx_df['gene_family'] = ctx_df['amr_gene'].apply(_gene_family)

    rows = []
    for fam, sub in ctx_df.groupby('gene_family'):
        n_comp = (sub['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
        n_tot  = len(sub)
        pct    = n_comp / n_tot * 100 if n_tot > 0 else 0
        k = n_comp
        ci_lo = beta_dist.ppf(0.025, k, n_tot - k + 1) * 100 if k > 0 else 0
        ci_hi = beta_dist.ppf(0.975, k + 1, n_tot - k) * 100 if k < n_tot else 100
        rows.append({'family': fam, 'pct': pct, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
                     'n_comp': n_comp, 'n_tot': n_tot})

    # Add overall
    k_all   = (ctx_df['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
    n_all   = len(ctx_df)
    pct_all = k_all / n_all * 100
    ci_lo_a = beta_dist.ppf(0.025, k_all, n_all - k_all + 1) * 100
    ci_hi_a = beta_dist.ppf(0.975, k_all + 1, n_all - k_all) * 100
    rows.append({'family': 'OVERALL', 'pct': pct_all, 'ci_lo': ci_lo_a, 'ci_hi': ci_hi_a,
                 'n_comp': k_all, 'n_tot': n_all})

    rows = sorted(rows, key=lambda r: r['family'] == 'OVERALL')

    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    y_pos = list(range(len(rows)))

    for i, r in enumerate(rows):
        color = '#6C757D' if r['family'] == 'OVERALL' else PALETTE.get(r['family'], '#ADB5BD')
        lw    = 2.5 if r['family'] == 'OVERALL' else 1.8
        ax.errorbar(r['pct'], i,
                    xerr=[[r['pct'] - r['ci_lo']], [r['ci_hi'] - r['pct']]],
                    fmt='o', color=color, capsize=4, linewidth=lw,
                    markersize=8 if r['family'] == 'OVERALL' else 6,
                    markerfacecolor=color)
        ax.text(r['ci_hi'] + 1, i,
                f"{r['pct']:.0f}%  [{r['ci_lo']:.0f}–{r['ci_hi']:.0f}%]  "
                f"({r['n_comp']}/{r['n_tot']})",
                va='center', fontsize=7.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([r['family'] for r in rows])
    ax.set_xlabel('Composite transposon rate (%)')
    ax.set_title('Composite transposon rate by carbapenemase family\n(Clopper–Pearson 95% CI)',
                 loc='left', fontweight='bold')
    ax.set_xlim(0, 130)
    ax.axvline(100, ls='--', color='#CED4DA', lw=0.8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    _save(fig, 'fig5_forest')
    return fig


# ── Combined figure (all panels) ─────────────────────────────────────────────
def fig_combined(manifest, amr_df, ctx_df, summ_df):
    fig = plt.figure(figsize=(12, 10))
    gs  = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.35)

    # ---- panel A: prevalence bar ----
    ax_a = fig.add_subplot(gs[0, 0])
    n_total  = len(manifest)
    carba    = amr_df[amr_df['is_carbapenem'].str.lower() == 'true']
    n_res    = carba['accession'].nunique()
    n_sus    = n_total - n_res
    from scipy.stats import beta as beta_dist
    k, n = n_res, n_total
    ci_lo = beta_dist.ppf(0.025, k, n - k + 1) * 100 if k > 0 else 0
    ci_hi = beta_dist.ppf(0.975, k + 1, n - k) * 100 if k < n else 100
    ax_a.bar(['Susceptible', 'Resistant'], [n_sus/n_total*100, n_res/n_total*100],
             color=['#ADB5BD', '#E63946'], width=0.5, zorder=2)
    ax_a.errorbar(['Resistant'], [n_res/n_total*100],
                  yerr=[[n_res/n_total*100 - ci_lo], [ci_hi - n_res/n_total*100]],
                  fmt='none', color='black', capsize=4)
    ax_a.set_ylabel('%')
    ax_a.set_title('A  Resistance prevalence', loc='left', fontweight='bold', fontsize=9)
    ax_a.set_ylim(0, 105)
    for sp in ('top', 'right'): ax_a.spines[sp].set_visible(False)

    # ---- panel B: gene family ----
    ax_b = fig.add_subplot(gs[0, 1])
    carba2 = carba.copy()
    carba2['gene_family'] = carba2['gene_name'].apply(_gene_family)
    fam = carba2.groupby('gene_family')['accession'].nunique().sort_values(ascending=False)
    ax_b.bar(fam.index, fam.values, color=[PALETTE.get(f,'#ADB5BD') for f in fam.index], width=0.5, zorder=2)
    for i, (f, v) in enumerate(fam.items()):
        ax_b.text(i, v + 0.3, f'{v}', ha='center', fontsize=8)
    ax_b.set_ylabel('Genomes')
    ax_b.set_title('B  Gene family distribution', loc='left', fontweight='bold', fontsize=9)
    for sp in ('top', 'right'): ax_b.spines[sp].set_visible(False)

    # ---- panel C: IS context stacked bar ----
    ax_c = fig.add_subplot(gs[1, 0])
    classes    = list(FLANK_PALETTE.keys())
    ctx2       = ctx_df.copy()
    ctx2['gene_family'] = ctx2['amr_gene'].apply(_gene_family)
    fam_order  = ctx2.groupby('gene_family')['accession'].nunique().sort_values(ascending=False).index.tolist()
    pivot      = pd.crosstab(ctx2['gene_family'], ctx2['flank_class'])
    for cls in classes:
        if cls not in pivot.columns: pivot[cls] = 0
    pivot      = pivot[classes].loc[fam_order]
    pivot_pct  = pivot.div(pivot.sum(axis=1), axis=0) * 100
    bottom     = np.zeros(len(fam_order))
    for cls in classes:
        vals = pivot_pct[cls].values
        ax_c.bar(fam_order, vals, bottom=bottom, color=FLANK_PALETTE[cls],
                 label=FLANK_LABELS[cls], width=0.55, zorder=2)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 5:
                ax_c.text(i, b + v/2, f'{v:.0f}%', ha='center', va='center',
                          fontsize=7.5, color='white', fontweight='bold')
        bottom += vals
    ax_c.set_ylabel('%')
    ax_c.set_title('C  IS flanking context', loc='left', fontweight='bold', fontsize=9)
    ax_c.legend(loc='lower right', frameon=False, fontsize=6.5)
    for sp in ('top', 'right'): ax_c.spines[sp].set_visible(False)

    # ---- panel D: IS families ----
    ax_d = fig.add_subplot(gs[1, 1])
    fam_counts = ctx_df[ctx_df['is_family'].notna() & (ctx_df['is_family'] != '')]['is_family'].value_counts().head(8)
    total_is   = fam_counts.sum()
    labels     = fam_counts.index.tolist()
    values     = fam_counts.values
    colors     = [IS_PALETTE.get(l, '#ADB5BD') for l in labels]
    ax_d.barh(labels[::-1], values[::-1], color=colors[::-1], zorder=2)
    for v, i in zip(values[::-1], range(len(values))):
        ax_d.text(v + 2, i, f'{v} ({v/total_is*100:.0f}%)', va='center', fontsize=7)
    ax_d.set_xlabel('IS–AMR pairs')
    ax_d.set_title('D  IS family composition', loc='left', fontweight='bold', fontsize=9)
    ax_d.set_xlim(0, values.max() * 1.5)
    for sp in ('top', 'right'): ax_d.spines[sp].set_visible(False)

    # ---- panel E: IS burden boxplot ----
    ax_e = fig.add_subplot(gs[2, 0])
    summ2 = summ_df.copy()
    summ2['n_is_features'] = pd.to_numeric(summ2['n_is_features'], errors='coerce')
    gene_fam_map = ctx2.groupby('accession')['gene_family'].agg(lambda x: x.mode()[0]).reset_index()
    merged2 = summ2.merge(gene_fam_map, on='accession', how='left')
    merged2['gene_family'] = merged2['gene_family'].fillna('OTHER')
    fam_order_b = merged2.groupby('gene_family')['n_is_features'].median().sort_values(ascending=False).index.tolist()
    data_b = [merged2[merged2['gene_family'] == f]['n_is_features'].dropna().values for f in fam_order_b]
    colors_b = [PALETTE.get(f, '#ADB5BD') for f in fam_order_b]
    bp = ax_e.boxplot(data_b, positions=range(len(fam_order_b)), patch_artist=True,
                      widths=0.4, showfliers=True,
                      medianprops={'color': 'black', 'linewidth': 1.5})
    for patch, color in zip(bp['boxes'], colors_b):
        patch.set_facecolor(color); patch.set_alpha(0.8)
    ax_e.set_xticks(range(len(fam_order_b)))
    ax_e.set_xticklabels(fam_order_b)
    ax_e.set_ylabel('IS features per genome')
    ax_e.set_title('E  IS element burden', loc='left', fontweight='bold', fontsize=9)
    for sp in ('top', 'right'): ax_e.spines[sp].set_visible(False)

    # ---- panel F: forest plot ----
    ax_f = fig.add_subplot(gs[2, 1])
    rows = []
    for fam2, sub2 in ctx2.groupby('gene_family'):
        n_comp = (sub2['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
        n_tot2  = len(sub2)
        pct2    = n_comp / n_tot2 * 100
        k2 = int(n_comp)
        ci_lo2 = beta_dist.ppf(0.025, k2, n_tot2 - k2 + 1) * 100 if k2 > 0 else 0
        ci_hi2 = beta_dist.ppf(0.975, k2 + 1, n_tot2 - k2) * 100 if k2 < n_tot2 else 100
        rows.append({'family': fam2, 'pct': pct2, 'ci_lo': ci_lo2, 'ci_hi': ci_hi2,
                     'n_comp': k2, 'n_tot': n_tot2})
    k_all  = (ctx2['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
    n_all  = len(ctx2)
    pct_all = k_all / n_all * 100
    ci_lo_a = beta_dist.ppf(0.025, int(k_all), n_all - int(k_all) + 1) * 100
    ci_hi_a = beta_dist.ppf(0.975, int(k_all) + 1, n_all - int(k_all)) * 100
    rows.append({'family': 'OVERALL', 'pct': pct_all, 'ci_lo': ci_lo_a, 'ci_hi': ci_hi_a,
                 'n_comp': k_all, 'n_tot': n_all})
    rows = sorted(rows, key=lambda r: r['family'] == 'OVERALL')
    for i, r in enumerate(rows):
        color = '#6C757D' if r['family'] == 'OVERALL' else PALETTE.get(r['family'], '#ADB5BD')
        ax_f.errorbar(r['pct'], i,
                      xerr=[[r['pct'] - r['ci_lo']], [r['ci_hi'] - r['pct']]],
                      fmt='o', color=color, capsize=4, linewidth=1.8,
                      markersize=7 if r['family'] == 'OVERALL' else 5,
                      markerfacecolor=color)
        ax_f.text(r['ci_hi'] + 1.5, i,
                  f"{r['pct']:.0f}%  ({r['n_comp']}/{r['n_tot']})",
                  va='center', fontsize=7)
    ax_f.set_yticks(range(len(rows)))
    ax_f.set_yticklabels([r['family'] for r in rows])
    ax_f.set_xlabel('Composite transposon rate (%)')
    ax_f.set_title('F  Forest plot: composite rate', loc='left', fontweight='bold', fontsize=9)
    ax_f.set_xlim(0, 130)
    ax_f.axvline(100, ls='--', color='#CED4DA', lw=0.8)
    for sp in ('top', 'right'): ax_f.spines[sp].set_visible(False)

    fig.suptitle('Carbapenem resistance and IS element context\nin Chinese clinical Klebsiella pneumoniae (n = 270)',
                 fontsize=11, fontweight='bold', y=1.01)
    _save(fig, 'fig_combined')
    return fig


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str)
    amr_df   = pd.read_csv(DATA_PROC / 'amr_hits.tsv', sep='\t', dtype=str)
    ctx_df   = pd.read_csv(DATA_PROC / 'is_context.tsv', sep='\t', dtype=str)
    summ_df  = pd.read_csv(DATA_PROC / 'context_summary.tsv', sep='\t')

    print('Generating Figure 1: Prevalence...')
    fig1_prevalence(manifest, amr_df)

    print('Generating Figure 2: IS context stacked bar...')
    fig2_is_context(ctx_df)

    print('Generating Figure 3: IS family distribution...')
    fig3_is_families(ctx_df)

    print('Generating Figure 4: IS burden boxplot...')
    fig4_is_burden(ctx_df, summ_df)

    print('Generating Figure 5: Forest plot...')
    fig5_forest(ctx_df)

    print('Generating combined figure...')
    fig_combined(manifest, amr_df, ctx_df, summ_df)

    print(f'All figures written to {FIGURES}/')
    for f in sorted(FIGURES.glob('*.png')):
        print(f'  {f.name}')


if __name__ == '__main__':
    main()
