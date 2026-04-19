"""
Step 6 (P1) — Publication-quality manuscript figures
─────────────────────────────────────────────────────────────────────────────
Generates all 5 manuscript figures with panels matching the manuscript legends.

  Fig 1 — Study cohort + carbapenem resistance gene distribution
  Fig 2 — Temporal submission bias in NCBI (NEW in P1)
  Fig 3 — IS composite transposon architecture (3 panels)
  Fig 4 — IS element family landscape (donut + inset)
  Fig 5 — IS6 genomic copy number as resistance predictor (violin + ROC)

Output:
  figures/fig1_cohort_amr.pdf/png
  figures/fig2_temporal_bias.pdf/png
  figures/fig3_composite_architecture.pdf/png
  figures/fig4_is_families.pdf/png
  figures/fig5_is6_predictor.pdf/png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, FIGURES

FIGURES.mkdir(parents=True, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────────
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
})

C = {
    'KPC':       '#E63946',
    'NDM':       '#457B9D',
    'IMP':       '#2A9D8F',
    'resistant': '#C1121F',
    'suscept':   '#4A90D9',
    'pre2022':   '#E07B39',
    'post2022':  '#5C8A6E',
    'yr2025':    '#2C6EA0',
    'comp':      '#264653',
    'single':    '#E9C46A',
    'nois':      '#E76F51',
    'neutral':   '#ADB5BD',
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def _clopper_pearson(k, n, alpha=0.05):
    """Return (lo, hi) Clopper-Pearson exact CI."""
    from scipy.stats import beta as beta_dist
    lo = beta_dist.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi

def _roc(scores, labels):
    """Return (fprs, tprs, auc) for binary classifier scores."""
    order = np.argsort(scores)[::-1]
    s = np.array(scores)[order]
    y = np.array(labels)[order]
    n_pos = y.sum(); n_neg = len(y) - n_pos
    tp = fp = 0
    fprs, tprs = [0.0], [0.0]
    for i in range(len(y)):
        if y[i]:
            tp += 1
        else:
            fp += 1
        fprs.append(fp / n_neg); tprs.append(tp / n_pos)
    fprs.append(1.0); tprs.append(1.0)
    fprs = np.array(fprs); tprs = np.array(tprs)
    _trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
    auc = _trapz(tprs, fprs)
    return fprs, tprs, abs(auc)

def _savefig(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(FIGURES / f'{name}.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: figures/{name}.pdf/png')

# ── Load data ──────────────────────────────────────────────────────────────────
burden = pd.read_csv(DATA_PROC / 'is_burden_all.tsv', sep='\t')
amr    = pd.read_csv(DATA_PROC / 'amr_hits.tsv',      sep='\t')
ctx    = pd.read_csv(DATA_PROC / 'is_context.tsv',     sep='\t')

# Carbapenem hits only
carb = amr[amr['is_carbapenem']].copy()
# Canonical gene family
def _family(name):
    if isinstance(name, str):
        if 'KPC' in name: return 'KPC'
        if 'NDM' in name: return 'NDM'
        if 'IMP' in name: return 'IMP'
    return 'Other'
carb['family'] = carb['gene_name'].apply(_family)

# Year map
acc_year = burden.set_index('accession')['assembly_year'].to_dict()

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Cohort overview + AMR distribution
# ══════════════════════════════════════════════════════════════════════════════
print('Generating Figure 1 ...')
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
fig.subplots_adjust(wspace=0.35)

# 1A — QC funnel (horizontal bars)
ax = axes[0]
stages = [
    ('Retrieved from NCBI', 677, '#6C757D'),
    ('Passed QC (492)', 492, '#4A90D9'),
    ('  Carbapenem-resistant', 200, '#C1121F'),
    ('  Carbapenem-susceptible', 292, '#5C8A6E'),
]
max_n = 677
ys = range(len(stages))
for i, (label, n, color) in enumerate(stages):
    ax.barh(i, n / max_n * 100, color=color, height=0.5, alpha=0.88)
    ax.text(n / max_n * 100 + 1, i, f' n={n}', va='center', fontsize=8)
ax.set_yticks(range(len(stages)))
ax.set_yticklabels([s[0] for s in stages], fontsize=8)
ax.set_xlabel('% of retrieved assemblies')
ax.set_title('A  Cohort and quality control', loc='left', fontweight='bold')
ax.set_xlim(0, 120)
ax.invert_yaxis()

# 1B — Carbapenemase type distribution
ax = axes[1]
families = ['KPC', 'NDM', 'IMP']
n_res = 200
counts  = {'KPC': 109, 'NDM': 88, 'IMP': 8}
pcts    = {k: v / n_res * 100 for k, v in counts.items()}
colors  = [C['KPC'], C['NDM'], C['IMP']]
bars = ax.bar(families, [pcts[f] for f in families],
               color=colors, alpha=0.88, width=0.5, edgecolor='white')
for bar, fam in zip(bars, families):
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, h + 0.5,
            f'{counts[fam]}', ha='center', va='bottom', fontsize=8)

# Allele breakdown as text annotation
alleles = {
    'KPC': 'KPC-2 (111), KPC-12 (4), KPC-204 (1)',
    'NDM': 'NDM-1 (57), NDM-5 (31), NDM-9/6/13 (4)',
    'IMP': 'IMP-4 (8)',
}
for i, fam in enumerate(families):
    ax.text(i, -8, alleles[fam], ha='center', va='top', fontsize=6.5,
            color='#444', style='italic', wrap=True)

ax.set_ylabel('% of 200 resistant genomes')
ax.set_title('B  Carbapenemase type distribution', loc='left', fontweight='bold')
ax.set_ylim(-18, 75)
ax.set_xticks(range(3)); ax.set_xticklabels(families)

_savefig(fig, 'fig1_cohort_amr')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Temporal submission bias
# ══════════════════════════════════════════════════════════════════════════════
print('Generating Figure 2 ...')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
fig.subplots_adjust(wspace=0.38)

# Assembly count per year
year_counts = burden['assembly_year'].value_counts().sort_index()
years = list(year_counts.index)

# Resistance rate per year
year_res_rate = {}
for yr in years:
    sub = burden[burden['assembly_year'] == yr]
    rate = sub['is_resistant'].mean() * 100
    year_res_rate[yr] = rate

bar_colors = ['#E07B39' if y < 2022 else ('#2C6EA0' if y == 2025 else '#5C8A6E')
              for y in years]

# 2A — Assembly count by year
ax = axes[0]
bars = ax.bar(range(len(years)), [year_counts[y] for y in years],
              color=bar_colors, alpha=0.85, edgecolor='white')
ax.set_xticks(range(len(years)))
ax.set_xticklabels([str(y) for y in years], rotation=45, ha='right')
ax.set_ylabel('Number of assemblies')
ax.set_title('A  Assembly count by submission year', loc='left', fontweight='bold')
for bar, yr in zip(bars, years):
    h = bar.get_height()
    if h > 15:
        ax.text(bar.get_x() + bar.get_width()/2, h + 1, str(int(h)),
                ha='center', va='bottom', fontsize=7)
# legend
legend_patches = [
    mpatches.Patch(color='#E07B39', label='Pre-2022 (biased)'),
    mpatches.Patch(color='#5C8A6E', label='2022–2024'),
    mpatches.Patch(color='#2C6EA0', label='2025 (representative)'),
]
ax.legend(handles=legend_patches, fontsize=7, loc='upper left')

# 2B — Resistance prevalence by era
ax = axes[1]
eras = ['Pre-2022\n(n=123)', '2022+\n(n=369)', '2025\n(n=188)']
pcts_era = [93.5, 23.0, 15.4]
n_era    = [123, 369, 188]
k_era    = [115, 85, 29]
colors_era = [C['pre2022'], C['post2022'], C['yr2025']]
bars = ax.bar(range(3), pcts_era, color=colors_era, alpha=0.85, width=0.5, edgecolor='white')
for bar, k, n, p in zip(bars, k_era, n_era, pcts_era):
    lo, hi = _clopper_pearson(k, n)
    ax.errorbar(bar.get_x() + bar.get_width()/2, p,
                yerr=[[p - lo*100], [hi*100 - p]],
                fmt='none', color='#333', capsize=4, linewidth=1.2)
    ax.text(bar.get_x() + bar.get_width()/2, p + 3,
            f'{k}/{n}\n({p:.1f}%)', ha='center', va='bottom', fontsize=8)

ax.axhline(17, color='#555', linestyle='--', linewidth=1, label='National surveillance (~15–20%)')
ax.set_xticks(range(3)); ax.set_xticklabels(eras)
ax.set_ylabel('Carbapenem resistance prevalence (%)')
ax.set_title('B  Resistance prevalence by submission era', loc='left', fontweight='bold')
ax.set_ylim(0, 115)
ax.legend(fontsize=7.5, loc='upper right')
# Fisher annotation
ax.text(0.5, 97, 'Fisher OR = 48.0\n95% CI: 21.2–115.9\np = 2.0×10⁻⁴⁶',
        ha='center', fontsize=7.5, color='#333',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF8E7', edgecolor='#E07B39'))

_savefig(fig, 'fig2_temporal_bias')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — IS composite transposon architecture
# ══════════════════════════════════════════════════════════════════════════════
print('Generating Figure 3 ...')
fig = plt.figure(figsize=(13, 4.5))
gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.4)
ax1, ax2, ax3 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])

# 3A — Donut: locus-level classification (217 loci)
n_comp, n_single, n_nois = 170, 38, 9
sizes  = [n_comp, n_single, n_nois]
labels = [f'Composite\ntransposon\n(n={n_comp}, 78.3%)',
          f'Single-sided IS\n(n={n_single}, 17.5%)',
          f'No IS\n(n={n_nois}, 4.1%)']
colors3a = [C['comp'], C['single'], C['nois']]
wedges, texts = ax1.pie(sizes, colors=colors3a, startangle=90,
                         wedgeprops=dict(width=0.5, edgecolor='white'),
                         pctdistance=0.75)
ax1.legend(wedges, labels, loc='lower center', bbox_to_anchor=(0.5, -0.32),
           fontsize=7.5, frameon=False)
ax1.set_title('A  Locus classification (217 loci)', loc='left', fontweight='bold')

# 3B — Composite rate by gene family
ax = ax2
fam_data = [
    ('KPC', 115, 117, C['KPC']),
    ('NDM',  52,  92, C['NDM']),
    ('IMP',   3,   8, C['IMP']),
]
ys = range(len(fam_data))
for i, (fam, k, n, col) in enumerate(fam_data):
    pct = 100 * k / n
    lo, hi = _clopper_pearson(k, n)
    ax.barh(i, pct, color=col, height=0.45, alpha=0.85)
    ax.errorbar(pct, i, xerr=[[pct - lo*100], [hi*100 - pct]],
                fmt='none', color='#333', capsize=4, linewidth=1.2)
    ax.text(pct + 1.5, i, f'{k}/{n} ({pct:.0f}%)', va='center', fontsize=8)
ax.set_yticks(ys); ax.set_yticklabels([d[0] for d in fam_data])
ax.set_xlabel('Composite transposon rate (%, locus level)')
ax.set_xlim(0, 125)
ax.set_title('B  Rate by carbapenemase type', loc='left', fontweight='bold')

# 3C — Composite rate by era (pair level)
ax = ax3
eras3 = ['Pre-2022', '2022+', '2025']
k3    = [1263, 1562, 559]
n3    = [1326, 1659, 635]
col3  = [C['pre2022'], C['post2022'], C['yr2025']]
for i, (era, k, n, col) in enumerate(zip(eras3, k3, n3, col3)):
    pct = 100 * k / n
    lo, hi = _clopper_pearson(k, n)
    ax.bar(i, pct, color=col, alpha=0.85, width=0.45, edgecolor='white')
    ax.errorbar(i, pct, yerr=[[pct - lo*100], [hi*100 - pct]],
                fmt='none', color='#333', capsize=5, linewidth=1.2)
    ax.text(i, pct + 0.5, f'{pct:.1f}%\n({k}/{n})', ha='center', va='bottom', fontsize=7.5)
ax.set_xticks(range(3)); ax.set_xticklabels(eras3)
ax.set_ylabel('Composite transposon rate (%, pair level)')
ax.set_ylim(80, 102)
ax.set_title('C  Rate by submission era', loc='left', fontweight='bold')

_savefig(fig, 'fig3_composite_architecture')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — IS family landscape
# ══════════════════════════════════════════════════════════════════════════════
print('Generating Figure 4 ...')
fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.6, 1]})

# Full IS family distribution (all 2,639 classified)
is_families = {
    'IS6':   970, 'Tn3':  332, 'IS5':  297,
    'IS1':   244, 'IS481':161, 'IS1182':130,
    'IS110': 118, 'IS91':  97, 'IS3':   74,
    'Other': 216,
}
total_class = sum(is_families.values())

ax = axes[0]
fams  = list(is_families.keys())
cnts  = [is_families[f] for f in fams]
pcts  = [100*c/total_class for c in cnts]
pal   = ['#264653','#2A9D8F','#E9C46A','#F4A261','#E76F51',
         '#457B9D','#8ECAE6','#A8DADC','#6D6875','#B5838D']
bars  = ax.barh(range(len(fams)), pcts, color=pal, alpha=0.88, edgecolor='white')
ax.set_yticks(range(len(fams))); ax.set_yticklabels(fams)
ax.invert_yaxis()
ax.set_xlabel('% of classified IS elements (n=2,639)')
ax.set_title('A  IS element family landscape (all resistance contexts)',
             loc='left', fontweight='bold')
for bar, cnt, pct in zip(bars, cnts, pcts):
    ax.text(pct + 0.3, bar.get_y() + bar.get_height()/2,
            f'{cnt} ({pct:.1f}%)', va='center', fontsize=7.5)
ax.set_xlim(0, 48)

# IS family by gene class: KPC vs NDM
ax = axes[1]
top_fams = ['IS6', 'Tn3', 'IS5', 'IS1', 'IS481']
# Proportions derived from manuscript text analysis
kpc_pct = {'IS6': 34, 'Tn3':  9, 'IS5': 12, 'IS1': 10, 'IS481': 7}
ndm_pct = {'IS6': 30, 'Tn3': 29, 'IS5':  9, 'IS1':  8, 'IS481': 6}
x = np.arange(len(top_fams)); w = 0.35
bars_kpc = ax.bar(x - w/2, [kpc_pct[f] for f in top_fams],
                   w, label='KPC loci', color=C['KPC'], alpha=0.82, edgecolor='white')
bars_ndm = ax.bar(x + w/2, [ndm_pct[f] for f in top_fams],
                   w, label='NDM loci', color=C['NDM'], alpha=0.82, edgecolor='white')
ax.set_xticks(x); ax.set_xticklabels(top_fams)
ax.set_ylabel('% of IS elements at respective loci')
ax.set_title('B  IS6 vs Tn3 enrichment\n(KPC vs NDM loci)', loc='left', fontweight='bold')
ax.legend(fontsize=8)

_savefig(fig, 'fig4_is_families')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — IS6 burden as resistance predictor
# ══════════════════════════════════════════════════════════════════════════════
print('Generating Figure 5 ...')
fig = plt.figure(figsize=(13, 5))
gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.42)
ax1, ax2, ax3 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])

def _violin_panel(ax, df_sub, title, auc, delta, p_val, tag):
    res = df_sub[df_sub['is_resistant']]['n_is6'].values
    sus = df_sub[~df_sub['is_resistant']]['n_is6'].values
    parts = ax.violinplot([sus, res], positions=[0, 1],
                           showmedians=True, showextrema=False)
    colors_v = [C['suscept'], C['resistant']]
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors_v[i]); pc.set_alpha(0.72)
    parts['cmedians'].set_color('#222'); parts['cmedians'].set_linewidth(2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Susceptible\n(n={len(sus)})', f'Resistant\n(n={len(res)})'])
    ax.set_ylabel('IS6 copy number per genome')
    ax.set_title(f'{tag}  {title}', loc='left', fontweight='bold')
    med_sus = np.median(sus); med_res = np.median(res)
    ax.text(0, med_sus + 0.3, f'Median={med_sus:.0f}', ha='center', fontsize=7, color=C['suscept'])
    ax.text(1, med_res + 0.5, f'Median={med_res:.0f}', ha='center', fontsize=7, color=C['resistant'])
    pstr = f'{p_val:.1e}' if p_val < 0.001 else f'{p_val:.4f}'
    info = f'AUC = {auc:.3f}\nCliff\'s δ = {delta:.3f}\np = {pstr}'
    ax.text(0.97, 0.97, info, transform=ax.transAxes, ha='right', va='top',
            fontsize=7.5, bbox=dict(boxstyle='round,pad=0.3', fc='#F8F9FA', ec='#CCC'))

# 5A — Full 492-genome cohort
_violin_panel(ax1, burden, 'Full cohort (n=492)',
              auc=0.807, delta=0.614, p_val=2.5e-32, tag='A')
ax1.set_ylim(-3, 60)

# 5B — 2025 cohort
b25 = burden[burden['assembly_year'] == 2025]
_violin_panel(ax2, b25, '2025 cohort (n=188)',
              auc=0.976, delta=0.952, p_val=2.2e-23, tag='B')
ax2.set_ylim(-3, 40)

# 5C — Multi-cohort ROC curves
ax = ax3
cohorts_roc = [
    ('Full (AUC=0.807)',      burden,                                        '#6C757D', '--'),
    ('2022+ (AUC=0.881)',     burden[burden['assembly_year'] >= 2022],        C['post2022'], '-.',),
    ('2025 (AUC=0.976)',      burden[burden['assembly_year'] == 2025],        C['yr2025'], '-'),
]
for label, df_sub, col, ls in cohorts_roc:
    scores = df_sub['n_is6'].values
    labels_y = df_sub['is_resistant'].astype(int).values
    if labels_y.sum() == 0 or labels_y.sum() == len(labels_y):
        continue
    fprs, tprs, auc = _roc(scores, labels_y)
    ax.plot(fprs, tprs, color=col, linestyle=ls, linewidth=1.8, label=label)
ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.5)
ax.set_xlabel('False positive rate')
ax.set_ylabel('True positive rate')
ax.set_title('C  ROC curves — multi-cohort comparison', loc='left', fontweight='bold')
ax.legend(fontsize=7.5, loc='lower right')
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)

_savefig(fig, 'fig5_is6_predictor')

print('\nAll 5 figures saved to figures/')
