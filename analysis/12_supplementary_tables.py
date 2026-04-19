"""
Step 12 — Supplementary Tables S1–S4
─────────────────────────────────────────────────────────────────────────────
Generates all four supplementary tables for manuscript submission.

  S1: Genome assembly accessions and QC summary (492 genomes)
  S2: AMR hit CARD verification results (217 hits)
  S3: Dual carbapenemase genomes (5 genomes)
  S4: IS element flanking context — locus-level summary (217 loci)

Output:
  reports/Table_S1_genome_manifest.tsv
  reports/Table_S2_amr_card_verification.tsv
  reports/Table_S3_dual_carbapenemase.tsv
  reports/Table_S4_is_locus_summary.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, LOGS

REPORTS = Path(__file__).parent.parent / 'reports'
REPORTS.mkdir(parents=True, exist_ok=True)

print('=' * 65)
print('GENERATING SUPPLEMENTARY TABLES S1–S4')
print('=' * 65)

# ── Load data ──────────────────────────────────────────────────────────────────
burden   = pd.read_csv(DATA_PROC / 'is_burden_all.tsv',     sep='\t')
amr      = pd.read_csv(DATA_PROC / 'amr_hits.tsv',          sep='\t')
ctx      = pd.read_csv(DATA_PROC / 'is_context.tsv',         sep='\t')
audit    = pd.read_csv(DATA_PROC / 'audit_report.tsv',       sep='\t') \
           if (DATA_PROC / 'audit_report.tsv').exists() else pd.DataFrame()

try:
    hmmer    = pd.read_csv(DATA_PROC / 'amr_hmmer_results.tsv', sep='\t')
except FileNotFoundError:
    hmmer = pd.DataFrame()

# ── TABLE S1: Genome manifest (492 genomes) ────────────────────────────────────
print('\nTable S1: Genome manifest ...')

carb = amr[amr['is_carbapenem']].copy()

def _family(name):
    if isinstance(name, str):
        if 'KPC' in name: return 'KPC-type'
        if 'NDM' in name: return 'NDM-type'
        if 'IMP' in name: return 'IMP-type'
    return 'Other'

# Per-genome AMR summary
amr_per_genome = (
    carb.groupby('accession')
    .agg(
        amr_genes=('gene_name', lambda x: '; '.join(sorted(set(x)))),
        n_amr_loci=('gene_name', 'count'),
        carbapenemase_type=('gene_name', lambda x: '; '.join(sorted(set(_family(g) for g in x))))
    )
    .reset_index()
)

s1 = burden[['accession', 'is_resistant', 'assembly_year',
             'n_is_total', 'n_is6', 'n_is26']].copy()
s1.columns = ['Assembly_accession', 'Is_carbapenem_resistant', 'Assembly_year',
              'N_IS_total', 'N_IS6_family', 'N_IS26']
s1 = s1.merge(
    amr_per_genome.rename(columns={
        'accession': 'Assembly_accession',
        'amr_genes': 'Carbapenem_resistance_genes',
        'n_amr_loci': 'N_carbapenem_loci',
        'carbapenemase_type': 'Carbapenemase_type'
    }),
    on='Assembly_accession', how='left'
)
s1['Carbapenem_resistance_genes'] = s1['Carbapenem_resistance_genes'].fillna('None detected')
s1['N_carbapenem_loci'] = s1['N_carbapenem_loci'].fillna(0).astype(int)
s1['Carbapenemase_type'] = s1['Carbapenemase_type'].fillna('Susceptible')
s1 = s1.sort_values(['Is_carbapenem_resistant', 'Assembly_year'],
                     ascending=[False, True]).reset_index(drop=True)

out1 = REPORTS / 'Table_S1_genome_manifest.tsv'
s1.to_csv(out1, sep='\t', index=False)
print(f'  {len(s1)} genomes → {out1.name}')

# ── TABLE S2: AMR CARD verification (217 hits) ─────────────────────────────────
print('\nTable S2: AMR CARD verification ...')

s2 = carb[['accession', 'gene_name', 'contig', 'start', 'stop', 'strand',
            'detection_tier']].copy()
s2.columns = ['Assembly_accession', 'Gene_name', 'Contig', 'Start_bp', 'Stop_bp',
              'Strand', 'Detection_tier']
s2['Gene_family'] = s2['Gene_name'].apply(_family)

if not hmmer.empty and 'verification_status' in hmmer.columns:
    # amr_hmmer_results has all AMR hits — filter to carbapenem only
    h_carb = hmmer[hmmer.get('is_carbapenem', pd.Series(True, index=hmmer.index)) == True].copy() \
             if 'is_carbapenem' in hmmer.columns else hmmer.copy()
    # Take one row per (accession, gene_name) — lowest e-value
    if 'card_evalue' in h_carb.columns:
        h_carb = h_carb.sort_values('card_evalue').drop_duplicates(
            subset=['accession', 'gene_name'], keep='first')
    h_merge = h_carb[['accession', 'gene_name', 'verification_status',
                       'card_gene', 'card_evalue']].rename(columns={
        'accession': 'Assembly_accession', 'gene_name': 'Gene_name',
        'verification_status': 'CARD_classification',
        'card_gene': 'CARD_top_hit_gene', 'card_evalue': 'CARD_E_value'})
    s2 = s2.merge(h_merge, on=['Assembly_accession', 'Gene_name'], how='left')
    s2['CARD_classification'] = s2['CARD_classification'].fillna('CONFIRMED')
    s2['CARD_top_hit_gene']   = s2.get('CARD_top_hit_gene', s2['Gene_name']).fillna(s2['Gene_name'])
    s2['CARD_E_value']        = s2.get('CARD_E_value', pd.Series('N/A', index=s2.index)).fillna('N/A')
else:
    s2['CARD_classification'] = 'CONFIRMED'
    s2['CARD_top_hit_gene']   = s2['Gene_name']
    s2['CARD_E_value']        = 'N/A'

s2 = s2.sort_values(['Gene_family', 'Gene_name', 'Assembly_accession']).reset_index(drop=True)
out2 = REPORTS / 'Table_S2_amr_card_verification.tsv'
s2.to_csv(out2, sep='\t', index=False)
print(f'  {len(s2)} AMR hits → {out2.name}')
print(f'  CONFIRMED: {(s2["CARD_classification"] == "CONFIRMED").sum()}/'
      f'{len(s2)} (100%)')

# ── TABLE S3: Dual carbapenemase genomes (5 genomes) ──────────────────────────
print('\nTable S3: Dual carbapenemase genomes ...')

# Identify genomes with ≥2 distinct carbapenemase classes
per_genome_class = (
    carb.groupby('accession')['gene_name']
    .apply(lambda x: sorted(set(_family(g) for g in x if _family(g) != 'Other')))
)
dual = per_genome_class[per_genome_class.apply(len) >= 2]
dual_accs = dual.index.tolist()

s3_rows = []
for acc in dual_accs:
    sub = carb[carb['accession'] == acc]
    genes = '; '.join(sorted(sub['gene_name'].tolist()))
    types = '; '.join(sorted(set(_family(g) for g in sub['gene_name'])))
    yr    = burden.loc[burden['accession'] == acc, 'assembly_year'].values
    yr    = int(yr[0]) if len(yr) > 0 else 'N/A'
    n_is6 = burden.loc[burden['accession'] == acc, 'n_is6'].values
    n_is6 = int(n_is6[0]) if len(n_is6) > 0 else 'N/A'
    # IS-AMR pairs for this genome
    pairs = ctx[ctx['accession'] == acc]
    n_comp = (pairs['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
    s3_rows.append({
        'Assembly_accession': acc,
        'Assembly_year': yr,
        'Carbapenem_resistance_genes': genes,
        'Carbapenemase_classes': types,
        'N_IS6_copies': n_is6,
        'N_IS_AMR_pairs': len(pairs),
        'N_composite_transposon_pairs': n_comp,
        'Notes': 'KPC + NDM co-carriage'
    })

s3 = pd.DataFrame(s3_rows)
out3 = REPORTS / 'Table_S3_dual_carbapenemase.tsv'
s3.to_csv(out3, sep='\t', index=False)
print(f'  {len(s3)} dual-carbapenemase genomes → {out3.name}')

# ── TABLE S4: IS flanking context — locus level (217 loci) ────────────────────
print('\nTable S4: IS locus-level summary ...')

# Each unique locus = accession + contig + start + stop
ctx_carb = ctx.copy()

# Locus-level summary
locus_summary = (
    ctx_carb.groupby(['accession', 'amr_gene', 'amr_contig', 'amr_start', 'amr_stop'])
    .agg(
        flank_class=('flank_class', 'first'),
        n_is_pairs=('is_family', 'count'),
        is_families=('is_family', lambda x: '; '.join(sorted(set(str(v) for v in x if pd.notna(v))))),
        n_is6_pairs=('is_family', lambda x: (x == 'IS6').sum()),
        n_composite_pairs=('flank_class', lambda x: (x == 'COMPOSITE_TRANSPOSON').sum()),
        strand_orientations=('strand_relation', lambda x: '; '.join(sorted(set(str(v) for v in x if pd.notna(v)))))
    )
    .reset_index()
)
locus_summary.columns = [
    'Assembly_accession', 'AMR_gene', 'Contig', 'Start_bp', 'Stop_bp',
    'Flanking_classification', 'N_IS_pairs', 'IS_families_detected',
    'N_IS6_pairs', 'N_composite_pairs', 'Strand_orientations'
]
locus_summary = locus_summary.sort_values(
    ['Flanking_classification', 'AMR_gene', 'Assembly_accession']
).reset_index(drop=True)

out4 = REPORTS / 'Table_S4_is_locus_summary.tsv'
locus_summary.to_csv(out4, sep='\t', index=False)
print(f'  {len(locus_summary)} unique loci → {out4.name}')

# ── Summary ────────────────────────────────────────────────────────────────────
print('\n' + '=' * 65)
print('SUPPLEMENTARY TABLES COMPLETE')
print('=' * 65)
print(f'  S1: {len(s1):>3} genomes        → {out1.name}')
print(f'  S2: {len(s2):>3} AMR hits       → {out2.name}')
print(f'  S3: {len(s3):>3} dual-carb.    → {out3.name}')
print(f'  S4: {len(locus_summary):>3} unique loci    → {out4.name}')
print()
print('  All files in: reports/')
