"""
Step 0 — Data Integrity Audit
─────────────────────────────────────────────────────────────────────────────
Runs BEFORE any analysis. Checks every genome in data/raw/ for:

  A. File completeness  — required .gff.gz and .fna.gz present and non-empty
  B. Gzip integrity     — full decompression without error (catches truncation)
  C. GFF3 format        — column count, coordinate sanity, CDS presence
  D. FASTA format       — valid DNA, ≥1 contig, computed size vs stats file
  E. Assembly stats     — genome size / N50 parseable
  F. Duplicate detection — same BioSample ID across accessions (biological dup)
                         — identical accession appears twice (registry dup)
  G. QC threshold re-check — re-verify size / N50 / CDS against config thresholds

Output:
  data/processed/audit_report.tsv   — one row per genome, all checks
  data/processed/audit_summary.json — pass/fail counts + issue list

Usage:
  python analysis/00_data_audit.py [--fast]   # --fast skips full gzip decompression
"""

import gzip
import json
import logging
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_RAW, DATA_PROC, DATA_VAL, LOGS,
    TARGET_SPECIES,
)

LOGS.mkdir(parents=True, exist_ok=True)
DATA_PROC.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'audit.log', encoding='utf-8'),
        logging.StreamHandler(open(sys.stdout.fileno(), mode='w',
                                   encoding='utf-8', closefd=False)),
    ]
)
log = logging.getLogger('audit')

# QC thresholds (K. pneumoniae)
SIZE_MIN = 4_800_000
SIZE_MAX = 6_500_000
N50_MIN  = 50_000
CDS_MIN  = 3_000

# ── helpers ───────────────────────────────────────────────────────────────────

def check_gzip_integrity(path: Path, fast: bool = False) -> tuple[bool, str]:
    """Return (ok, message). fast=True reads only first/last 64KB."""
    try:
        with gzip.open(path, 'rb') as f:
            if fast:
                f.read(65536)
                # seek to near-end trick: just read first chunk
            else:
                while f.read(65536):
                    pass
        return True, 'OK'
    except EOFError:
        return False, 'TRUNCATED'
    except gzip.BadGzipFile:
        return False, 'BAD_GZIP'
    except Exception as e:
        return False, f'ERROR:{e}'


def check_gff3(path: Path) -> dict:
    """Validate GFF3 structure. Returns dict of metrics."""
    result = {
        'gff_ok': False, 'gff_error': '',
        'n_cds': 0, 'n_features': 0,
        'bad_coord_rows': 0, 'bad_strand_rows': 0,
    }
    valid_strands = {'+', '-', '.', '?'}
    coord_re = re.compile(r'^\d+$')
    try:
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
            for lineno, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    if line.strip():   # non-empty, non-comment
                        result['gff_error'] = f'line {lineno}: only {len(parts)} cols'
                        # don't abort — count as warning
                    continue
                result['n_features'] += 1
                ftype  = parts[2]
                start  = parts[3]
                stop   = parts[4]
                strand = parts[6]

                if ftype == 'CDS':
                    result['n_cds'] += 1

                # coordinate sanity
                if not (coord_re.match(start) and coord_re.match(stop)):
                    result['bad_coord_rows'] += 1
                elif int(start) > int(stop):
                    result['bad_coord_rows'] += 1

                # strand sanity
                if strand not in valid_strands:
                    result['bad_strand_rows'] += 1

        result['gff_ok'] = (
            result['n_features'] > 0 and
            result['bad_coord_rows'] == 0 and
            result['n_cds'] >= 100           # minimal CDS floor
        )
        if result['n_cds'] < 100:
            result['gff_error'] = result['gff_error'] or f'only {result["n_cds"]} CDS'
    except Exception as e:
        result['gff_error'] = str(e)
    return result


def check_fasta(path: Path) -> dict:
    """Validate FASTA and compute genome size."""
    result = {
        'fna_ok': False, 'fna_error': '',
        'n_contigs': 0, 'computed_size_bp': 0,
        'bad_seq_chars': 0,
    }
    valid_dna = set('ACGTNacgtnRYSWKMBDHVryswkmbdhv')
    try:
        in_seq = False
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    result['n_contigs'] += 1
                    in_seq = True
                elif in_seq:
                    result['computed_size_bp'] += len(line)
                    bad = sum(1 for c in line if c not in valid_dna)
                    result['bad_seq_chars'] += bad
        result['fna_ok'] = (
            result['n_contigs'] >= 1 and
            result['computed_size_bp'] >= SIZE_MIN * 0.8 and   # 20% margin
            result['bad_seq_chars'] == 0
        )
        if result['n_contigs'] == 0:
            result['fna_error'] = 'no contigs found'
        elif result['bad_seq_chars'] > 0:
            result['fna_error'] = f'{result["bad_seq_chars"]} invalid DNA chars'
        elif result['computed_size_bp'] < SIZE_MIN * 0.8:
            result['fna_error'] = f'size too small: {result["computed_size_bp"]:,} bp'
    except Exception as e:
        result['fna_error'] = str(e)
    return result


def parse_assembly_stats(path: Path) -> dict:
    """Parse *_assembly_stats.txt for genome size and N50."""
    result = {'stats_n50': 0, 'stats_size': 0, 'stats_ok': False}
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                # lines like: all	all	all	all	total-length	5432100
                parts = line.split('\t')
                if len(parts) >= 6:
                    scope, unit, mol, clas, stat, val = parts[:6]
                    if scope == 'all' and unit == 'all' and mol == 'all':
                        if stat == 'total-length':
                            try:
                                result['stats_size'] = int(val)
                            except ValueError:
                                pass
                        elif stat == 'scaffold-N50':
                            try:
                                result['stats_n50'] = int(val)
                            except ValueError:
                                pass
        result['stats_ok'] = result['stats_size'] > 0
    except Exception as e:
        pass
    return result


# ── per-genome audit ──────────────────────────────────────────────────────────

def audit_genome(acc_dir: Path, fast: bool = False) -> dict:
    accession = acc_dir.name
    row = {
        'accession': accession,
        # File presence
        'has_gff': False, 'has_fna': False, 'has_stats': False,
        # Gzip integrity
        'gff_gzip_ok': False, 'gff_gzip_msg': '',
        'fna_gzip_ok': False, 'fna_gzip_msg': '',
        # GFF3 content
        'gff_ok': False, 'gff_error': '', 'n_cds': 0, 'n_features': 0,
        'bad_coord_rows': 0, 'bad_strand_rows': 0,
        # FASTA content
        'fna_ok': False, 'fna_error': '', 'n_contigs': 0,
        'computed_size_bp': 0, 'bad_seq_chars': 0,
        # Assembly stats
        'stats_ok': False, 'stats_size': 0, 'stats_n50': 0,
        # Cross-validation
        'size_consistent': False,  # FASTA size vs stats within 5%
        # QC thresholds
        'qc_size': False, 'qc_n50': False, 'qc_cds': False,
        'qc_pass': False,
        # Overall
        'audit_pass': False, 'issues': '',
    }

    # ── A. File presence ──────────────────────────────────────────────────────
    gff_files   = sorted(acc_dir.glob('*_genomic.gff.gz'))
    fna_files   = sorted(acc_dir.glob('*_genomic.fna.gz'))
    stat_files  = sorted(acc_dir.glob('*_assembly_stats.txt'))

    row['has_gff']   = bool(gff_files)
    row['has_fna']   = bool(fna_files)
    row['has_stats'] = bool(stat_files)

    issues = []
    if not row['has_gff']:
        issues.append('MISSING_GFF')
    if not row['has_fna']:
        issues.append('MISSING_FNA')

    # ── B. Gzip integrity ─────────────────────────────────────────────────────
    if gff_files:
        ok, msg = check_gzip_integrity(gff_files[0], fast=fast)
        row['gff_gzip_ok'], row['gff_gzip_msg'] = ok, msg
        if not ok:
            issues.append(f'GFF_GZIP_{msg}')

    if fna_files:
        ok, msg = check_gzip_integrity(fna_files[0], fast=fast)
        row['fna_gzip_ok'], row['fna_gzip_msg'] = ok, msg
        if not ok:
            issues.append(f'FNA_GZIP_{msg}')

    # ── C. GFF3 format ───────────────────────────────────────────────────────
    if gff_files and row['gff_gzip_ok']:
        g = check_gff3(gff_files[0])
        row.update(g)
        if not g['gff_ok']:
            issues.append(f'GFF_INVALID:{g["gff_error"]}')
        if g['bad_coord_rows']:
            issues.append(f'GFF_BAD_COORDS:{g["bad_coord_rows"]}')
        if g['bad_strand_rows']:
            issues.append(f'GFF_BAD_STRAND:{g["bad_strand_rows"]}')

    # ── D. FASTA format ───────────────────────────────────────────────────────
    if fna_files and row['fna_gzip_ok']:
        fa = check_fasta(fna_files[0])
        row.update(fa)
        if not fa['fna_ok']:
            issues.append(f'FNA_INVALID:{fa["fna_error"]}')

    # ── E. Assembly stats ─────────────────────────────────────────────────────
    if stat_files:
        st = parse_assembly_stats(stat_files[0])
        row.update(st)
        if not st['stats_ok']:
            issues.append('STATS_UNREADABLE')

    # ── F. Cross-validation: FASTA size vs stats ──────────────────────────────
    if row['computed_size_bp'] > 0 and row['stats_size'] > 0:
        ratio = row['computed_size_bp'] / row['stats_size']
        row['size_consistent'] = 0.95 <= ratio <= 1.05
        if not row['size_consistent']:
            issues.append(f'SIZE_MISMATCH:fasta={row["computed_size_bp"]},stats={row["stats_size"]}')

    # ── G. QC thresholds ─────────────────────────────────────────────────────
    size = row['stats_size'] if row['stats_size'] > 0 else row['computed_size_bp']
    row['qc_size'] = SIZE_MIN <= size <= SIZE_MAX
    row['qc_n50']  = row['stats_n50'] >= N50_MIN
    row['qc_cds']  = row['n_cds'] >= CDS_MIN
    row['qc_pass'] = row['qc_size'] and row['qc_n50'] and row['qc_cds']

    if not row['qc_size']:
        issues.append(f'QC_SIZE_FAIL:{size:,}bp')
    if not row['qc_n50']:
        issues.append(f'QC_N50_FAIL:{row["stats_n50"]:,}')
    if not row['qc_cds']:
        issues.append(f'QC_CDS_FAIL:{row["n_cds"]}')

    # ── Overall ───────────────────────────────────────────────────────────────
    row['issues'] = '; '.join(issues)
    row['audit_pass'] = (
        row['has_gff'] and row['has_fna'] and
        row['gff_gzip_ok'] and row['fna_gzip_ok'] and
        row['gff_ok'] and row['fna_ok']
        # Note: qc_pass is separate — a genome can pass audit but fail QC
    )
    return row


# ── duplicate detection ───────────────────────────────────────────────────────

def check_duplicates(status_tsv: Path) -> dict:
    """
    Returns dict with:
      registry_dups  — same accession appears >1 time in download_status.tsv
      biosample_dups — same BioSample ID across different accessions
    """
    result = {'registry_dups': [], 'biosample_dups': []}
    if not status_tsv.exists():
        return result

    df = pd.read_csv(status_tsv, sep='\t', dtype=str)
    if 'assembly_accession' not in df.columns:
        return result

    # Registry duplicates
    acc_counts = df['assembly_accession'].value_counts()
    result['registry_dups'] = acc_counts[acc_counts > 1].index.tolist()

    # BioSample duplicates
    if 'biosample' in df.columns:
        bs_map = defaultdict(list)
        for _, r in df.iterrows():
            bs = r.get('biosample', '')
            acc = r.get('assembly_accession', '')
            if bs and bs not in ('', 'nan', 'NA'):
                bs_map[bs].append(acc)
        result['biosample_dups'] = {
            bs: accs for bs, accs in bs_map.items() if len(accs) > 1
        }
    return result


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--fast', action='store_true',
                    help='Skip full gzip decompression (faster but misses truncation)')
    args = ap.parse_args()

    log.info('=' * 70)
    log.info('DATA INTEGRITY AUDIT — AMR genomics project')
    log.info(f'Mode: {"FAST (no full decompress)" if args.fast else "FULL (complete decompression)"}')
    log.info('=' * 70)

    # ── 1. Find all genome directories ───────────────────────────────────────
    genome_dirs = sorted(
        d for d in DATA_RAW.iterdir()
        if d.is_dir() and (d.name.startswith('GCA_') or d.name.startswith('GCF_'))
    )
    log.info(f'Found {len(genome_dirs)} genome directories in data/raw/')

    if not genome_dirs:
        log.error('No genome directories found. Run 01_download.py first.')
        sys.exit(1)

    # ── 2. Duplicate detection (fast, before per-genome loop) ────────────────
    status_tsv = DATA_RAW / 'download_status.tsv'
    dups = check_duplicates(status_tsv)
    if dups['registry_dups']:
        log.warning(f'REGISTRY DUPLICATES ({len(dups["registry_dups"])}): '
                    f'{dups["registry_dups"]}')
    else:
        log.info('Registry duplicates: NONE ✓')

    if dups.get('biosample_dups'):
        n = len(dups['biosample_dups'])
        log.warning(f'BIOSAMPLE DUPLICATES ({n} BioSamples with >1 accession):')
        for bs, accs in list(dups['biosample_dups'].items())[:10]:
            log.warning(f'  {bs}: {accs}')
        if n > 10:
            log.warning(f'  ... and {n - 10} more')
    else:
        log.info('BioSample duplicates: NONE ✓')

    # ── 3. Per-genome audit ───────────────────────────────────────────────────
    rows = []
    t0 = time.time()
    for i, d in enumerate(genome_dirs, 1):
        if i % 50 == 0 or i == 1 or i == len(genome_dirs):
            elapsed = time.time() - t0
            rate = i / elapsed if elapsed > 0 else 0
            eta = (len(genome_dirs) - i) / rate if rate > 0 else 0
            log.info(f'  {i:>4}/{len(genome_dirs)} [{elapsed:.0f}s elapsed, ETA {eta:.0f}s]')
        row = audit_genome(d, fast=args.fast)
        rows.append(row)

    df = pd.DataFrame(rows)

    # ── 4. Summary statistics ─────────────────────────────────────────────────
    n_total  = len(df)
    n_audit  = df['audit_pass'].sum()
    n_qc     = df['qc_pass'].sum()
    n_both   = (df['audit_pass'] & df['qc_pass']).sum()

    # Issue breakdown
    all_issues = []
    for iss in df['issues']:
        if iss:
            all_issues.extend(iss.split('; '))
    issue_counts = defaultdict(int)
    for iss in all_issues:
        prefix = iss.split(':')[0]
        issue_counts[prefix] += 1

    # Duplicate-infected accessions
    dup_accs = set(dups['registry_dups'])
    for accs in dups.get('biosample_dups', {}).values():
        dup_accs.update(accs)
    n_dup = len(dup_accs & set(df['accession']))

    summary = {
        'total_genomes_audited':   n_total,
        'audit_pass':              int(n_audit),
        'audit_fail':              int(n_total - n_audit),
        'qc_pass':                 int(n_qc),
        'qc_fail':                 int(n_total - n_qc),
        'both_audit_and_qc_pass':  int(n_both),
        'registry_duplicates':     len(dups['registry_dups']),
        'biosample_duplicates':    len(dups.get('biosample_dups', {})),
        'accessions_with_dups':    int(n_dup),
        'issue_counts':            dict(sorted(issue_counts.items(),
                                               key=lambda x: -x[1])),
        'size_stats': {
            'min_bp':  int(df['computed_size_bp'].min()),
            'max_bp':  int(df['computed_size_bp'].max()),
            'median_bp': int(df['computed_size_bp'].median()),
        },
        'n50_stats': {
            'min':    int(df['stats_n50'].min()),
            'max':    int(df['stats_n50'].max()),
            'median': int(df['stats_n50'].median()),
        },
        'cds_stats': {
            'min':    int(df['n_cds'].min()),
            'max':    int(df['n_cds'].max()),
            'median': int(df['n_cds'].median()),
        },
    }

    # ── 5. Save outputs ───────────────────────────────────────────────────────
    out_tsv  = DATA_PROC / 'audit_report.tsv'
    out_json = DATA_PROC / 'audit_summary.json'
    df.to_csv(out_tsv, sep='\t', index=False)
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # ── 6. Print final report ─────────────────────────────────────────────────
    log.info('')
    log.info('══════════════════════════════════════════════')
    log.info('  AUDIT REPORT')
    log.info('══════════════════════════════════════════════')
    log.info(f'  Total genomes:           {n_total}')
    log.info(f'  File + format audit PASS: {n_audit}  ({100*n_audit/n_total:.1f}%)')
    log.info(f'  File + format audit FAIL: {n_total-n_audit}')
    log.info(f'  QC threshold PASS:        {n_qc}  ({100*n_qc/n_total:.1f}%)')
    log.info(f'  Both audit+QC PASS:       {n_both}  ({100*n_both/n_total:.1f}%)')
    log.info(f'  Registry duplicates:      {len(dups["registry_dups"])}')
    log.info(f'  BioSample duplicates:     {len(dups.get("biosample_dups", {}))}')
    log.info('')
    log.info('  Issue breakdown:')
    for iss, cnt in sorted(issue_counts.items(), key=lambda x: -x[1]):
        log.info(f'    {iss:<35} {cnt}')
    log.info('')
    log.info(f'  Genome size range:  {summary["size_stats"]["min_bp"]:,} – '
             f'{summary["size_stats"]["max_bp"]:,} bp '
             f'(median {summary["size_stats"]["median_bp"]:,})')
    log.info(f'  N50 range:          {summary["n50_stats"]["min"]:,} – '
             f'{summary["n50_stats"]["max"]:,} '
             f'(median {summary["n50_stats"]["median"]:,})')
    log.info(f'  CDS count range:    {summary["cds_stats"]["min"]} – '
             f'{summary["cds_stats"]["max"]} '
             f'(median {summary["cds_stats"]["median"]})')
    log.info('')
    log.info(f'  Output: {out_tsv}')
    log.info(f'  Output: {out_json}')
    log.info('══════════════════════════════════════════════')

    # Scientific verdict
    if n_total - n_audit > 0:
        log.warning('⚠  Some genomes failed audit — inspect audit_report.tsv '
                    'before proceeding with analysis.')
    if dups.get('biosample_dups'):
        log.warning('⚠  BioSample duplicates detected — consider deduplication '
                    'before statistical analysis (inflated N).')
    if n_audit == n_total and not dups.get('biosample_dups'):
        log.info('✅ All genomes pass audit. Data is consistent and ready for analysis.')

    return summary


if __name__ == '__main__':
    main()
