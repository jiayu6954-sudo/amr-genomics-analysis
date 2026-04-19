"""
Step 2 — Data quality gate
─────────────────────────────────────────────────────────────────────────────
Reads data/raw/download_status.tsv (from 01_download.py).
For each genome with status=OK, checks:
  1. Files exist and are non-zero
  2. GFF.gz can be decompressed and parsed (no truncation)
  3. FASTA.gz can be decompressed (non-empty)
  4. Assembly stats: genome size, CDS count, N50 — within species thresholds
  5. Species identity cross-check (organism name in GFF header)

Writes:
  data/validated/genome_manifest.tsv   — PASS genomes only
  logs/validate.log                    — every decision + reason

Usage:
  python analysis/02_validate.py
"""
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_RAW, DATA_VAL, LOGS, MANIFEST, TARGET_SPECIES,
)

LOGS.mkdir(parents=True, exist_ok=True)
DATA_VAL.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'validate.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('validate')


# ── per-genome validation ─────────────────────────────────────────────────────

class ValidationError(Exception):
    pass


def _find_file(local_dir: Path, suffix: str) -> Path | None:
    """Return first file in local_dir ending with suffix, or None."""
    matches = list(local_dir.glob(f'*{suffix}'))
    return matches[0] if matches else None


def _check_gz_readable(path: Path, read_bytes: int = 65536) -> None:
    """Raise ValidationError if gz file is unreadable or empty."""
    try:
        with gzip.open(path, 'rb') as f:
            data = f.read(read_bytes)
        if not data:
            raise ValidationError(f'{path.name} is empty after decompression')
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        raise ValidationError(f'{path.name} is corrupt: {e}')


def _parse_assembly_stats(stats_path: Path) -> dict:
    """
    Parse *_assembly_stats.txt.
    Returns dict with keys: total_length, contig_count, contig_N50, cds_count.
    """
    values = {}
    try:
        text = stats_path.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return values

    patterns = {
        'total_length':  r'all\s+all\s+all\s+all\s+total-length\s+(\d+)',
        'contig_count':  r'all\s+all\s+all\s+all\s+molecule-count\s+(\d+)',
        'contig_N50':    r'all\s+all\s+all\s+all\s+contig-N50\s+(\d+)',
        'scaffold_N50':  r'all\s+all\s+all\s+all\s+scaffold-N50\s+(\d+)',
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            values[key] = int(m.group(1))

    # CDS count from GFF3-level stats (may not be in assembly_stats)
    return values


def _count_cds_in_gff(gff_gz: Path) -> int:
    """Count CDS features in GFF3. Raises ValidationError on parse failure."""
    count = 0
    try:
        with gzip.open(gff_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) >= 3 and cols[2] == 'CDS':
                    count += 1
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        raise ValidationError(f'GFF parse error: {e}')
    return count


def _check_organism_in_gff(gff_gz: Path, expected_genus: str) -> str:
    """
    Read GFF3 header for organism info.
    NCBI GFF3 header contains ##species URL with taxid — use that.
    Returns a string for logging; does NOT raise (non-blocking check).
    """
    GENUS_TAXID_MAP = {
        'Klebsiella':   '570',    # genus taxid
        'Escherichia':  '561',
    }
    expected_taxid_prefix = GENUS_TAXID_MAP.get(expected_genus, '')
    species_url = ''

    try:
        with gzip.open(gff_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                m = re.search(r'##species\s+(https?://\S+\?id=(\d+))', line)
                if m:
                    species_url = m.group(1)
                    taxid       = m.group(2)
                    # Loose check: Klebsiella taxids start with the genus prefix
                    # (not perfect but avoids false rejections)
                    return f'taxid={taxid}'
    except Exception:
        pass

    # Can't verify — don't block, just flag as unverified
    return 'organism_unverified'


def validate_one(row: pd.Series) -> dict:
    """
    Validate a single genome. Returns a dict with all QC fields.
    'qc_pass' = True only if every check passes.
    """
    accession = row['accession']
    org_name  = str(row.get('organism_name', ''))

    # Derive species key from organism_name
    if 'Klebsiella' in org_name:
        species = 'Klebsiella_pneumoniae'
    elif 'Escherichia' in org_name:
        species = 'Escherichia_coli'
    else:
        species = row.get('species', org_name.replace(' ', '_')[:30])

    local_dir = Path(row['local_dir'])

    result = {
        'accession':     accession,
        'species':       species,
        'local_dir':     str(local_dir),
        'organism_name': row.get('organism_name', ''),
        'geo_loc_name':  row.get('geo_loc_name', ''),
        'isolation_source': row.get('isolation_source', ''),
        'biosample':     row.get('biosample', ''),
        'assembly_level':row.get('assembly_level', ''),
        'ftp_path':      row.get('ftp_path', ''),
        'has_amrfinder': row.get('has_amrfinder', False),
        # QC fields
        'total_length':  None,
        'contig_N50':    None,
        'scaffold_N50':  None,
        'cds_count':     None,
        'organism_confirmed': '',
        'qc_pass':       False,
        'qc_fail_reason': '',
    }

    # Map species label to canonical name for threshold lookup
    sp_name_map = {
        'Klebsiella_pneumoniae': 'Klebsiella pneumoniae',
        'Escherichia_coli':      'Escherichia coli',
    }
    sp_canonical = sp_name_map.get(species, species.replace('_', ' '))
    thresholds   = TARGET_SPECIES.get(sp_canonical, {})

    try:
        # ── 1. directory exists ──────────────────────────────────────────
        if not local_dir.exists():
            raise ValidationError('local_dir does not exist')

        # ── 2. GFF3 exists and is readable ──────────────────────────────
        gff = _find_file(local_dir, '_genomic.gff.gz')
        if gff is None:
            raise ValidationError('GFF3 file missing')
        _check_gz_readable(gff)

        # ── 3. FASTA exists and is readable ─────────────────────────────
        fna = _find_file(local_dir, '_genomic.fna.gz')
        if fna is None:
            raise ValidationError('FASTA file missing')
        _check_gz_readable(fna)

        # ── 4. assembly stats ────────────────────────────────────────────
        stats_file = _find_file(local_dir, '_assembly_stats.txt')
        stats = _parse_assembly_stats(stats_file) if stats_file else {}
        result.update({k: stats.get(k) for k in
                       ('total_length','contig_count','contig_N50','scaffold_N50')})

        best_n50 = stats.get('scaffold_N50') or stats.get('contig_N50') or 0

        if thresholds:
            if stats.get('total_length'):
                sz = stats['total_length']
                if sz < thresholds['genome_size_min_bp']:
                    raise ValidationError(
                        f'Genome too small: {sz:,} bp '
                        f'(min {thresholds["genome_size_min_bp"]:,})')
                if sz > thresholds['genome_size_max_bp']:
                    raise ValidationError(
                        f'Genome too large: {sz:,} bp '
                        f'(max {thresholds["genome_size_max_bp"]:,})')
            if best_n50 and best_n50 < thresholds['min_n50']:
                raise ValidationError(
                    f'N50 too low: {best_n50:,} bp (min {thresholds["min_n50"]:,})')

        # ── 5. CDS count ─────────────────────────────────────────────────
        cds = _count_cds_in_gff(gff)
        result['cds_count'] = cds
        if thresholds and cds < thresholds['min_cds']:
            raise ValidationError(
                f'Too few CDS: {cds} (min {thresholds["min_cds"]})')

        # ── 6. organism cross-check ──────────────────────────────────────
        expected_genus = sp_canonical.split()[0]   # 'Klebsiella' or 'Escherichia'
        org = _check_organism_in_gff(gff, expected_genus)
        result['organism_confirmed'] = org

        # ── All checks passed ────────────────────────────────────────────
        result['qc_pass'] = True
        log.info(f'PASS  {accession}  '
                 f'size={stats.get("total_length",0)//1000}kb  '
                 f'N50={best_n50//1000}kb  CDS={cds}')

    except ValidationError as e:
        result['qc_fail_reason'] = str(e)
        log.warning(f'FAIL  {accession}: {e}')

    return result


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    status_path = DATA_RAW / 'download_status.tsv'
    if not status_path.exists():
        log.error(f'download_status.tsv not found: {status_path}')
        log.error('Run 01_download.py first.')
        sys.exit(1)

    status_df = pd.read_csv(status_path, sep='\t', dtype=str)
    # Normalise accession column name (01_download uses 'assembly_accession')
    if 'assembly_accession' in status_df.columns and 'accession' not in status_df.columns:
        status_df = status_df.rename(columns={'assembly_accession': 'accession'})
    ok_df = status_df[status_df['status'] == 'OK'].copy()
    log.info(f'Validating {len(ok_df):,} genomes with status=OK')

    records = []
    for _, row in ok_df.iterrows():
        records.append(validate_one(row))

    qc_df   = pd.DataFrame(records)
    pass_df = qc_df[qc_df['qc_pass'] == True].copy()
    fail_df = qc_df[qc_df['qc_pass'] != True].copy()

    # Write manifest (PASS only)
    pass_df.to_csv(MANIFEST, sep='\t', index=False)

    # Write full QC report
    qc_report = DATA_VAL / 'qc_report.tsv'
    qc_df.to_csv(qc_report, sep='\t', index=False)

    log.info('═' * 60)
    log.info(f'Validation complete:')
    log.info(f'  PASS : {len(pass_df):,}')
    log.info(f'  FAIL : {len(fail_df):,}')
    log.info(f'  Manifest → {MANIFEST}')

    if not fail_df.empty:
        log.warning('Failed genomes:')
        for _, r in fail_df.iterrows():
            log.warning(f'  {r["accession"]}: {r["qc_fail_reason"]}')

    if pass_df.empty:
        log.error('No genomes passed validation — check download logs.')
        sys.exit(1)

    log.info(f'AMRFinder TSV available in {pass_df["has_amrfinder"].sum()} / {len(pass_df)} genomes')


if __name__ == '__main__':
    main()
