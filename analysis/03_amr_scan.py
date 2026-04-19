"""
Step 3 — AMR gene detection
─────────────────────────────────────────────────────────────────────────────
Two-tier approach (rigour before speed):

  Tier 1 — NCBI AMRFinder TSV (if present)
    NCBI runs AMRFinderPlus on every genome they process. The TSV file has
    curated gene names, exact coordinates, percent identity, coverage.
    This is the gold standard.

  Tier 2 — GFF3 regex scan (fallback for genomes without AMRFinder TSV)
    Search product/gene annotations with validated regex patterns.
    Results are flagged as 'GFF_KEYWORD' confidence.

Output:
  data/processed/amr_hits.tsv   — one row per AMR gene hit
  logs/amr_scan.log

Usage:
  python analysis/03_amr_scan.py
"""
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    AMR_GFF_KEYWORDS, CARBAPENEM_GENES, DATA_PROC, LOGS, MANIFEST,
)

LOGS.mkdir(parents=True, exist_ok=True)
DATA_PROC.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'amr_scan.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('amr_scan')

# Pre-compile GFF keyword patterns
_GFF_PATTERNS = [re.compile(p, re.IGNORECASE) for p in AMR_GFF_KEYWORDS]

# Carbapenem gene names as case-insensitive set for fast lookup
_CARBA_LOWER = {g.lower() for g in CARBAPENEM_GENES}


def _is_carbapenem_gene(gene_name: str) -> bool:
    """True if gene_name matches any known carbapenem resistance gene."""
    n = gene_name.lower().strip()
    if n in _CARBA_LOWER:
        return True
    # Prefix match: handles bla-NDM-1 → NDM-1, blaNDM → NDM
    for target in _CARBA_LOWER:
        if n.endswith(target) or n.startswith('bla' + target):
            return True
    return False


# ── Tier 1: parse NCBI AMRFinder TSV ─────────────────────────────────────────

def parse_amrfinder_tsv(tsv_path: Path, accession: str) -> list[dict]:
    """
    Parse *_amrfinderplus.tsv. AMRFinder columns (NCBI format, 2023+):
    Protein identifier, Gene symbol, Sequence name, Scope, Element type,
    Element subtype, Class, Subclass, Method, Target length, Reference seq length,
    % Coverage of reference sequence, % Identity to reference sequence,
    Alignment length, Accession of closest sequence, Name of closest sequence,
    HMM id, HMM description, Contig id, Start, Stop, Strand
    """
    hits = []
    try:
        df = pd.read_csv(tsv_path, sep='\t', dtype=str, comment='#')
    except Exception as e:
        log.error(f'{accession}: AMRFinder TSV parse error: {e}')
        return hits

    # Normalise column names
    df.columns = [c.strip() for c in df.columns]

    for _, row in df.iterrows():
        gene_sym = str(row.get('Gene symbol', '')).strip()
        elem_type = str(row.get('Element type', '')).strip().upper()
        elem_class= str(row.get('Class', '')).strip().upper()
        subclass  = str(row.get('Subclass', '')).strip().upper()
        contig    = str(row.get('Contig id', '')).strip()
        start_str = str(row.get('Start', '')).strip()
        stop_str  = str(row.get('Stop', '')).strip()
        strand    = str(row.get('Strand', '+')).strip()
        pct_id    = row.get('% Identity to reference sequence', '')
        pct_cov   = row.get('% Coverage of reference sequence', '')
        method    = str(row.get('Method', '')).strip()

        # Include only AMR elements (not stress, virulence)
        if elem_type not in ('AMR', 'POINT'):
            continue

        # Parse coordinates
        try:
            start = int(start_str) - 1   # convert to 0-based
            stop  = int(stop_str)
        except ValueError:
            start = stop = -1

        hit = {
            'accession':    accession,
            'gene_name':    gene_sym,
            'element_type': elem_type,
            'drug_class':   elem_class,
            'drug_subclass':subclass,
            'contig':       contig,
            'start':        start,
            'stop':         stop,
            'strand':       strand,
            'pct_identity': pct_id,
            'pct_coverage': pct_cov,
            'method':       method,
            'detection_tier':'AMRFINDER',
            'is_carbapenem': _is_carbapenem_gene(gene_sym)
                             or 'CARBAPENEM' in subclass
                             or 'BETA-LACTAM' in elem_class,
        }
        hits.append(hit)

    log.debug(f'{accession}: AMRFinder → {len(hits)} AMR hits')
    return hits


# ── Tier 2: GFF3 keyword scan (fallback) ─────────────────────────────────────

_ATTRS_RE = re.compile(r'(\w+)=([^;]+)')

_GENE_ALLELE_RE = re.compile(
    r'\b(NDM|KPC|OXA|IMP|VIM|CTX-M|TEM|SHV|CMY|ACC|DHA|FOX|MOX|ACT|MIR|EBC|CMV)\b'
    r'[-_]?(\d+[a-z]?)?', re.IGNORECASE
)


def _canonical_gene(gene_attr: str, product_attr: str) -> str:
    """Extract short canonical AMR gene name like KPC-2, NDM-5, IMP-4."""
    # blaKPC-2 → KPC-2
    m = re.search(r'bla([A-Z]{2,5}-\d+)', gene_attr, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    # blaKPC (no allele) → KPC
    m = re.search(r'bla([A-Z]{2,5})\b', gene_attr, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    # Extract allele from product: "...KPC-2..." → "KPC-2"
    m = _GENE_ALLELE_RE.search(product_attr)
    if m:
        base = m.group(1).upper()
        num  = m.group(2) or ''
        return f'{base}-{num}' if num else base
    return (gene_attr or product_attr[:40]).strip()


def _parse_gff_attrs(attr_str: str) -> dict:
    return dict(m.groups() for m in _ATTRS_RE.finditer(attr_str))


def parse_gff_keywords(gff_gz: Path, accession: str) -> list[dict]:
    """
    Scan GFF3 for AMR gene keywords. Captures CDS / gene features only.
    Returns list of hit dicts.
    """
    hits = []
    try:
        with gzip.open(gff_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                feat_type = cols[2]
                if feat_type not in ('CDS', 'gene'):
                    continue

                attrs   = _parse_gff_attrs(cols[8])
                product = attrs.get('product', '')
                gene    = attrs.get('gene', '')
                note    = attrs.get('Note', '')
                combined = f'{product} {gene} {note}'

                matched_pattern = None
                for pat in _GFF_PATTERNS:
                    if pat.search(combined):
                        matched_pattern = pat.pattern
                        break

                if not matched_pattern:
                    continue

                # Prefer short canonical name; extract allele from product if needed
                gene_name = _canonical_gene(gene, product)

                try:
                    start = int(cols[3]) - 1
                    stop  = int(cols[4])
                except ValueError:
                    start = stop = -1

                hit = {
                    'accession':    accession,
                    'gene_name':    gene_name,
                    'element_type': 'AMR',
                    'drug_class':   'BETA-LACTAM',
                    'drug_subclass':'CARBAPENEM',
                    'contig':       cols[0],
                    'start':        start,
                    'stop':         stop,
                    'strand':       cols[6],
                    'pct_identity': '',
                    'pct_coverage': '',
                    'method':       'GFF_KEYWORD',
                    'detection_tier':'GFF_KEYWORD',
                    'is_carbapenem': _is_carbapenem_gene(gene_name)
                                     or 'carbapenem' in combined.lower()
                                     or 'NDM' in combined
                                     or 'KPC' in combined,
                }
                hits.append(hit)

    except Exception as e:
        log.error(f'{accession}: GFF scan error: {e}')

    log.debug(f'{accession}: GFF keyword → {len(hits)} AMR hits')
    return hits


# ── main ──────────────────────────────────────────────────────────────────────

def scan_genome(row: pd.Series) -> list[dict]:
    """Scan one genome. Tier 1 if AMRFinder available, else Tier 2."""
    accession = row['accession']
    local_dir = Path(row['local_dir'])
    has_amrf  = str(row.get('has_amrfinder', '')).lower() in ('true','1','yes')

    if has_amrf:
        amrf_files = list(local_dir.glob('*_amrfinderplus.tsv'))
        if amrf_files:
            return parse_amrfinder_tsv(amrf_files[0], accession)

    # Fallback to GFF keyword scan
    gff_files = list(local_dir.glob('*_genomic.gff.gz'))
    if gff_files:
        hits = parse_gff_keywords(gff_files[0], accession)
        if hits:
            log.warning(
                f'{accession}: no AMRFinder TSV, using GFF_KEYWORD '
                f'(lower confidence)')
        return hits

    log.error(f'{accession}: no GFF or AMRFinder file found')
    return []


def main():
    if not MANIFEST.exists():
        log.error(f'Manifest not found: {MANIFEST}')
        log.error('Run 02_validate.py first.')
        sys.exit(1)

    manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str)
    log.info(f'Scanning {len(manifest):,} validated genomes for AMR genes')

    all_hits = []
    tier_counts = {'AMRFINDER': 0, 'GFF_KEYWORD': 0, 'NONE': 0}

    for _, row in manifest.iterrows():
        hits = scan_genome(row)
        if hits:
            tier = hits[0].get('detection_tier', 'UNKNOWN')
            tier_counts[tier] = tier_counts.get(tier, 0) + 1
        else:
            tier_counts['NONE'] += 1
        all_hits.extend(hits)

    if not all_hits:
        log.warning('No AMR hits found across all genomes.')
        pd.DataFrame(columns=['accession','gene_name','element_type',
                               'drug_class','drug_subclass','contig',
                               'start','stop','strand','pct_identity',
                               'pct_coverage','method','detection_tier',
                               'is_carbapenem']
                     ).to_csv(DATA_PROC / 'amr_hits.tsv', sep='\t', index=False)
        return

    hits_df = pd.DataFrame(all_hits)

    # Deduplicate: same accession+contig+start+stop may appear if CDS and gene
    # features both match (GFF3 convention where gene and CDS share coordinates).
    # Keep the row with the more informative gene_name (shorter = more canonical).
    hits_df['_locus'] = (hits_df['accession'] + '|' + hits_df['contig'].astype(str)
                         + '|' + hits_df['start'].astype(str)
                         + '|' + hits_df['stop'].astype(str))
    hits_df = (hits_df
               .sort_values('gene_name', key=lambda s: s.str.len())
               .drop_duplicates(subset='_locus', keep='first')
               .drop(columns='_locus')
               .reset_index(drop=True))

    out_path = DATA_PROC / 'amr_hits.tsv'
    hits_df.to_csv(out_path, sep='\t', index=False)

    n_carba  = hits_df['is_carbapenem'].apply(lambda x: str(x).lower() == 'true').sum()
    n_genomes_with_carba = hits_df[
        hits_df['is_carbapenem'].apply(lambda x: str(x).lower() == 'true')
    ]['accession'].nunique()

    log.info('═' * 60)
    log.info(f'AMR scan complete:')
    log.info(f'  Total AMR hits        : {len(hits_df):,}')
    log.info(f'  Carbapenem-resistance : {n_carba:,} hits in {n_genomes_with_carba} genomes')
    log.info(f'  Detection tier breakdown:')
    for k, v in tier_counts.items():
        log.info(f'    {k:<15}: {v:,} genomes')
    log.info(f'  Output → {out_path}')


if __name__ == '__main__':
    main()
