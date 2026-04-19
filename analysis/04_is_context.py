"""
Step 4 — IS element genomic context of AMR genes
─────────────────────────────────────────────────────────────────────────────
For each carbapenem resistance gene from amr_hits.tsv:
  1. Find IS elements on the same contig within IS_FLANK_WINDOW_BP
  2. Record distance, orientation, IS family name
  3. Classify flanking structure:
       COMPOSITE_TRANSPOSON   — IS on both sides (classic composite)
       SINGLE_IS_UPSTREAM     — one IS, upstream
       SINGLE_IS_DOWNSTREAM   — one IS, downstream
       IS_COMPLEX             — multiple IS, mixed orientation
       NO_IS                  — no IS element within window

Output:
  data/processed/is_context.tsv      — one row per AMR gene × IS element pair
  data/processed/context_summary.tsv — one row per genome, summary counts
  logs/is_context.log

Usage:
  python analysis/04_is_context.py
"""
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_PROC, IS_FLANK_WINDOW_BP, IS_GFF_KEYWORDS, LOGS, MANIFEST,
)

LOGS.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'is_context.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('is_context')

_IS_PATTERNS = [re.compile(p, re.IGNORECASE) for p in IS_GFF_KEYWORDS]
_ATTRS_RE    = re.compile(r'(\w+)=([^;]+)')


def _parse_attrs(attr_str: str) -> dict:
    return dict(m.groups() for m in _ATTRS_RE.finditer(attr_str))


def _extract_is_family(text: str) -> str:
    """
    Try to extract IS family name from product/gene annotation.
    Examples: 'IS26 transposase' → 'IS26'
              'ISKpn26' → 'ISKpn26'
              'Tn4401 transposase' → 'Tn4401'
    """
    m = re.search(r'\b(Tn\d+)\b', text, re.IGNORECASE)
    if m: return m.group(1)
    m = re.search(r'\b(IS[A-Za-z0-9]+)\b', text, re.IGNORECASE)
    if m: return m.group(1)
    return 'IS_unknown'


def load_is_features(gff_gz: Path, accession: str) -> pd.DataFrame:
    """
    Parse GFF3 and return a DataFrame of IS element features.
    Columns: contig, start, stop, strand, is_family, product
    """
    records = []
    try:
        with gzip.open(gff_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                feat_type = cols[2]
                if feat_type not in ('CDS', 'gene', 'mobile_element',
                                     'repeat_region'):
                    continue
                attrs   = _parse_attrs(cols[8])
                product = attrs.get('product', '')
                gene    = attrs.get('gene', '')
                note    = attrs.get('Note', '')
                combined= f'{product} {gene} {note}'

                matched = any(p.search(combined) for p in _IS_PATTERNS)
                if not matched:
                    continue

                try:
                    start = int(cols[3]) - 1
                    stop  = int(cols[4])
                except ValueError:
                    continue

                records.append({
                    'contig':    cols[0],
                    'start':     start,
                    'stop':      stop,
                    'strand':    cols[6],
                    'is_family': _extract_is_family(combined),
                    'product':   product[:80],
                })

    except Exception as e:
        log.error(f'{accession}: IS feature parse error: {e}')

    return pd.DataFrame(records)


def classify_flank(upstream_count: int, downstream_count: int) -> str:
    """Classify IS element flanking structure."""
    if upstream_count > 0 and downstream_count > 0:
        return 'COMPOSITE_TRANSPOSON'
    elif upstream_count > 0:
        return 'SINGLE_IS_UPSTREAM'
    elif downstream_count > 0:
        return 'SINGLE_IS_DOWNSTREAM'
    else:
        return 'NO_IS'


def analyse_context(
        amr_row: pd.Series,
        is_df: pd.DataFrame,
        window: int = IS_FLANK_WINDOW_BP
) -> list[dict]:
    """
    Find IS elements flanking one AMR gene.
    Returns list of dicts (one per IS element found in window).
    Returns [{'no_is': True}] placeholder if none found.
    """
    contig = str(amr_row.get('contig', ''))
    try:
        amr_start = int(amr_row.get('start', -1))
        amr_stop  = int(amr_row.get('stop', -1))
    except (ValueError, TypeError):
        return []

    if amr_start < 0 or not contig:
        return []

    # Filter IS elements on same contig
    is_contig = is_df[is_df['contig'] == contig]
    if is_contig.empty:
        return [{
            'accession':      amr_row['accession'],
            'amr_gene':       amr_row['gene_name'],
            'amr_contig':     contig,
            'amr_start':      amr_start,
            'amr_stop':       amr_stop,
            'amr_strand':     str(amr_row.get('strand', '?')),
            'is_family':      '',
            'is_start':       -1,
            'is_stop':        -1,
            'is_strand':      '',
            'strand_relation': '',
            'distance_bp':    -1,
            'position':       'NO_IS',
            'flank_class':    'NO_IS',
        }]

    # Flanking IS elements within window
    results = []
    upstream_count = downstream_count = 0

    for _, is_row in is_contig.iterrows():
        is_start = int(is_row['start'])
        is_stop  = int(is_row['stop'])

        # Distance: gap between IS and AMR gene (0 if overlapping)
        gap_right = max(0, is_start - amr_stop)   # IS downstream
        gap_left  = max(0, amr_start - is_stop)   # IS upstream
        distance  = min(gap_right, gap_left)

        if distance > window:
            continue

        position = 'UPSTREAM' if gap_left <= gap_right else 'DOWNSTREAM'
        if position == 'UPSTREAM':
            upstream_count += 1
        else:
            downstream_count += 1

        # is_strand relative to amr_strand: 'same' / 'opposite' / 'unknown'
        amr_strand = str(amr_row.get('strand', '?'))
        is_strand   = str(is_row.get('strand', '?'))
        if amr_strand in ('+', '-') and is_strand in ('+', '-'):
            strand_relation = 'same' if amr_strand == is_strand else 'opposite'
        else:
            strand_relation = 'unknown'

        results.append({
            'accession':      amr_row['accession'],
            'amr_gene':       amr_row['gene_name'],
            'amr_contig':     contig,
            'amr_start':      amr_start,
            'amr_stop':       amr_stop,
            'amr_strand':     amr_strand,
            'is_family':      is_row['is_family'],
            'is_start':       is_start,
            'is_stop':        is_stop,
            'is_strand':      is_strand,
            'strand_relation': strand_relation,
            'distance_bp':    distance,
            'position':       position,
            'flank_class':    '',   # filled below
        })

    if not results:
        results.append({
            'accession':      amr_row['accession'],
            'amr_gene':       amr_row['gene_name'],
            'amr_contig':     contig,
            'amr_start':      amr_start,
            'amr_stop':       amr_stop,
            'amr_strand':     str(amr_row.get('strand', '?')),
            'is_family':      '',
            'is_start':       -1,
            'is_stop':        -1,
            'is_strand':      '',
            'strand_relation': '',
            'distance_bp':    -1,
            'position':       'NO_IS',
            'flank_class':    'NO_IS',
        })
    else:
        flank_class = classify_flank(upstream_count, downstream_count)
        for r in results:
            r['flank_class'] = flank_class

    return results


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    amr_path = DATA_PROC / 'amr_hits.tsv'
    if not amr_path.exists():
        log.error(f'amr_hits.tsv not found: {amr_path}')
        log.error('Run 03_amr_scan.py first.')
        sys.exit(1)
    if not MANIFEST.exists():
        log.error(f'Manifest not found: {MANIFEST}')
        sys.exit(1)

    amr_df   = pd.read_csv(amr_path, sep='\t', dtype=str)
    manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str)

    # Focus on carbapenem resistance genes
    carba_df = amr_df[
        amr_df['is_carbapenem'].apply(lambda x: str(x).lower() == 'true')
    ].copy()
    log.info(f'Analysing IS context for {len(carba_df):,} carbapenem resistance hits '
             f'across {carba_df["accession"].nunique()} genomes')

    all_context = []
    summary_rows = []

    for accession in carba_df['accession'].unique():
        man_row = manifest[manifest['accession'] == accession]
        if man_row.empty:
            log.warning(f'{accession}: not in manifest, skipping')
            continue

        local_dir = Path(man_row.iloc[0]['local_dir'])
        gff_files = list(local_dir.glob('*_genomic.gff.gz'))
        if not gff_files:
            log.warning(f'{accession}: no GFF found, skipping IS context')
            continue

        is_df  = load_is_features(gff_files[0], accession)
        subset = carba_df[carba_df['accession'] == accession]
        genome_contexts = []

        for _, amr_row in subset.iterrows():
            ctx = analyse_context(amr_row, is_df)
            genome_contexts.extend(ctx)
            all_context.extend(ctx)

        # Genome-level summary
        ctx_pdf = pd.DataFrame(genome_contexts)
        carba_genes = subset['gene_name'].unique().tolist()
        n_composite  = (ctx_pdf['flank_class'] == 'COMPOSITE_TRANSPOSON').sum()
        n_single_up  = (ctx_pdf['flank_class'] == 'SINGLE_IS_UPSTREAM').sum()
        n_single_dn  = (ctx_pdf['flank_class'] == 'SINGLE_IS_DOWNSTREAM').sum()
        n_no_is      = (ctx_pdf['flank_class'] == 'NO_IS').sum()
        is_families  = ','.join(ctx_pdf['is_family'].dropna().unique())

        summary_rows.append({
            'accession':       accession,
            'species':         man_row.iloc[0].get('species', ''),
            'carbapenem_genes':','.join(carba_genes),
            'n_carba_hits':    len(subset),
            'n_is_features':   len(is_df),
            'n_composite':     n_composite,
            'n_single_up':     n_single_up,
            'n_single_dn':     n_single_dn,
            'n_no_is':         n_no_is,
            'is_families_seen':is_families,
        })

        log.info(f'{accession}: {",".join(carba_genes)} | '
                 f'IS features={len(is_df)} | composite={n_composite} | '
                 f'no_IS={n_no_is}')

    # Write outputs
    ctx_df = pd.DataFrame(all_context)
    ctx_df.to_csv(DATA_PROC / 'is_context.tsv', sep='\t', index=False)

    sum_df = pd.DataFrame(summary_rows)
    sum_df.to_csv(DATA_PROC / 'context_summary.tsv', sep='\t', index=False)

    log.info('═' * 60)
    log.info(f'IS context analysis complete:')
    log.info(f'  Context pairs: {len(ctx_df):,}')
    if not sum_df.empty:
        log.info(f'  Composite transposons: {sum_df["n_composite"].sum()}')
        log.info(f'  No flanking IS: {sum_df["n_no_is"].sum()}')
    log.info(f'  Output → {DATA_PROC / "is_context.tsv"}')
    log.info(f'           {DATA_PROC / "context_summary.tsv"}')


if __name__ == '__main__':
    main()
