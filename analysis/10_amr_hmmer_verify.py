"""
Step 10 — AMR gene sequence verification against CARD protein database
─────────────────────────────────────────────────────────────────────────────
Addresses L-001: all 51 AMR hits detected via GFF3 keyword scan (Tier 2);
no sequence-level validation performed. This script provides the missing
sequence-level evidence using pyhmmer phmmer vs CARD protein homolog models.

Method:
  1. Load 51 AMR hits from amr_hits.tsv
  2. For each hit: extract nucleotide sequence from *_genomic.fna.gz
                  (contig[start:stop], reverse-complement if strand=-)
                  translate to protein
  3. Run pyhmmer.phmmer — each AMR protein vs CARD protein database
  4. Record best hit: CARD gene name, ARO, E-value, bitscore, identity
  5. Classify:
       CONFIRMED      — E ≤ 1e-5 and CARD gene name matches our gene name
       NAME_MISMATCH  — E ≤ 1e-5 but gene name differs (different gene!)
       NO_HIT         — no significant hit in CARD (gene absent/truncated)
       EXTRACT_FAIL   — could not extract protein sequence from FNA

Database required:
  data/db/card/protein_fasta_protein_homolog_model.fasta
  Download instructions (printed if file absent):
    1. Visit https://card.mcmaster.ca/download
    2. Download "CARD Data" (broadstreet-vX.X.X.tar.bz2)
    3. Extract: tar -xf broadstreet-vX.X.X.tar.bz2
    4. Copy protein_fasta_protein_homolog_model.fasta to:
       e:/miniconda3/envs/llama-env/amr_project/data/db/card/

Outputs:
  data/processed/amr_hmmer_results.tsv   per-hit HMMER verification results
  logs/amr_hmmer.log

Usage:
  python analysis/10_amr_hmmer_verify.py
"""
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

try:
    import pyhmmer.easel
    import pyhmmer.hmmer
except ImportError:
    print('ERROR: pyhmmer not installed. Run: pip install pyhmmer', file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, LOGS, MANIFEST

# ── constants ─────────────────────────────────────────────────────────────────
HMMER_EVALUE_THR = 1e-5
MIN_PROT_LEN     = 50       # skip fragments < 50 aa
CARD_DIR         = Path(__file__).resolve().parents[1] / 'data' / 'db' / 'card'
CARD_FASTA       = CARD_DIR / 'protein_fasta_protein_homolog_model.fasta'

LOGS.mkdir(parents=True, exist_ok=True)
CARD_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'amr_hmmer.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger('amr_hmmer')

# ── sequence utilities ────────────────────────────────────────────────────────
_CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
_RC = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _domain_identity(dom) -> float:
    """Fraction of identical residues in aligned region."""
    aln = dom.alignment
    total = matched = 0
    for h, t in zip(aln.hmm_sequence, aln.target_sequence):
        if h != '-' and t != '-':
            total += 1
            if h.upper() == t.upper():
                matched += 1
    return matched / total if total else 0.0


def translate(nuc: str) -> str:
    aa = []
    for i in range(0, len(nuc) - 2, 3):
        codon = nuc[i:i + 3].upper()
        aa.append(_CODON.get(codon, 'X'))
    return ''.join(aa).rstrip('*')


def load_contig(fna_gz: Path, target: str) -> str | None:
    """Stream-read one specific contig from compressed FASTA."""
    cid, parts, found = None, [], False
    try:
        with gzip.open(fna_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('>'):
                    if found and parts:
                        return ''.join(parts)
                    cid   = line[1:].split()[0]
                    found = (cid == target)
                    parts = []
                elif found:
                    parts.append(line)
    except Exception as e:
        log.error(f'FNA stream error ({fna_gz.name}): {e}')
        return None
    return ''.join(parts) if found and parts else None


def extract_protein(fna_gz: Path, contig: str, start: int, stop: int,
                    strand: str) -> tuple[str | None, str]:
    """Return (protein_str, status_note)."""
    seq = load_contig(fna_gz, contig)
    if seq is None:
        return None, f'contig {contig} not found in FNA'
    if stop > len(seq):
        return None, f'coordinates {start}:{stop} exceed contig length {len(seq)}'
    nuc = seq[start:stop]
    if strand == '-':
        nuc = revcomp(nuc)
    if len(nuc) < MIN_PROT_LEN * 3:
        return None, f'gene too short ({len(nuc)} bp)'
    prot = translate(nuc)
    if len(prot) < MIN_PROT_LEN:
        return None, f'protein too short after translation ({len(prot)} aa)'
    return prot, 'OK'


# ── CARD header parsing ───────────────────────────────────────────────────────

# CARD protein header formats (multiple versions):
#   >gb|AAM53049.1|KPC-2|ARO:3000589|[Klebsiella...]|BETA-LACTAM|...
#   >gnl|BL_ORD_ID|0|hsp_id:1|...
#   >gb|AAA17477.1|TEM-1|ARO:3000215|...
_ARO_RE  = re.compile(r'ARO:(\d+)')
_GENE_RE = re.compile(
    r'(?:^|[|>])'
    r'((?:KPC|NDM|OXA|IMP|VIM|CTX-M|TEM|SHV|CMY|ACC|DHA|FOX|MOX|ACT|MIR|EBC)'
    r'[-_]?\d+[a-z]?)',
    re.IGNORECASE,
)


def parse_card_header(header: str) -> dict:
    """Extract gene name and ARO from CARD FASTA header."""
    aro = _ARO_RE.search(header)
    gene = _GENE_RE.search(header)
    # Fallback: 3rd pipe-separated field (common in CARD v3)
    parts = header.lstrip('>').split('|')
    pipe_gene = parts[2].strip() if len(parts) >= 3 else ''
    return {
        'card_gene':   (gene.group(1).upper() if gene else pipe_gene),
        'card_aro':    (aro.group(0) if aro else ''),
        'card_header': header[:120],
    }


def gene_names_match(query_gene: str, card_gene: str) -> bool:
    """True if gene names refer to the same resistance gene family."""
    if not query_gene or not card_gene:
        return False
    q = query_gene.upper().strip()
    c = card_gene.upper().strip()
    if q == c:
        return True
    # Prefix match: KPC matches KPC-2; NDM matches NDM-5
    q_base = re.sub(r'[-_]?\d+[A-Z]?$', '', q)
    c_base = re.sub(r'[-_]?\d+[A-Z]?$', '', c)
    return q_base == c_base and len(q_base) >= 2


# ── CARD database check ───────────────────────────────────────────────────────

def check_card_database() -> bool:
    if CARD_FASTA.exists() and CARD_FASTA.stat().st_size > 100_000:
        log.info(f'CARD database: {CARD_FASTA} ({CARD_FASTA.stat().st_size:,} bytes)')
        return True
    log.error(
        f'\nCARD protein database not found at:\n'
        f'  {CARD_FASTA}\n\n'
        f'To download:\n'
        f'  1. Visit https://card.mcmaster.ca/download\n'
        f'  2. Click "CARD Data" → download broadstreet-vX.X.X.tar.bz2\n'
        f'  3. Extract the archive\n'
        f'  4. Copy protein_fasta_protein_homolog_model.fasta to:\n'
        f'     {CARD_DIR}/\n\n'
        f'Alternative (direct download, ~100MB):\n'
        f'  curl -L "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2" '
        f'-o card.tar.bz2 && tar -xf card.tar.bz2 protein_fasta_protein_homolog_model.fasta\n'
        f'  mv protein_fasta_protein_homolog_model.fasta "{CARD_DIR}/"\n'
    )
    return False


# ── CARD protein database loading ─────────────────────────────────────────────

def load_card_proteins(
    alphabet: 'pyhmmer.easel.Alphabet',
) -> tuple['pyhmmer.easel.DigitalSequenceBlock', dict[bytes, dict]]:
    """Load CARD protein FASTA → DigitalSequenceBlock + header metadata dict."""
    log.info(f'Loading CARD proteins from {CARD_FASTA.name} …')
    seqs: list['pyhmmer.easel.DigitalSequence'] = []
    meta: dict[bytes, dict]                      = {}

    with open(CARD_FASTA, 'r', encoding='utf-8', errors='replace') as f:
        header, parts = '', []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header and parts:
                    prot = ''.join(parts)
                    # Strip stop codons and gaps
                    prot = prot.replace('*', '').replace('-', '').replace('X', 'A')
                    if len(prot) >= MIN_PROT_LEN:
                        seq_id = header[1:].split()[0]   # str key — matches hit.name
                        try:
                            ts = pyhmmer.easel.TextSequence(
                                name=seq_id, sequence=prot)
                            ds = ts.digitize(alphabet)
                            seqs.append(ds)
                            meta[seq_id] = parse_card_header(header)
                        except Exception:
                            pass
                header, parts = line, []
            else:
                parts.append(line)
        # last record
        if header and parts:
            prot = ''.join(parts).replace('*', '').replace('-', '').replace('X', 'A')
            if len(prot) >= MIN_PROT_LEN:
                seq_id = header[1:].split()[0]   # str key
                try:
                    ts = pyhmmer.easel.TextSequence(name=seq_id, sequence=prot)
                    seqs.append(ts.digitize(alphabet))
                    meta[seq_id] = parse_card_header(header)
                except Exception:
                    pass

    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
    log.info(f'Loaded {len(seqs):,} CARD proteins')
    return block, meta


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if not check_card_database():
        sys.exit(1)
    if not MANIFEST.exists():
        log.error(f'Manifest not found: {MANIFEST}'); sys.exit(1)
    amr_path = DATA_PROC / 'amr_hits.tsv'
    if not amr_path.exists():
        log.error(f'amr_hits.tsv not found: {amr_path}'); sys.exit(1)

    manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str)
    amr_df   = pd.read_csv(amr_path, sep='\t', dtype=str)
    carba_df = amr_df[amr_df['is_carbapenem'].str.lower() == 'true'].copy()
    log.info(f'AMR hits to verify: {len(carba_df):,} across '
             f'{carba_df["accession"].nunique()} genomes')

    # ── Extract protein sequences ────────────────────────────────────────────
    alphabet = pyhmmer.easel.Alphabet.amino()
    query_seqs:  list = []
    query_meta:  dict[str, dict] = {}   # seq_label → row metadata

    results: list[dict] = []

    for _, row in carba_df.iterrows():
        accession = row['accession']
        gene_name = row['gene_name']
        contig    = str(row['contig'])
        strand    = str(row['strand'])

        try:
            start = int(row['start'])
            stop  = int(row['stop'])
        except (ValueError, TypeError):
            results.append({**row.to_dict(),
                            'protein_len': 0,
                            'card_gene': '', 'card_aro': '', 'card_header': '',
                            'card_evalue': None, 'card_bitscore': None,
                            'card_identity': None,
                            'verification_status': 'EXTRACT_FAIL',
                            'notes': 'invalid coordinates'})
            continue

        man_row = manifest[manifest['accession'] == accession]
        if man_row.empty:
            results.append({**row.to_dict(),
                            'protein_len': 0,
                            'card_gene': '', 'card_aro': '', 'card_header': '',
                            'card_evalue': None, 'card_bitscore': None,
                            'card_identity': None,
                            'verification_status': 'EXTRACT_FAIL',
                            'notes': 'accession not in manifest'})
            continue

        local_dir = Path(man_row.iloc[0]['local_dir'])
        fna_files = list(local_dir.glob('*_genomic.fna.gz'))
        if not fna_files:
            results.append({**row.to_dict(),
                            'protein_len': 0,
                            'card_gene': '', 'card_aro': '', 'card_header': '',
                            'card_evalue': None, 'card_bitscore': None,
                            'card_identity': None,
                            'verification_status': 'EXTRACT_FAIL',
                            'notes': 'FNA file not found'})
            continue

        prot, note = extract_protein(fna_files[0], contig, start, stop, strand)
        if prot is None:
            results.append({**row.to_dict(),
                            'protein_len': 0,
                            'card_gene': '', 'card_aro': '', 'card_header': '',
                            'card_evalue': None, 'card_bitscore': None,
                            'card_identity': None,
                            'verification_status': 'EXTRACT_FAIL',
                            'notes': note})
            log.warning(f'{accession}/{gene_name}: extraction failed — {note}')
            continue

        seq_label = f'{accession}|{gene_name}|{contig}|{start}|{stop}'
        try:
            ts = pyhmmer.easel.TextSequence(name=seq_label, sequence=prot)
            query_seqs.append(ts.digitize(alphabet))
            query_meta[seq_label] = {**row.to_dict(), 'protein_len': len(prot)}
        except Exception as e:
            results.append({**row.to_dict(),
                            'protein_len': len(prot),
                            'card_gene': '', 'card_aro': '', 'card_header': '',
                            'card_evalue': None, 'card_bitscore': None,
                            'card_identity': None,
                            'verification_status': 'EXTRACT_FAIL',
                            'notes': f'digitize error: {e}'})
            continue

        log.debug(f'{accession}/{gene_name}: {len(prot)} aa extracted')

    log.info(f'Proteins extracted: {len(query_seqs):,} / {len(carba_df):,}')

    if not query_seqs:
        log.error('No proteins to search — check FNA files')
        sys.exit(1)

    # ── Load CARD and run phmmer ─────────────────────────────────────────────
    card_block, card_meta = load_card_proteins(alphabet)
    log.info(f'Running phmmer: {len(query_seqs)} AMR proteins vs '
             f'{len(card_block)} CARD sequences …')

    # phmmer: one iteration per query protein
    best_hits: dict[str, dict] = {}   # seq_label → best hit info

    for hits in pyhmmer.hmmer.phmmer(iter(query_seqs), card_block, cpus=0):
        seq_label = hits.query.name      # pyhmmer 0.12: name is str
        best: dict | None = None

        for hit in hits:
            if not hit.included:
                continue
            evalue   = float(hit.evalue)
            bitscore = float(hit.score)
            identity = 0.0
            for dom in hit.domains:
                if dom.included:
                    d = _domain_identity(dom)
                    if d > identity:
                        identity = d

            if best is None or evalue < best['card_evalue']:
                target_meta = card_meta.get(hit.name, {})
                best = {
                    'card_gene':     target_meta.get('card_gene', ''),
                    'card_aro':      target_meta.get('card_aro', ''),
                    'card_header':   target_meta.get('card_header', hit.name),
                    'card_evalue':   evalue,
                    'card_bitscore': bitscore,
                    'card_identity': round(identity, 4),
                }

        best_hits[seq_label] = best or {}

    # ── Compile final results ─────────────────────────────────────────────────
    for seq_label, meta in query_meta.items():
        hit = best_hits.get(seq_label, {})
        if not hit:
            status = 'NO_HIT'
            notes  = 'no significant hit in CARD (E > 1e-5)'
        elif gene_names_match(meta['gene_name'], hit.get('card_gene', '')):
            status = 'CONFIRMED'
            notes  = (f'gene name matches; E={hit["card_evalue"]:.1e}  '
                      f'bitscore={hit["card_bitscore"]:.1f}  '
                      f'identity={hit["card_identity"]:.1%}')
        else:
            status = 'NAME_MISMATCH'
            notes  = (f'our={meta["gene_name"]}  CARD={hit.get("card_gene","")}; '
                      f'E={hit["card_evalue"]:.1e}')

        results.append({
            **meta,
            'card_gene':            hit.get('card_gene', ''),
            'card_aro':             hit.get('card_aro', ''),
            'card_header':          hit.get('card_header', ''),
            'card_evalue':          hit.get('card_evalue'),
            'card_bitscore':        hit.get('card_bitscore'),
            'card_identity':        hit.get('card_identity'),
            'verification_status':  status,
            'notes':                notes,
        })

    # ── Write output ─────────────────────────────────────────────────────────
    out_df = pd.DataFrame(results)
    out_path = DATA_PROC / 'amr_hmmer_results.tsv'
    out_df.to_csv(out_path, sep='\t', index=False)

    # ── Summary ──────────────────────────────────────────────────────────────
    status_counts = out_df['verification_status'].value_counts().to_dict()
    n_confirmed     = status_counts.get('CONFIRMED', 0)
    n_mismatch      = status_counts.get('NAME_MISMATCH', 0)
    n_no_hit        = status_counts.get('NO_HIT', 0)
    n_fail          = status_counts.get('EXTRACT_FAIL', 0)
    n_total         = len(out_df)

    log.info('═' * 65)
    log.info(f'AMR HMMER VERIFICATION COMPLETE')
    log.info(f'  Total AMR hits       : {n_total}')
    log.info(f'  CONFIRMED            : {n_confirmed} '
             f'({100*n_confirmed/n_total:.1f}%)  ← sequence-verified')
    log.info(f'  NAME_MISMATCH        : {n_mismatch} '
             f'(hit found but gene name differs — investigate)')
    log.info(f'  NO_HIT               : {n_no_hit} '
             f'(no CARD hit; gene absent, truncated, or novel)')
    log.info(f'  EXTRACT_FAIL         : {n_fail} '
             f'(FNA/coord issue — check logs/amr_hmmer.log)')

    if n_mismatch > 0:
        mm = out_df[out_df['verification_status'] == 'NAME_MISMATCH']
        log.info('  NAME_MISMATCH detail:')
        for _, r in mm.iterrows():
            log.info(f'    {r["accession"]}  {r["gene_name"]} → CARD: {r["card_gene"]} '
                     f'(E={r["card_evalue"]:.1e})')

    tier_upgrade = (
        'HMMER_VERIFIED' if n_confirmed == n_total - n_fail
        else f'PARTIAL ({n_confirmed}/{n_total - n_fail} confirmed)'
    )
    log.info(f'  Detection tier upgrade: GFF_KEYWORD → {tier_upgrade}')
    log.info(f'  → {out_path}')


if __name__ == '__main__':
    main()
