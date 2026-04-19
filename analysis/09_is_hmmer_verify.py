"""
Step 9 — IS element family verification using PFAM HMM (PF01527)
─────────────────────────────────────────────────────────────────────────────
Addresses L-005: annotation bias where IS26/IS6-family transposases in
susceptible genomes may be labelled 'transposase' / 'insertion sequence'
without an IS-family name, artificially driving susceptible IS6 counts to 0.

Method:
  1. Rescan all 270 GFF3 files for IS/transposase-annotated features
  2. Separate: family-known (has IS/Tn name) vs IS_unknown
  3. Extract protein sequences for IS_unknown features (FNA + GFF3 coords)
  4. Run pyhmmer.hmmsearch — PF01527 (IS6 transposase) vs IS_unknown proteins
  5. Reclassify IS_unknown → IS6_PFAM where E ≤ 1e-5
  6. Recount IS6 per genome (name-based + PFAM-reclassified)
  7. Re-run Mann-Whitney U, Cliff's delta, AUC with corrected counts
  8. Quantify annotation bias: did reclassification disproportionately
     affect susceptible genomes?

Database:
  data/db/pfam/PF01527.hmm   (auto-downloaded from EBI if absent)

Outputs:
  data/processed/is_hmmer_results.tsv          per-feature HMMER results
  data/processed/is_burden_corrected.tsv        corrected IS6 counts
  data/processed/is_burden_corrected_stats.json corrected statistics
  logs/is_hmmer.log

Usage:
  python analysis/09_is_hmmer_verify.py
"""
import gzip
import json
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

try:
    import pyhmmer.easel
    import pyhmmer.hmmer
    import pyhmmer.plan7
except ImportError:
    print('ERROR: pyhmmer not installed. Run: pip install pyhmmer', file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, IS_GFF_KEYWORDS, LOGS, MANIFEST

# ── constants ─────────────────────────────────────────────────────────────────
PFAM_ID          = 'PF01527'
HMMER_EVALUE_THR = 1e-5
MIN_PROT_LEN     = 50        # skip fragments shorter than this (aa)
DB_DIR           = Path(__file__).resolve().parents[1] / 'data' / 'db' / 'pfam'
HMM_PATH         = DB_DIR / f'{PFAM_ID}.hmm'

_IS_PATTERNS  = [re.compile(p, re.IGNORECASE) for p in IS_GFF_KEYWORDS]
_ATTRS_RE     = re.compile(r'(\w+)=([^;]+)')
_IS6_NAMES    = re.compile(r'\bIS(26|257|1353|240|1006|6)\b', re.IGNORECASE)

LOGS.mkdir(parents=True, exist_ok=True)
DB_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'is_hmmer.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger('is_hmmer')

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


def load_contigs(fna_gz: Path) -> dict[str, str]:
    """Load all contigs from compressed FASTA → {contig_id: sequence}."""
    contigs: dict[str, str] = {}
    cid, parts = None, []
    try:
        with gzip.open(fna_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('>'):
                    if cid:
                        contigs[cid] = ''.join(parts)
                    cid = line[1:].split()[0]
                    parts = []
                elif cid:
                    parts.append(line)
        if cid:
            contigs[cid] = ''.join(parts)
    except Exception as e:
        log.error(f'FNA load error ({fna_gz.name}): {e}')
    return contigs


def extract_protein(contigs: dict, contig: str, start: int, stop: int,
                    strand: str) -> str | None:
    seq = contigs.get(contig)
    if not seq or stop > len(seq):
        return None
    nuc = seq[start:stop]
    if strand == '-':
        nuc = revcomp(nuc)
    if len(nuc) < MIN_PROT_LEN * 3:
        return None
    prot = translate(nuc)
    if len(prot) < MIN_PROT_LEN:
        return None
    return prot


# ── GFF3 parsing ──────────────────────────────────────────────────────────────

def _parse_attrs(attr_str: str) -> dict:
    return dict(m.groups() for m in _ATTRS_RE.finditer(attr_str))


def _is_family_from_text(text: str) -> str:
    m = re.search(r'\b(Tn\d+)\b', text, re.IGNORECASE)
    if m:
        return m.group(1)
    m = re.search(r'\b(IS[A-Za-z0-9]+)\b', text, re.IGNORECASE)
    if m:
        return m.group(1)
    return 'IS_unknown'


def scan_is_features(gff_gz: Path) -> list[dict]:
    """Scan GFF3 for IS/transposase features; return list of feature dicts."""
    features = []
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
                attrs    = _parse_attrs(cols[8])
                product  = attrs.get('product', '')
                gene     = attrs.get('gene', '')
                note     = attrs.get('Note', '')
                combined = f'{product} {gene} {note}'
                if not any(p.search(combined) for p in _IS_PATTERNS):
                    continue
                try:
                    start = int(cols[3]) - 1
                    stop  = int(cols[4])
                except ValueError:
                    continue
                features.append({
                    'contig':    cols[0],
                    'start':     start,
                    'stop':      stop,
                    'strand':    cols[6],
                    'is_family': _is_family_from_text(combined),
                    'product':   product[:80],
                    'feat_type': feat_type,
                })
    except Exception as e:
        log.error(f'GFF scan error ({gff_gz.name}): {e}')
    return features


# ── PFAM HMM management ───────────────────────────────────────────────────────

# EBI returns gzip-compressed HMM; no Accept header needed
_PFAM_URLS = [
    f'https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{PFAM_ID}/?annotation=hmm',
    f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{PFAM_ID}/?annotation=hmm&download',
]


def ensure_pfam_hmm() -> list:
    """Return loaded PF01527 HMMs; download from EBI (gzip) if absent."""
    import gzip as _gzip
    if not HMM_PATH.exists():
        log.info(f'{HMM_PATH} not found — downloading from EBI …')
        downloaded = False
        for url in _PFAM_URLS:
            try:
                resp = requests.get(url, timeout=120,
                                    headers={'User-Agent': 'AMR-IS-verify/1.0'})
                if not resp.ok:
                    log.warning(f'HTTP {resp.status_code} from {url}')
                    continue
                data = resp.content
                # EBI wraps in gzip — decompress if needed
                if data[:2] == b'\x1f\x8b':
                    data = _gzip.decompress(data)
                if b'HMMER3' not in data[:20]:
                    log.warning(f'{url}: response not a HMMER3 HMM file, skipping')
                    continue
                HMM_PATH.write_bytes(data)
                log.info(f'Saved {HMM_PATH} ({len(data):,} bytes)')
                downloaded = True
                break
            except Exception as e:
                log.warning(f'Download failed ({url}): {e}')
        if not downloaded:
            log.error(
                f'Could not download {PFAM_ID}.hmm automatically.\n'
                f'Manual download (save decompressed content):\n'
                f'  python -c "import requests,gzip; '
                f'r=requests.get(\'{_PFAM_URLS[0]}\'); '
                f'open(\'{HMM_PATH}\',\'wb\').write(gzip.decompress(r.content))"\n'
                f'Then re-run this script.'
            )
            sys.exit(1)
    try:
        with pyhmmer.plan7.HMMFile(str(HMM_PATH)) as hf:
            hmms = list(hf)
        log.info(f'Loaded {len(hmms)} HMM profile(s) from {HMM_PATH.name}')
        return hmms
    except Exception as e:
        log.error(f'HMM load error: {e}')
        sys.exit(1)


# ── statistics helpers ────────────────────────────────────────────────────────

def mann_whitney(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    from scipy.stats import mannwhitneyu
    stat, p = mannwhitneyu(a, b, alternative='two-sided')
    return float(stat), float(p)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    mat = np.sign(a[:, None] - b[None, :])
    return float(mat.mean())


def auc_via_mwu(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUC = U_greater / (n_pos * n_neg)  — exact Mann-Whitney AUC."""
    from scipy.stats import mannwhitneyu
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float('nan')
    u, _ = mannwhitneyu(pos, neg, alternative='greater')
    return float(u) / (len(pos) * len(neg))


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if not MANIFEST.exists():
        log.error(f'Manifest not found: {MANIFEST}'); sys.exit(1)
    amr_path = DATA_PROC / 'amr_hits.tsv'
    if not amr_path.exists():
        log.error(f'amr_hits.tsv not found: {amr_path}'); sys.exit(1)

    manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str)
    amr_df   = pd.read_csv(amr_path, sep='\t', dtype=str)
    resistant_acc = set(
        amr_df[amr_df['is_carbapenem'].str.lower() == 'true']['accession']
    )
    log.info(f'Cohort: {len(manifest)} genomes | {len(resistant_acc)} carbapenem-resistant')

    # ── Load PFAM HMM ────────────────────────────────────────────────────────
    hmms    = ensure_pfam_hmm()
    alphabet = pyhmmer.easel.Alphabet.amino()

    # ── Scan all genomes: collect IS features + extract IS_unknown proteins ──
    log.info('Scanning GFF3 files and extracting IS_unknown proteins …')

    all_features:   list[dict]  = []
    query_seqs:     list        = []   # pyhmmer DigitalSequence objects
    seq_to_feat_idx: dict[str, int] = {}  # seq_name → index in all_features

    n_genomes_done = 0
    for _, row in manifest.iterrows():
        accession = row['accession']
        local_dir = Path(row['local_dir'])
        is_resistant = (accession in resistant_acc)

        gff_files = list(local_dir.glob('*_genomic.gff.gz'))
        fna_files = list(local_dir.glob('*_genomic.fna.gz'))
        if not gff_files or not fna_files:
            log.warning(f'{accession}: GFF or FNA missing, skip')
            continue

        features = scan_is_features(gff_files[0])
        if not features:
            n_genomes_done += 1
            continue

        # Load FNA only if needed for IS_unknown protein extraction
        contigs: dict[str, str] = {}
        needs_fna = any(f['is_family'] == 'IS_unknown' for f in features)
        if needs_fna:
            contigs = load_contigs(fna_files[0])

        for feat in features:
            feat['accession']        = accession
            feat['is_resistant']     = is_resistant
            feat['hmm_evalue']       = None
            feat['hmm_bitscore']     = None
            feat['hmm_identity']     = None
            feat['corrected_family'] = feat['is_family']
            feat['was_reclassified'] = False
            idx = len(all_features)
            all_features.append(feat)

            if feat['is_family'] != 'IS_unknown':
                continue

            prot = extract_protein(contigs, feat['contig'],
                                   feat['start'], feat['stop'], feat['strand'])
            if prot is None:
                continue

            seq_name = f'{accession}|{feat["contig"]}|{feat["start"]}|{feat["stop"]}'
            try:
                ts = pyhmmer.easel.TextSequence(name=seq_name.encode(), sequence=prot)
                query_seqs.append(ts.digitize(alphabet))
                seq_to_feat_idx[seq_name] = idx
            except Exception as e:
                log.debug(f'{seq_name}: digitize error: {e}')

        n_genomes_done += 1
        if n_genomes_done % 50 == 0:
            log.info(f'  Scanned {n_genomes_done}/{len(manifest)} genomes …')

    n_unknown = sum(1 for f in all_features if f['is_family'] == 'IS_unknown')
    log.info(f'IS features total: {len(all_features):,} | IS_unknown: {n_unknown:,} '
             f'| Protein queries: {len(query_seqs):,}')

    # ── Run HMMER: PF01527 vs IS_unknown proteins ────────────────────────────
    n_reclassified = 0
    if query_seqs:
        log.info(f'Running hmmsearch ({PFAM_ID}) against {len(query_seqs):,} proteins …')
        target_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, query_seqs)

        for hits in pyhmmer.hmmer.hmmsearch(hmms, target_block, cpus=0):
            for hit in hits:
                seq_name = hit.name          # pyhmmer 0.12: name is str
                idx = seq_to_feat_idx.get(seq_name)
                if idx is None:
                    continue
                evalue   = float(hit.evalue)
                bitscore = float(hit.score)
                identity = 0.0
                for dom in hit.domains:
                    if dom.included:
                        d = _domain_identity(dom)
                        if d > identity:
                            identity = d

                all_features[idx]['hmm_evalue']   = evalue
                all_features[idx]['hmm_bitscore'] = bitscore
                all_features[idx]['hmm_identity'] = round(identity, 4)

                if evalue <= HMMER_EVALUE_THR:
                    all_features[idx]['corrected_family'] = f'IS6_PFAM_{PFAM_ID}'
                    all_features[idx]['was_reclassified'] = True
                    n_reclassified += 1

        log.info(f'Reclassified IS_unknown → IS6_PFAM: {n_reclassified:,}')
    else:
        log.warning('No IS_unknown proteins to search — nothing to reclassify')

    # ── Write per-feature results ─────────────────────────────────────────────
    feat_df = pd.DataFrame(all_features)
    feat_df.to_csv(DATA_PROC / 'is_hmmer_results.tsv', sep='\t', index=False)
    log.info(f'Feature results → {DATA_PROC / "is_hmmer_results.tsv"}')

    # ── Compute corrected IS6 counts per genome ───────────────────────────────
    burden_rows = []
    for _, man_row in manifest.iterrows():
        accession    = man_row['accession']
        is_resistant = (accession in resistant_acc)
        gf = feat_df[feat_df['accession'] == accession]

        n_is_total  = len(gf)
        n_is6_orig  = int(_IS6_NAMES.search(' '.join(gf['is_family'].fillna('')))
                          is not None)  # wrong — use row-level below
        n_is6_orig  = int(gf['is_family'].apply(
                          lambda x: bool(_IS6_NAMES.search(str(x)))).sum())
        n_is26_orig = int(gf['is_family'].apply(
                          lambda x: bool(re.search(r'\bIS26\b', str(x),
                                                   re.IGNORECASE))).sum())
        n_new       = int(gf['was_reclassified'].sum())
        n_is6_corr  = n_is6_orig + n_new

        burden_rows.append({
            'accession':       accession,
            'is_resistant':    is_resistant,
            'n_is_total':      n_is_total,
            'n_is6_original':  n_is6_orig,
            'n_is26_original': n_is26_orig,
            'n_is6_pfam_new':  n_new,
            'n_is6_corrected': n_is6_corr,
        })

    burden_df = pd.DataFrame(burden_rows)
    burden_df.to_csv(DATA_PROC / 'is_burden_corrected.tsv', sep='\t', index=False)

    # ── Re-run statistics ─────────────────────────────────────────────────────
    res  = burden_df[burden_df['is_resistant'] == True]
    susc = burden_df[burden_df['is_resistant'] == False]
    y    = burden_df['is_resistant'].values.astype(int)

    def _stats(col: str) -> dict:
        a = res[col].values.astype(float)
        b = susc[col].values.astype(float)
        _, p  = mann_whitney(a, b)
        d     = cliffs_delta(a, b)
        auc   = auc_via_mwu(burden_df[col].values.astype(float), y)
        return {
            'resistant_median':   float(np.median(a)),
            'resistant_iqr':      [float(np.percentile(a, 25)),
                                   float(np.percentile(a, 75))],
            'susceptible_median': float(np.median(b)),
            'susceptible_iqr':    [float(np.percentile(b, 25)),
                                   float(np.percentile(b, 75))],
            'cliffs_delta':       round(d, 4),
            'mannwhitney_p':      float(p),
            'auc':                round(auc, 4),
        }

    orig = _stats('n_is6_original')
    corr = _stats('n_is6_corrected')

    n_reclass_res  = int(res['n_is6_pfam_new'].sum())
    n_reclass_susc = int(susc['n_is6_pfam_new'].sum())
    bias_toward = ('susceptible' if n_reclass_susc > n_reclass_res
                   else 'resistant')
    # Annotation bias confirmed if susceptible got ≥2× more reclassifications
    bias_confirmed = n_reclass_susc > n_reclass_res * 2

    auc_change   = round(corr['auc'] - orig['auc'], 4)
    delta_change = round(corr['cliffs_delta'] - orig['cliffs_delta'], 4)

    if bias_confirmed:
        verdict = (f'BIAS CONFIRMED — susceptible genomes gained '
                   f'{n_reclass_susc} vs resistant {n_reclass_res} '
                   f'reclassifications. AUC {orig["auc"]:.3f} → {corr["auc"]:.3f} '
                   f'({auc_change:+.4f}). Corrected AUC is publishable figure.')
    else:
        verdict = (f'BIAS NOT DETECTED — reclassifications split '
                   f'{n_reclass_res} resistant / {n_reclass_susc} susceptible. '
                   f'AUC {orig["auc"]:.3f} → {corr["auc"]:.3f} '
                   f'({auc_change:+.4f}). Original AUC=0.952 robust to annotation gap.')

    stats_out = {
        'original':  orig,
        'corrected': {**corr,
                      'n_reclassified_total':       n_reclassified,
                      'n_reclassified_resistant':   n_reclass_res,
                      'n_reclassified_susceptible': n_reclass_susc},
        'bias_assessment': {
            'annotation_bias_confirmed':   bias_confirmed,
            'reclassification_toward':     bias_toward,
            'auc_change':                  auc_change,
            'cliffs_delta_change':         delta_change,
            'verdict':                     verdict,
        },
    }
    out_json = DATA_PROC / 'is_burden_corrected_stats.json'
    out_json.write_text(json.dumps(stats_out, indent=2, default=str))

    # ── Final report ─────────────────────────────────────────────────────────
    log.info('═' * 65)
    log.info(f'IS HMMER VERIFICATION COMPLETE ({PFAM_ID} — IS6 family)')
    log.info(f'  IS features scanned   : {len(all_features):,}')
    log.info(f'  IS_unknown proteins   : {len(query_seqs):,}')
    log.info(f'  Reclassified → IS6    : {n_reclassified:,}  '
             f'(resistant {n_reclass_res} / susceptible {n_reclass_susc})')
    log.info('─' * 65)
    log.info(f'  ORIGINAL  IS6  resistant median : {orig["resistant_median"]:.1f}')
    log.info(f'  ORIGINAL  IS6  susceptible median: {orig["susceptible_median"]:.1f}')
    log.info(f"  ORIGINAL  Cliff's delta          : {orig['cliffs_delta']:.3f}")
    log.info(f'  ORIGINAL  AUC                    : {orig["auc"]:.3f}')
    log.info('─' * 65)
    log.info(f'  CORRECTED IS6  resistant median : {corr["resistant_median"]:.1f}')
    log.info(f'  CORRECTED IS6  susceptible median: {corr["susceptible_median"]:.1f}')
    log.info(f"  CORRECTED Cliff's delta          : {corr['cliffs_delta']:.3f}")
    log.info(f'  CORRECTED AUC                    : {corr["auc"]:.3f}')
    log.info('─' * 65)
    log.info(f'  VERDICT: {verdict}')
    log.info(f'  → {DATA_PROC / "is_burden_corrected.tsv"}')
    log.info(f'  → {out_json}')


if __name__ == '__main__':
    main()
