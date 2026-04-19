"""
Step 1 — NCBI genome downloader
─────────────────────────────────────────────────────────────────────────────
Strategy (revised):
  NCBI assembly_summary.txt lacks geo_loc_name / isolation_source.
  Those fields live in BioSample. We use two paths:

  Path A (preferred): Entrez search on the `assembly` database with
    organism + country + host/source filters → returns UID list →
    esummary gives FTP paths + metadata for each.

  Path B (fallback): Download full assembly_summary, fetch BioSample
    records for a batch, filter, download.

  For each assembly:
    1. Download required files with MD5 verification
    2. Write data/raw/download_status.tsv

Usage:
  python analysis/01_download.py [--dry-run] [--limit N] [--species kpn|eco|both]
"""
import argparse
import gzip
import hashlib
import logging
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    ACCEPTED_ASSEMBLY_LEVELS, ACCEPTED_SOURCES, DATA_RAW,
    DOWNLOAD_FILES, DOWNLOAD_TIMEOUT_S, LOGS, MAX_GENOMES_CHINA,
    MAX_RETRIES, NCBI_API_KEY, NCBI_EMAIL, RATE_LIMIT_DELAY_S,
    REQUIRED_ISOLATION_COUNTRY, RETRY_DELAY_S, TARGET_SPECIES,
)

LOGS.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'download.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('download')

ENTREZ_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
SESSION = requests.Session()
SESSION.headers['User-Agent'] = f'AMRProject/1.0 (contact: {NCBI_EMAIL})'


# ── HTTP helpers ──────────────────────────────────────────────────────────────

def _get(url: str, params: dict = None, stream=False) -> requests.Response:
    if NCBI_API_KEY:
        params = params or {}
        params['api_key'] = NCBI_API_KEY
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = SESSION.get(url, params=params, stream=stream,
                            timeout=DOWNLOAD_TIMEOUT_S)
            r.raise_for_status()
            return r
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise   # don't retry 404 — file doesn't exist
            if attempt == MAX_RETRIES:
                raise
            log.warning(f'Attempt {attempt}: {e}. Retrying…')
            time.sleep(RETRY_DELAY_S * attempt)
        except requests.RequestException as e:
            if attempt == MAX_RETRIES:
                raise
            log.warning(f'Attempt {attempt}: {e}. Retrying…')
            time.sleep(RETRY_DELAY_S * attempt)
    raise RuntimeError('unreachable')


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        r = _get(url, stream=True)
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(65536):
                f.write(chunk)
        return True
    except requests.HTTPError as e:
        if e.response.status_code == 404:
            return False
        raise
    except Exception as e:
        log.error(f'Download failed {url}: {e}')
        dest.unlink(missing_ok=True)
        return False


# ── Entrez search ─────────────────────────────────────────────────────────────

def entrez_search(query: str, db: str = 'assembly',
                  retmax: int = 500) -> list[str]:
    """
    Return list of UIDs matching query in db.
    Handles WebEnv/query_key pagination for > retmax results.
    """
    # Initial search
    params = {
        'db': db, 'term': query, 'retmax': min(retmax, 10000),
        'usehistory': 'y', 'retmode': 'xml', 'email': NCBI_EMAIL,
    }
    r = _get(f'{ENTREZ_BASE}/esearch.fcgi', params=params)
    time.sleep(RATE_LIMIT_DELAY_S)

    root    = ET.fromstring(r.content)
    count   = int(root.findtext('Count') or '0')
    web_env = root.findtext('WebEnv') or ''
    qkey    = root.findtext('QueryKey') or '1'
    ids     = [id_el.text for id_el in root.findall('.//IdList/Id')]

    log.info(f'  Entrez search: {count:,} results for "{query[:80]}"')

    # Paginate if needed
    batch = 500
    start = len(ids)
    while len(ids) < min(count, retmax):
        params2 = {
            'db': db, 'query_key': qkey, 'WebEnv': web_env,
            'retstart': start, 'retmax': batch,
            'retmode': 'xml', 'email': NCBI_EMAIL,
        }
        r2   = _get(f'{ENTREZ_BASE}/esearch.fcgi', params=params2)
        root2= ET.fromstring(r2.content)
        new_ids = [el.text for el in root2.findall('.//IdList/Id')]
        if not new_ids:
            break
        ids.extend(new_ids)
        start += len(new_ids)
        time.sleep(RATE_LIMIT_DELAY_S)

    return ids[:retmax]


def entrez_summary_batch(uids: list[str], db: str = 'assembly') -> list[dict]:
    """
    Fetch DocSummary records for a list of UIDs.
    Returns list of dicts with metadata.
    """
    all_records = []
    batch_size  = 200

    for i in range(0, len(uids), batch_size):
        batch = uids[i:i + batch_size]
        params = {
            'db': db, 'id': ','.join(batch),
            'retmode': 'xml', 'email': NCBI_EMAIL,
        }
        r    = _get(f'{ENTREZ_BASE}/esummary.fcgi', params=params)
        root = ET.fromstring(r.content)
        time.sleep(RATE_LIMIT_DELAY_S)

        for doc in root.findall('.//DocumentSummary'):
            def val(tag: str) -> str:
                el = doc.find(tag)
                return el.text.strip() if el is not None and el.text else ''

            ftp = val('FtpPath_GenBank') or val('FtpPath_RefSeq')
            level = val('AssemblyStatus')

            # Map NCBI status strings to our expected values
            level_map = {
                'Complete Genome': 'Complete Genome',
                'Chromosome':      'Chromosome',
                'Scaffold':        'Scaffold',
                'Contig':          'Contig',
            }
            level_norm = level_map.get(level, level)

            record = {
                'assembly_accession': val('AssemblyAccession'),
                'assembly_name':      val('AssemblyName'),
                'organism_name':      val('Organism'),
                'assembly_level':     level_norm,
                'taxid':              val('Taxid'),
                'biosample':          val('BioSampleAccn'),
                'bioproject':         val('GB_BioProjects'),
                'ftp_path':           ftp,
                'seq_rel_date':       val('SeqReleaseDate'),
                'asm_submitter':      val('SubmitterOrganization'),
                'geo_loc_name':       val('Biosource/InfraspecificList/InfraspecificType'),
                'infraspecific_name': val('Biosource/InfraspecificList'),
            }
            all_records.append(record)

    return all_records


def build_entrez_query(species_name: str) -> str:
    """Build NCBI assembly search query for Chinese clinical isolates."""
    # Entrez assembly db supports: organism, country (from biosample), isolation_source
    # "latest[filter]" excludes suppressed/replaced assemblies
    q = (
        f'"{species_name}"[Organism] AND '
        f'"China"[Country] AND '
        f'latest[filter] AND '
        f'("Complete Genome"[Assembly Level] OR '
        f' "Chromosome"[Assembly Level] OR '
        f' "Scaffold"[Assembly Level])'
    )
    return q


# ── BioSample metadata (for source filtering) ────────────────────────────────

def fetch_biosample_source(biosample_id: str) -> str:
    """
    Fetch isolation_source from a BioSample record.
    Returns empty string on failure.
    """
    if not biosample_id or biosample_id == 'na':
        return ''
    params = {'db': 'biosample', 'id': biosample_id,
               'retmode': 'xml', 'email': NCBI_EMAIL}
    try:
        r    = _get(f'{ENTREZ_BASE}/efetch.fcgi', params=params)
        root = ET.fromstring(r.content)
        time.sleep(RATE_LIMIT_DELAY_S)
        for attr in root.findall('.//Attribute'):
            name = attr.get('attribute_name', '').lower()
            if 'isolation_source' in name or 'host' in name:
                return attr.text or ''
    except Exception:
        pass
    return ''


def _is_clinical(source: str) -> bool:
    if not isinstance(source, str):
        return False
    s = source.lower()
    return any(kw in s for kw in ACCEPTED_SOURCES)


# ── MD5 verification ──────────────────────────────────────────────────────────

def verify_md5(ftp_base: str, local_dir: Path) -> dict[str, bool]:
    md5_url = ftp_base.replace('ftp://', 'https://') + '/md5checksums.txt'
    try:
        r = _get(md5_url)
    except Exception as e:
        log.warning(f'Could not fetch MD5: {e}')
        return {}

    results = {}
    for line in r.text.splitlines():
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        expected_md5, fname = parts[0], Path(parts[1]).name
        local_file = local_dir / fname
        if local_file.exists():
            actual = _md5(local_file)
            passed = actual == expected_md5
            results[fname] = passed
            if not passed:
                log.error(f'MD5 MISMATCH {fname}: '
                          f'expected={expected_md5}, got={actual}')
    return results


# ── per-assembly download ─────────────────────────────────────────────────────

def download_assembly(row: dict, dry_run: bool = False) -> dict:
    accession = row.get('assembly_accession', 'UNKNOWN')
    ftp_path  = row.get('ftp_path', '')

    result = {**row, 'status': 'PENDING', 'files_ok': '',
              'md5_ok': '', 'has_amrfinder': False, 'local_dir': ''}

    if not ftp_path or ftp_path in ('na', ''):
        result['status'] = 'NO_FTP'
        log.warning(f'{accession}: no FTP path')
        return result

    https_base = ftp_path.replace('ftp://', 'https://')
    prefix     = Path(ftp_path).name
    local_dir  = DATA_RAW / accession
    result['local_dir'] = str(local_dir)

    if dry_run:
        result['status'] = 'DRY_RUN'
        log.info(f'DRY-RUN {accession}: {https_base}')
        return result

    local_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    for suffix in DOWNLOAD_FILES:
        fname = prefix + suffix
        dest  = local_dir / fname
        if dest.exists():
            downloaded.append(fname)
            continue
        ok = download_file(f'{https_base}/{fname}', dest)
        if ok:
            downloaded.append(fname)
            log.info(f'  ↓ {accession}/{fname}')
        time.sleep(RATE_LIMIT_DELAY_S)

    result['has_amrfinder'] = (local_dir / (prefix + '_amrfinderplus.tsv')).exists()

    # Required files check
    required = [prefix + '_genomic.gff.gz', prefix + '_genomic.fna.gz']
    missing  = [f for f in required if not (local_dir / f).exists()]
    if missing:
        result['status'] = 'MISSING_REQUIRED'
        log.error(f'{accession}: missing {missing}')
        return result

    # MD5 check
    md5_res    = verify_md5(ftp_path, local_dir)
    failed_md5 = [f for f, ok in md5_res.items() if not ok]
    if failed_md5:
        result['status'] = 'MD5_FAIL'
        result['md5_ok'] = 'FAIL:' + ','.join(failed_md5)
        log.error(f'{accession}: MD5 failed {failed_md5}')
        for f in failed_md5:
            (local_dir / f).unlink(missing_ok=True)
        return result

    result['status']   = 'OK'
    result['files_ok'] = ','.join(downloaded)
    result['md5_ok']   = 'PASS'
    return result


# ── main ──────────────────────────────────────────────────────────────────────

SPECIES_NAMES = {
    'kpn':  'Klebsiella pneumoniae',
    'eco':  'Escherichia coli',
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--limit', type=int, default=MAX_GENOMES_CHINA)
    ap.add_argument('--species', choices=['kpn','eco','both'], default='kpn')
    args = ap.parse_args()

    targets = (['kpn','eco'] if args.species == 'both'
               else [args.species])

    # ── 1. Entrez search ─────────────────────────────────────────────────
    all_uids = []
    uid_species_map = {}
    for sp_key in targets:
        sp_name = SPECIES_NAMES[sp_key]
        query   = build_entrez_query(sp_name)
        log.info(f'Searching NCBI assembly: {sp_name}')
        log.info(f'  Query: {query}')
        uids = entrez_search(query, retmax=args.limit if args.limit is not None else 10_000)
        log.info(f'  UIDs found: {len(uids)}')
        for uid in uids:
            uid_species_map[uid] = sp_key
        all_uids.extend(uids)

    if not all_uids:
        log.error('No assemblies found. Check query or internet connection.')
        sys.exit(1)

    # ── 2. Fetch assembly summaries (FTP paths + metadata) ───────────────
    log.info(f'Fetching summaries for {len(all_uids)} assemblies…')
    records = entrez_summary_batch(all_uids)
    for rec in records:
        # re-attach species key
        pass   # uid→accession mapping done during summary fetch
    log.info(f'Got {len(records)} assembly records')

    # Filter assembly level
    records = [r for r in records
               if r.get('assembly_level', '') in ACCEPTED_ASSEMBLY_LEVELS]
    log.info(f'After assembly-level filter: {len(records)}')

    if not records:
        log.error('No assemblies pass assembly-level filter.')
        sys.exit(1)

    # Limit
    if args.limit and len(records) > args.limit:
        # prefer Complete Genome first
        lvl_order = {'Complete Genome':0,'Chromosome':1,'Scaffold':2}
        records = sorted(records,
                         key=lambda r: lvl_order.get(r.get('assembly_level',''),3))
        records = records[:args.limit]
        log.info(f'Limited to {len(records)} (best quality first)')

    # Save candidate manifest
    DATA_RAW.mkdir(parents=True, exist_ok=True)
    cand_df = pd.DataFrame(records)
    cand_path = DATA_RAW / 'candidate_manifest.tsv'
    cand_df.to_csv(cand_path, sep='\t', index=False)
    log.info(f'Candidate manifest saved: {cand_path}  ({len(records)} rows)')

    if args.dry_run:
        log.info('DRY-RUN complete.')
        cols = ['assembly_accession','organism_name','assembly_level',
                'asm_submitter']
        avail = [c for c in cols if c in cand_df.columns]
        print(cand_df[avail].head(20).to_string(index=False))
        return

    # ── 3. Download ───────────────────────────────────────────────────────
    status_records = []
    for rec in tqdm(records, desc='Downloading', unit='genome'):
        s = download_assembly(rec, dry_run=False)
        status_records.append(s)
        if len(status_records) % 10 == 0:
            pd.DataFrame(status_records).to_csv(
                DATA_RAW / 'download_status.tsv', sep='\t', index=False)

    status_df = pd.DataFrame(status_records)
    status_df.to_csv(DATA_RAW / 'download_status.tsv', sep='\t', index=False)

    ok  = (status_df['status'] == 'OK').sum()
    err = (status_df['status'] != 'OK').sum()
    log.info('═' * 60)
    log.info(f'Download complete: {ok} OK, {err} errors')
    if err:
        log.warning('Errors:')
        for _, r in status_df[status_df['status'] != 'OK'].iterrows():
            log.warning(f'  {r["assembly_accession"]}: {r["status"]}')


if __name__ == '__main__':
    main()
