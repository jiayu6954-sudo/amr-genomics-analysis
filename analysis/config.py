"""
Single source of truth for all parameters.
Change thresholds here only — never hard-code them in pipeline scripts.
"""
from pathlib import Path

# ── project paths ─────────────────────────────────────────────────────────────
ROOT         = Path(__file__).resolve().parents[1]
DATA_RAW     = ROOT / 'data' / 'raw'
DATA_PROC    = ROOT / 'data' / 'processed'
DATA_VAL     = ROOT / 'data' / 'validated'
LOGS         = ROOT / 'logs'
FIGURES      = ROOT / 'figures'
REPORTS      = ROOT / 'reports'

# manifest: validated genome registry — single source of truth
MANIFEST     = DATA_VAL / 'genome_manifest.tsv'

# ── NCBI credentials (required for Entrez API) ────────────────────────────────
NCBI_EMAIL   = 'jiayu6954@gmail.com'   # used in HTTP User-Agent headers
NCBI_API_KEY = ''                       # optional: get at ncbi.nlm.nih.gov/account

# NCBI assembly summary URL templates
ASSEMBLY_SUMMARY_URLS = {
    'Klebsiella_pneumoniae':
        'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/'
        'Klebsiella_pneumoniae/assembly_summary.txt',
    'Escherichia_coli':
        'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/'
        'Escherichia_coli/assembly_summary.txt',
}

# ── target species ────────────────────────────────────────────────────────────
TARGET_SPECIES = {
    'Klebsiella pneumoniae': {
        'genome_size_min_bp': 4_800_000,
        'genome_size_max_bp': 6_500_000,
        'min_cds': 3_000,
        'min_n50': 50_000,
    },
    'Escherichia coli': {
        'genome_size_min_bp': 4_200_000,
        'genome_size_max_bp': 5_800_000,
        'min_cds': 3_000,
        'min_n50': 50_000,
    },
}

# ── inclusion / exclusion filters ─────────────────────────────────────────────
REQUIRED_ISOLATION_COUNTRY = 'China'

# accepted isolation_source keywords (case-insensitive substring match)
ACCEPTED_SOURCES = [
    'clinical', 'patient', 'blood', 'urine', 'sputum', 'wound',
    'hospital', 'human', 'infant', 'neonatal', 'respiratory',
]

# accepted assembly levels
ACCEPTED_ASSEMBLY_LEVELS = {'Complete Genome', 'Chromosome', 'Scaffold'}

# RefSeq exclusion: suppress/withdrawn assemblies
EXCLUDED_REFSEQ_STATUSES = {'suppressed', 'replaced', 'withdrawn'}

# ── carbapenem resistance gene names (CARD/AMRFinder canonical names) ─────────
CARBAPENEM_GENES = {
    # NDM family
    'NDM-1', 'NDM-2', 'NDM-3', 'NDM-4', 'NDM-5', 'NDM-6', 'NDM-7',
    'NDM-9', 'NDM-12', 'NDM-13', 'NDM-14', 'NDM-15', 'NDM-16', 'NDM-17',
    # KPC family
    'KPC-2', 'KPC-3', 'KPC-4', 'KPC-5', 'KPC-6', 'KPC-7', 'KPC-8',
    # OXA carbapenemases
    'OXA-48', 'OXA-181', 'OXA-232', 'OXA-162', 'OXA-204', 'OXA-244',
    'OXA-23', 'OXA-24', 'OXA-25', 'OXA-26', 'OXA-40', 'OXA-58',
    # IMP family
    'IMP-1', 'IMP-4', 'IMP-6', 'IMP-8', 'IMP-26',
    # VIM family
    'VIM-1', 'VIM-2', 'VIM-4',
}

# GFF3 keyword patterns for annotation-based fallback detection
# (used when AMRFinder TSV is unavailable)
AMR_GFF_KEYWORDS = [
    r'\bNDM-\d',                    # NDM-type MBL alleles (NDM-1, NDM-5...)
    r'\bKPC-\d|blaKPC',             # KPC carbapenemase (gene or allele)
    r'OXA-48', r'OXA-181', r'OXA-232',  # OXA carbapenemases
    r'\bIMP-\d|blaIMP',             # IMP-type MBL alleles (IMP-1, IMP-4...)
    r'\bVIM-\d|blaVIM',             # VIM-type MBL alleles (VIM-1, VIM-2...)
    r'carbapenem.{0,30}(resistance|beta-lactamase)',
    r'bla[A-Z]{2,5}-\d',           # bla gene naming pattern with allele number
]

# ── IS element detection keywords (against GFF product annotations) ───────────
IS_GFF_KEYWORDS = [
    r'IS\d+',           # IS26, IS5, IS1, IS3, etc.
    r'ISKpn\d*',        # Klebsiella-specific IS elements
    r'transposase',
    r'insertion sequence',
    r'Tn\d+',           # transposons
    r'Tn4401',          # KPC-bearing transposon
]

# flanking window: look for IS elements within this distance of an AMR gene
IS_FLANK_WINDOW_BP = 10_000

# ── download settings ─────────────────────────────────────────────────────────
DOWNLOAD_TIMEOUT_S   = 60       # per-file HTTP timeout (seconds)
MAX_RETRIES          = 3        # retry failed downloads
RETRY_DELAY_S        = 5
MAX_GENOMES_CHINA    = None     # P1: fetch all ~677 China clinical genomes
RATE_LIMIT_DELAY_S   = 0.35     # delay between NCBI requests (< 3/s without API key)

# files to download per assembly
DOWNLOAD_FILES = [
    '_genomic.gff.gz',       # GFF3 annotation
    '_genomic.fna.gz',       # genome FASTA
    '_amrfinderplus.tsv',    # NCBI AMRFinder results (may not exist for all)
    '_assembly_stats.txt',   # N50, contig count, genome size
]
