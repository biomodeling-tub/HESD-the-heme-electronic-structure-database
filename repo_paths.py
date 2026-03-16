from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
TABLES_DIR = REPO_ROOT / "tables"
DERIVED_TABLES_DIR = TABLES_DIR / "derived_tables"
PLOTS_DIR = REPO_ROOT / "plots"
LOGFILES_DIR = REPO_ROOT / "logfiles"
JSONS_DIR = REPO_ROOT / "jsons"
PDB_DIR = REPO_ROOT / "PDB"
PRIOR_ANALYSIS_DIR = REPO_ROOT / "prior_analysis"
QM_ANALYSIS_TABLES_DIR = PRIOR_ANALYSIS_DIR / "qm_analysis_tables"

CANONICAL_ROOT_TABLES = {
    "pyDISH.csv",
    "HESD.csv",
    "HESD_unfiltered.csv",
    "redox_potentials_literature.csv",
}


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def canonical_table_path(filename: str) -> Path:
    return TABLES_DIR / filename


def derived_table_path(filename: str) -> Path:
    ensure_directory(DERIVED_TABLES_DIR)
    return DERIVED_TABLES_DIR / filename


def resolve_table_input(filename: str) -> Path:
    if filename in CANONICAL_ROOT_TABLES:
        return canonical_table_path(filename)

    candidates = [
        derived_table_path(filename),
        canonical_table_path(filename),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def resolve_existing(*paths: Path) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None
