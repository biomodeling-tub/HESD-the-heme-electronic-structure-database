# HESD: The Heme Electronic Structure Database

This repository combines structure preparation, quantum-chemistry parsing, and downstream analysis for heme proteins.

The current codebase supports two practical reproducibility levels:

1. Rebuild the analysis tables and plots from the checked-in intermediates in `logfiles/`, `jsons/`, `PDB/`, and `tables/pyDISH.csv`.
2. Re-run the expensive structure-preparation and quantum-chemistry workflow, which additionally requires external software such as CHARMM, XTB, Gaussian16, Open Babel, and a SLURM-based HPC environment.

Run all commands from the repository root unless stated otherwise.

## What Is In This Checkout

In the checkout inspected on March 16, 2026, the repository already contains:

- `tables/pyDISH.csv`: 1515 metadata/distortion rows used as the main heme annotation table.
- `logfiles/`: 2892 Gaussian log files.
- `jsons/`: 2890 parsed JSON files derived from Gaussian logs.
- `tables/HESD_unfiltered.csv`: 2890 rows x 6522 columns.
- `tables/HESD.csv`: 1800 rows x 449 columns.
- `tables/derived_tables/processed_output.csv`: 728 rows x 81 columns.
- `tables/derived_tables/reduction_potentials.csv`: 197 rows x 7 columns.
- `PDB/`: 939 per-structure directories plus extracted PDB/XYZ assets.

## Environment

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate heme
```

The Python analysis stack is defined in `environment.yml`. The full end-to-end workflow also expects external programs that are not installed by conda here:

- `CHARMM`
- `xtb`
- `Gaussian16`
- `Open Babel` or `obabel`
- `SLURM` commands such as `sbatch` and `squeue`

## Reproducibility Map

### Rebuildable from checked-in data

These outputs can be regenerated from files already present in this repository:

- `tables/HESD_unfiltered.csv` from `jsons/`
- `tables/derived_tables/processed_output.csv` and related processed tables from `tables/HESD_unfiltered.csv`
- `tables/derived_tables/reduction_potentials.csv`
- `prior_analysis/qm_analysis_tables/*.csv`
- `tables/derived_tables/rmsd*.csv`
- `tables/derived_tables/iron_axial_distances.csv`
- `tables/derived_tables/iron_plane_distances.csv`
- plots generated from those tables

### Requires external QC/HPC software

These stages depend on external programs and cluster-specific setup:

- heme extraction and protonation in `scripts/heme_struct.py`
- unconstrained XTB optimization in `scripts/run_xtb_unconstrained.py`
- Gaussian input generation in `automation/beb00042_folders_files_xtb.py`
- Gaussian job orchestration in `automation/pipeline_xtb.sh` or `scripts/gaussian_pipeline_manager.sh`

### Archived outputs with incomplete current provenance

These files are present, but the current top-level scripts do not emit them exactly as checked in:

- `tables/HESD.csv`: preserved snapshot, not written by a current top-level script
- `tables/derived_tables/*`: archived generated outputs; the active code writes to `tables/`
- `tables/pyDISH.csv`: consumed by the pipeline, but the repo does not expose a complete raw-to-`pyDISH.csv` generator

## Important Path Caveats

The active pipeline modules now resolve paths relative to the repository root. Legacy absolute paths still exist in some archived or debug-only scripts, but the main rebuild commands below no longer depend on `/home/pbuser/Desktop/PhD_WORK/heme`.

Two practical consequences:

- Prefer passing explicit paths to constructors such as `g16parser(...)`, `JSONProcessor(...)`, and `RMSDAnalyzer(...)`.
- Treat `tables/` as the home of canonical input/main-product tables and `tables/derived_tables/` as the home of generated downstream analysis outputs.

## Minimal Rebuild From Checked-In Data

### 1. Recreate the flattened database from JSON

This rebuilds `tables/HESD_unfiltered.csv` from `jsons/`.

```bash
python - <<'PY'
from jsonprocessor import JSONProcessor

JSONProcessor(
    json_dir="jsons",
    output_csv="tables/HESD_unfiltered.csv",
).process_files()
PY
```

### 2. Recreate the processed analysis tables

This is the actively maintained downstream table generator. It writes:

- `tables/derived_tables/processed_output.csv`
- `tables/derived_tables/processed_output_unsorted_axial.csv`
- `tables/derived_tables/processed_output_sorted_axial.csv`
- `tables/derived_tables/processed_output_homo_lumo_all.csv` if requested
- `tables/derived_tables/reduction_potentials.csv`
- `prior_analysis/qm_analysis_tables/*.csv`

```bash
mkdir -p prior_analysis

python - <<'PY'
import pandas as pd
from preprocessor import DataPreprocessor

df = pd.read_csv("tables/HESD_unfiltered.csv", low_memory=False)
processor = DataPreprocessor(
    df,
    write_file=True,
    skip_plots=True,
    analysis_folder="prior_analysis",
    homo_lumo_all=True,
)
processor.process()
PY
```

### 3. Recreate structural comparison tables

This step uses the prepared structures in `PDB/` and the freshly written `tables/derived_tables/processed_output.csv`.

```bash
python - <<'PY'
from calculate_rmsd import RMSDAnalyzer

analyzer = RMSDAnalyzer(pdb_base_dir="PDB", verbose=True)
analyzer.run_all_analyses(
    rmsd=True,
    axial=True,
    iron_plane=True,
    absolute_iron_plane=True,
    distortion=False,
)
PY
```

Expected runtime outputs include:

- `tables/derived_tables/rmsd_results.csv`
- `tables/derived_tables/rmsd.csv`
- `tables/derived_tables/iron_axial_distances.csv`
- `tables/derived_tables/iron_plane_distances.csv`
- `tables/derived_tables/iron_plane_distances_plots.csv`
- `tables/derived_tables/suspicious_distances.log`
- plot files in `plots/`

### 4. Recreate the LaTeX summary tables

```bash
python create_latex_table.py
```

This script reads `tables/derived_tables/processed_output.csv` by default and writes LaTeX tables into `tables/derived_tables/`.

### 5. Recreate selected analysis plots

Exploratory QM analysis tables are created in step 2. Additional plot scripts can then be run, for example:

```bash
python plot_qm_analysis.py
python scripts/plot_energy_vs_homo_lumo_gap.py
python scripts/print_homo_lumo_statistics.py
```

## Re-Parsing Gaussian Logs

If you want to rebuild `jsons/` directly from `logfiles/`:

```bash
python - <<'PY'
from g16parser import g16parser

parser = g16parser(
    log_dir="logfiles",
    json_dir="jsons",
    csv_path="tables/pyDISH.csv",
)
parser.parse()
PY
```

This is the correct entry point for the current parser. Do not rely on the class defaults unless your checkout also lives at the legacy `.../heme` path.

## Full End-to-End Regeneration

Use this only if you want to regenerate structures and quantum calculations, not just the final tables.

### A. Prepare heme structures and CHARMM/XTB inputs

```bash
python scripts/heme_struct.py
```

This script downloads PDB entries, isolates hemes, creates protonated systems, and prepares XTB inputs. It depends on `tables/pyDISH.csv`, internet access for PDB retrieval, and external programs.

Subset runs are now explicit via `--row-start`; the default behavior is the full table.

### B. Run unconstrained XTB optimization

```bash
python scripts/run_xtb_unconstrained.py
```

For charge-multiplicity exploration:

```bash
python scripts/run_xtb_unconstrained.py --mode charge-mult --pdb-list pdb_ids.txt
```

### C. Generate Gaussian input directories from XTB results

Create `pdb_ids.txt` with one PDB ID per line, then run:

```bash
cd automation
python beb00042_folders_files_xtb.py
```

This creates `automation_xtb/<pdb>/<pdb><state>/` job folders and Gaussian input scripts based on the latest XTB-optimized geometries found under `PDB/<pdb>/`.

### D. Submit and monitor Gaussian jobs

For the XTB-based workflow:

```bash
cd automation
bash pipeline_xtb.sh
```

For the older non-XTB automation layout:

```bash
bash scripts/gaussian_pipeline_manager.sh
```

Both orchestration scripts assume SLURM, `sbatch`, `squeue`, and cluster-specific job scripts.

### E. Parse logs and rebuild analysis tables

After Gaussian jobs finish, return to the parser and analysis steps:

1. `g16parser.py`
2. `jsonprocessor.py`
3. `preprocessor.py`
4. `calculate_rmsd.py`
5. plotting and table scripts

## Working With The Checked-In Snapshots

The four canonical/input tables remain in `tables/`:

- `tables/pyDISH.csv`
- `tables/HESD.csv`
- `tables/HESD_unfiltered.csv`
- `tables/redox_potentials_literature.csv`

Fresh analysis runs now write derived outputs to `tables/derived_tables/` by default. Active scripts resolve those derived files automatically, so no manual copying should be needed.

## Tests

The repository contains many script-style regression tests in `tests/`. The path- and table-related tests that previously failed during collection were updated alongside the reproducibility fixes. Some of the broader plotting tests can still be slow because they generate real figures and touch large data files.

## File Guide

- `scripts/heme_struct.py`: structure download, heme isolation, CHARMM/XTB preparation
- `scripts/orient_heme.py`: porphyrin orientation logic
- `scripts/run_xtb_unconstrained.py`: XTB optimization driver
- `automation/beb00042_folders_files_xtb.py`: Gaussian input/job generation from XTB results
- `automation/pipeline_xtb.sh`: SLURM job orchestration for the XTB-based workflow
- `g16parser.py`: Gaussian log parser
- `jsonprocessor.py`: JSON flattening into tabular data
- `preprocessor.py`: feature filtering, encoding, derived properties, reduction potentials
- `calculate_rmsd.py`: RMSD, Fe-axial, and Fe-plane analyses
- `create_latex_table.py`: LaTeX table generation
- `plot_qm_analysis.py` and `plots.py`: plotting utilities

## Recommended Reproduction Order

If your goal is to reproduce the published analysis from files already present in the repo, use this order:

1. `jsons/` -> `tables/HESD_unfiltered.csv`
2. `tables/HESD_unfiltered.csv` -> `tables/derived_tables/processed_output.csv` and `tables/derived_tables/reduction_potentials.csv`
3. `PDB/` + `tables/derived_tables/processed_output.csv` -> RMSD and Fe distance tables
4. processed tables -> LaTeX tables and plots

If your goal is to regenerate the raw quantum-chemistry data as well, prepend the CHARMM/XTB/Gaussian workflow before step 1.
