# Python-based End-to-End Protein Thermostability Design via MD

This repository contains an end-to-end, notebook-driven pipeline to **design and screen protein variants with improved thermostability** using molecular dynamics (MD), with **T4 lysozyme** as the case study.

The core value of the project is that the pipeline is:
- **fully rerunnable from start to end**,
- **data-generating at every stage** (files are saved to structured folders),
- **report-ready**, with results, figures, and methodological rationale already presented in the notebook.
- **portable to other proteins** by changing the protein input (PDB ID or structure path) and output paths.

## What This Project Does

The notebook implements a complete workflow:
1. Download and prepare WT structure.
2. Run WT MD simulations.
3. Compute stability-oriented metrics (RMSD, RMSF, native contacts Q).
4. Propose mutation candidates using RMSF + SASA + environment-aware rules.
5. Build mutants (Modeller), run short screening MD, compute paired DeltaQ.
6. Validate top candidates with longer MD runs and produce final ranking.

All generated simulation artifacts are saved under `data/` and reused by downstream steps.

Although this repository showcases the workflow on T4 lysozyme, the same pipeline logic can be reused for other proteins with minimal input changes.

## Repository Structure

- `notebooks/`: main end-to-end notebook (`T4_lysozyme_thermostability.ipynb`), runnable from start to end and already populated with narrative results.
- `src/`: reusable Python modules for MD setup, analysis, and mutation utilities.
- `scripts/`: helper scripts (Modeller-based mutant generation).
- `data/`: generated inputs/outputs for WT and mutant simulations.
- `pictures/`: workflow figures and static project images.
- `requirements.txt`: Python dependencies for the notebook and modules.

## Environment Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Notes:
- Some steps depend on heavy scientific packages (OpenMM, MDTraj, MDAnalysis, Biopython).
- Mutant generation uses **Modeller** via `scripts/modeller_mutate.py` (often run in a dedicated conda environment, as shown in the notebook).

## How To Run

From repository root:

```bash
jupyter notebook
```

Open:

`notebooks/T4_lysozyme_thermostability.ipynb`

Then run all cells from top to bottom.

## Output and Data Flow

The pipeline writes intermediate and final artifacts into `data/`, including:
- prepared structures,
- trajectories (`traj.dcd`),
- topology/state files,
- mutation screening tables,
- long-run validation summaries.

Downstream notebook cells read these files back automatically, so the workflow is reproducible and traceable.

## Reproducibility Notes

- The notebook is designed to run sequentially and regenerate the full analysis.
- Existing outputs can be reused to avoid recomputation.
- Runtime depends heavily on hardware (CPU/GPU) and MD simulation length.
