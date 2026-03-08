from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple


def _require_openmm():
    try:
        from openmm.app import (
            PDBFile,
            ForceField,
            Modeller,
            Simulation,
            PME,
            HBonds,
            DCDReporter,
            StateDataReporter,
        )
        from openmm import LangevinIntegrator, Platform
        from openmm.unit import kelvin, picosecond, femtoseconds, nanometer, molar
    except Exception as exc:  # pragma: no cover - import guard
        raise ImportError(
            "OpenMM is required for system preparation and MD. Install openmm before use."
        ) from exc

    return {
        "PDBFile": PDBFile,
        "ForceField": ForceField,
        "Modeller": Modeller,
        "Simulation": Simulation,
        "PME": PME,
        "HBonds": HBonds,
        "DCDReporter": DCDReporter,
        "StateDataReporter": StateDataReporter,
        "LangevinIntegrator": LangevinIntegrator,
        "Platform": Platform,
        "kelvin": kelvin,
        "picosecond": picosecond,
        "femtoseconds": femtoseconds,
        "nanometer": nanometer,
        "molar": molar,
    }


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    try:
        from Bio.PDB import PDBList
    except Exception as exc:  # pragma: no cover - import guard
        raise ImportError("Biopython is required for PDB downloads.") from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = output_dir / f"{pdb_id}.pdb"
    if pdb_path.exists():
        return pdb_path
    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=str(output_dir))
    Path(pdb_file).rename(pdb_path)
    return pdb_path


@dataclass
class PreparedSystem:
    pdb_path: Path
    topology: object
    positions: object
    system: object
    forcefield_files: Tuple[str, ...]



def prepare_system(
    pdb_path: Path,
    forcefield_files: Iterable[str],
    water_model: str = "tip3p",
    padding_nm: float = 1.0,
    ionic_strength_m: float = 0.15,
    output_pdb_path: Path | None = None,
    ph: float = 7.0,
) -> PreparedSystem:
    omm = _require_openmm()
    pdb = omm["PDBFile"](str(pdb_path))
    forcefield = omm["ForceField"](*forcefield_files)
    modeller = omm["Modeller"](pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=ph)
    modeller.addSolvent(
        forcefield,
        model=water_model,
        padding=padding_nm * omm["nanometer"],
        ionicStrength=ionic_strength_m * omm["molar"],
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=omm["PME"],
        constraints=omm["HBonds"],
    )
    if output_pdb_path is not None:
        output_pdb_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_pdb_path, "w", encoding="utf-8") as handle:
            omm["PDBFile"].writeFile(modeller.topology, modeller.positions, handle)
    return PreparedSystem(
        pdb_path=Path(pdb_path),
        topology=modeller.topology,
        positions=modeller.positions,
        system=system,
        forcefield_files=tuple(forcefield_files),
    )


def run_md(
    prep: PreparedSystem,
    output_dir: Path,
    temperature_k: float = 300.0,
    friction_ps: float = 1.0,
    timestep_fs: float = 2.0,
    n_steps: int = 100000,
    report_interval: int = 1000,
    platform_name: str | None = None,
) -> dict:
    omm = _require_openmm()
    output_dir.mkdir(parents=True, exist_ok=True)
    topology_path = output_dir / "topology.pdb"
    trajectory_path = output_dir / "traj.dcd"
    state_path = output_dir / "state.csv"

    integrator = omm["LangevinIntegrator"](
        temperature_k * omm["kelvin"],
        friction_ps / omm["picosecond"],
        timestep_fs * omm["femtoseconds"],
    )
    if platform_name:
        platform = omm["Platform"].getPlatformByName(platform_name)
        simulation = omm["Simulation"](
            prep.topology,
            prep.system,
            integrator,
            platform,
        )
    else:
        simulation = omm["Simulation"](prep.topology, prep.system, integrator)

    simulation.context.setPositions(prep.positions)
    simulation.minimizeEnergy()

    with open(topology_path, "w", encoding="utf-8") as handle:
        omm["PDBFile"].writeFile(prep.topology, prep.positions, handle)

    simulation.reporters.append(
        omm["DCDReporter"](str(trajectory_path), report_interval)
    )
    simulation.reporters.append(
        omm["StateDataReporter"](
            str(state_path),
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=n_steps,
            separator=",",
        )
    )
    simulation.step(n_steps)

    return {
        "topology": topology_path,
        "trajectory": trajectory_path,
        "state": state_path,
    }
