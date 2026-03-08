from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np


def _require_mdtraj():
    try:
        import mdtraj as md
    except Exception as exc:  # pragma: no cover - import guard
        raise ImportError("MDTraj is required for trajectory analysis.") from exc
    return md


def load_trajectory(trajectory_path: Path, topology_path: Path):
    md = _require_mdtraj()
    return md.load(str(trajectory_path), top=str(topology_path))


def _protein_ca_indices(traj):
    ca_indices = traj.topology.select("protein and name CA")
    if ca_indices.size == 0:
        raise ValueError("No protein C-alpha atoms found for analysis.")
    return ca_indices


def _protein_residue_pairs(traj, min_seq_sep: int = 3):
    residues = [res.index for res in traj.topology.residues if res.is_protein]
    pairs = []
    for i, res_i in enumerate(residues):
        for res_j in residues[i + 1 :]:
            if abs(res_j - res_i) > min_seq_sep:
                pairs.append((res_i, res_j))
    if not pairs:
        raise ValueError("No protein residue pairs found for native contacts.")
    return pairs


def _protein_residue_indices(traj):
    return [res.index for res in traj.topology.residues if res.is_protein]


def compute_rmsd(traj) -> np.ndarray:
    md = _require_mdtraj()
    ca_indices = _protein_ca_indices(traj)
    ca_traj = traj.atom_slice(ca_indices)
    return md.rmsd(ca_traj, ca_traj, 0)


def compute_rmsf(traj) -> np.ndarray:
    md = _require_mdtraj()
    ca_indices = _protein_ca_indices(traj)
    ca_traj = traj.atom_slice(ca_indices)
    return md.rmsf(ca_traj, ca_traj, 0)


def compute_native_contacts_q(traj, cutoff_nm: float = 0.45) -> np.ndarray:
    md = _require_mdtraj()
    pairs = _protein_residue_pairs(traj)
    distances, _ = md.compute_contacts(traj, contacts=pairs, scheme="closest-heavy")
    native_mask = distances[0] < cutoff_nm
    if not np.any(native_mask):
        raise ValueError("No native contacts found at the given cutoff.")
    native_distances = distances[:, native_mask]
    return (native_distances < cutoff_nm).mean(axis=1)


def compute_residue_contact_persistence(
    traj, cutoff_nm: float = 0.45, min_seq_sep: int = 3
) -> tuple[list[int], np.ndarray]:
    """Compute per-residue contact persistence from native contacts."""
    md = _require_mdtraj()
    residue_indices = _protein_residue_indices(traj)
    pairs = _protein_residue_pairs(traj, min_seq_sep=min_seq_sep)
    distances, _ = md.compute_contacts(traj, contacts=pairs, scheme="closest-heavy")
    contact_freq = (distances < cutoff_nm).mean(axis=0)

    contact_map = {idx: [] for idx in residue_indices}
    for (res_i, res_j), freq in zip(pairs, contact_freq):
        contact_map[res_i].append(freq)
        contact_map[res_j].append(freq)

    persistence = np.array(
        [float(np.mean(contact_map[idx])) if contact_map[idx] else 0.0 for idx in residue_indices]
    )
    return residue_indices, persistence


def summarize_metrics(
    rmsd: np.ndarray,
    rmsf: np.ndarray,
    q_values: np.ndarray,
) -> Dict[str, float]:
    return {
        "rmsd_mean": float(np.mean(rmsd)),
        "rmsd_std": float(np.std(rmsd)),
        "rmsf_mean": float(np.mean(rmsf)),
        "rmsf_max": float(np.max(rmsf)),
        "q_mean": float(np.mean(q_values)),
        "q_std": float(np.std(q_values)),
    }
