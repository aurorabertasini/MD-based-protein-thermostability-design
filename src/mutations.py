from __future__ import annotations

from typing import Iterable, List

import numpy as np


def select_candidate_sites(
    rmsf: np.ndarray,
    rmsf_threshold: float = 0.25,
    known_sites: Iterable[int] | None = None,
) -> List[int]:
    """Select residue indices (0-based) with elevated RMSF.

    known_sites can be used to seed the list with literature-supported positions.
    """
    high_motion = np.where(rmsf >= rmsf_threshold)[0].tolist()
    if known_sites:
        combined = sorted(set(high_motion).union(set(known_sites)))
        return combined
    return high_motion
