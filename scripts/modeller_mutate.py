from __future__ import annotations

import argparse
from pathlib import Path

from modeller import Environ, selection
from modeller.scripts import complete_pdb


ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def parse_mutation(code: str, mode: str, default_chain: str) -> tuple[str, int, str]:
    code = code.strip()
    if len(code) < 3:
        raise ValueError(f"Invalid mutation code: {code}")
    if mode == "wt":
        wt = code[0].upper()
        if wt not in ONE_TO_THREE:
            raise ValueError(f"Unknown WT residue code: {wt}")
        rest = code[1:]
        chain_id = default_chain
    else:
        if code[0].isalpha():
            chain_id = code[0].upper()
            rest = code[1:]
        else:
            chain_id = default_chain
            rest = code
    resseq = int(rest[:-1])
    mut = rest[-1].upper()
    if mut not in ONE_TO_THREE:
        raise ValueError(f"Unknown residue code: {mut}")
    return chain_id, resseq, mut


def build_mutant(
    pdb_path: Path,
    mutation: str,
    output_dir: Path,
    mode: str,
    default_chain: str,
) -> Path:
    chain_id, resseq, mut = parse_mutation(mutation, mode, default_chain)

    env = Environ()
    env.io.hetatm = True
    env.io.water = True
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")

    mdl = complete_pdb(env, str(pdb_path))
    target = selection(mdl.residue_range(f"{resseq}:{chain_id}", f"{resseq}:{chain_id}"))
    target.mutate(residue_type=ONE_TO_THREE[mut])

    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = output_dir / "mutated_tmp.pdb"
    mdl.write(file=str(tmp_path))

    # Rebuild topology/atoms for the mutated residue
    mdl2 = complete_pdb(env, str(tmp_path))
    output_path = output_dir / "structure.pdb"
    mdl2.write(file=str(output_path))
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate mutant PDBs with Modeller.")
    parser.add_argument("--pdb", required=True, help="Path to WT PDB file")
    parser.add_argument(
        "--mutations",
        required=True,
        help="Comma-separated mutations like A117V,B23A",
    )
    parser.add_argument(
        "--mutation-format",
        choices=["chain", "wt"],
        default="chain",
        help="Interpret codes as chain+resseq+mut (chain) or wt+resseq+mut (wt).",
    )
    parser.add_argument(
        "--default-chain",
        default="A",
        help="Chain ID to use when --mutation-format=wt or when chain is omitted.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output base directory (mutants will be placed in subfolders).",
    )
    args = parser.parse_args()

    pdb_path = Path(args.pdb).expanduser().resolve()
    out_dir = Path(args.outdir).expanduser().resolve()
    mutations = [m.strip() for m in args.mutations.split(",") if m.strip()]

    for mut in mutations:
        mut_dir = out_dir / mut
        out_path = build_mutant(
            pdb_path,
            mut,
            mut_dir,
            args.mutation_format,
            args.default_chain,
        )
        print(f"{mut}: {out_path}")


if __name__ == "__main__":
    main()
