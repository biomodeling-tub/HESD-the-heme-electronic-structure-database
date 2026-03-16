#!/usr/bin/env python3
"""Create planar and isolated distorted heme conformations from a crystal heme PDB."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


CORE_NAME_SET = {
    "NA",
    "NB",
    "NC",
    "ND",
    "CHA",
    "CHB",
    "CHC",
    "CHD",
    *{f"C{i}{ring}" for i in range(1, 5) for ring in "ABCD"},
}
CORE_WITH_FE = set(CORE_NAME_SET) | {"FE"}


@dataclass
class AtomRecord:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    occ: float
    bfac: float
    elem: str
    charge: str
    coord: np.ndarray


@dataclass
class PSFAtom:
    idx: int
    segid: str
    resid: str
    resname: str
    atom_name: str


def parse_pdb(path: Path) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        elem = line[76:78].strip() or name[:1]
        atoms.append(
            AtomRecord(
                record=line[0:6].strip() or "HETATM",
                serial=int(line[6:11]),
                name=name,
                altloc=line[16:17],
                resname=line[17:20].strip() or "HEM",
                chain=line[21:22] or "A",
                resseq=int(line[22:26]),
                icode=line[26:27],
                occ=float(line[54:60] or 1.0),
                bfac=float(line[60:66] or 0.0),
                elem=elem.upper(),
                charge=line[78:80].strip(),
                coord=np.array(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    dtype=float,
                ),
            )
        )
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return atoms


def parse_psf(path: Path) -> tuple[list[PSFAtom], list[tuple[int, int]]]:
    lines = path.read_text(encoding="ascii", errors="ignore").splitlines()

    natom_idx = next((i for i, line in enumerate(lines) if "!NATOM" in line), None)
    if natom_idx is None:
        raise ValueError(f"Could not find !NATOM in PSF file: {path}")
    natom = int(lines[natom_idx].split()[0])

    psf_atoms: list[PSFAtom] = []
    for line in lines[natom_idx + 1 : natom_idx + 1 + natom]:
        parts = line.split()
        if len(parts) < 6:
            continue
        psf_atoms.append(
            PSFAtom(
                idx=int(parts[0]),
                segid=parts[1],
                resid=parts[2],
                resname=parts[3],
                atom_name=parts[4],
            )
        )

    nbond_idx = next((i for i, line in enumerate(lines) if "!NBOND" in line), None)
    if nbond_idx is None:
        raise ValueError(f"Could not find !NBOND in PSF file: {path}")
    nbond = int(lines[nbond_idx].split()[0])

    bond_ints: list[int] = []
    i = nbond_idx + 1
    while len(bond_ints) < 2 * nbond and i < len(lines):
        bond_ints.extend(int(v) for v in lines[i].split())
        i += 1
    if len(bond_ints) < 2 * nbond:
        raise ValueError(f"Incomplete bond list in PSF file: {path}")

    bonds = [(bond_ints[j], bond_ints[j + 1]) for j in range(0, 2 * nbond, 2)]
    return psf_atoms, bonds


def infer_system_paths(input_path: Path) -> tuple[Path, Path]:
    pdb_id = input_path.name.split("_")[0]
    base_dir = input_path.parent

    pdb_candidates = [
        base_dir / f"{pdb_id}_system.pdb",
        base_dir / f"{pdb_id}_system_protonated.pdb",
    ]
    psf_candidates = [
        base_dir / f"{pdb_id}_system.psf",
        base_dir / f"{pdb_id}_system_protonated.psf",
    ]

    system_pdb = next((p for p in pdb_candidates if p.exists()), None)
    system_psf = next((p for p in psf_candidates if p.exists()), None)

    if system_pdb is None or system_psf is None:
        raise FileNotFoundError(
            "Could not infer system PDB/PSF from input path. Please pass --system-pdb and --system-psf."
        )
    return system_pdb, system_psf


def extract_core_bonds_from_system(system_pdb: Path, system_psf: Path) -> list[tuple[str, str]]:
    # Parse system PDB as requested to ensure the same system is loaded in-memory.
    _ = parse_pdb(system_pdb)

    psf_atoms, psf_bonds = parse_psf(system_psf)
    by_idx = {a.idx: a for a in psf_atoms}

    heme_ids = {
        a.idx
        for a in psf_atoms
        if (a.resname.upper() == "HEM" or a.segid.upper() == "HEME") and a.atom_name in CORE_WITH_FE
    }
    if not heme_ids:
        raise ValueError("No HEM core atoms found in system PSF.")

    # In-memory pruning: remove sidechains and axial ligands by only retaining porphyrin core + Fe.
    bond_names: set[tuple[str, str]] = set()
    for i, j in psf_bonds:
        if i not in heme_ids or j not in heme_ids:
            continue
        ni = by_idx[i].atom_name
        nj = by_idx[j].atom_name
        if ni in CORE_WITH_FE and nj in CORE_WITH_FE:
            bond_names.add(tuple(sorted((ni, nj))))

    return sorted(bond_names)


def write_pdb(path: Path, atoms: Iterable[AtomRecord], title: str) -> None:
    with path.open("w", encoding="ascii") as f:
        f.write(f"HEADER    {title[:66]}\n")
        for serial, atom in enumerate(atoms, start=1):
            f.write(
                f"{atom.record:<6}{serial:5d} "
                f"{atom.name:<4}{atom.altloc:1}"
                f"{atom.resname:>3} {atom.chain:1}"
                f"{atom.resseq:4d}{atom.icode:1}   "
                f"{atom.coord[0]:8.3f}{atom.coord[1]:8.3f}{atom.coord[2]:8.3f}"
                f"{atom.occ:6.2f}{atom.bfac:6.2f}          "
                f"{atom.elem:>2}{atom.charge:>2}\n"
            )
        f.write("END\n")


def atom_lookup(atoms: list[AtomRecord]) -> Dict[str, int]:
    return {atom.name: i for i, atom in enumerate(atoms)}


def core_indices(atoms: list[AtomRecord]) -> tuple[list[int], int]:
    idx = atom_lookup(atoms)
    fe_idx = idx.get("FE")
    if fe_idx is None:
        raise ValueError("Input heme does not contain atom name FE.")
    core = [i for i, atom in enumerate(atoms) if atom.name in CORE_NAME_SET]
    if len(core) < 20:
        raise ValueError(
            f"Found only {len(core)} porphyrin core atoms (need >= 20). Atom naming may be non-standard."
        )
    return sorted(set(core + [fe_idx])), fe_idx


def fit_plane(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = coords.mean(axis=0)
    centered = coords - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    normal = normal / np.linalg.norm(normal)
    return center, normal


def plane_basis(center: np.ndarray, normal: np.ndarray, ref_point: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    v = ref_point - center
    v = v - np.dot(v, normal) * normal
    if np.linalg.norm(v) < 1e-8:
        trial = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(trial, normal)) > 0.95:
            trial = np.array([0.0, 1.0, 0.0])
        v = trial - np.dot(trial, normal) * normal
    e1 = v / np.linalg.norm(v)
    e2 = np.cross(normal, e1)
    e2 = e2 / np.linalg.norm(e2)
    return e1, e2


def project_to_plane(atoms: list[AtomRecord], center: np.ndarray, normal: np.ndarray) -> list[AtomRecord]:
    projected: list[AtomRecord] = []
    for atom in atoms:
        rel = atom.coord - center
        planar = atom.coord - np.dot(rel, normal) * normal
        projected.append(AtomRecord(**{**atom.__dict__, "coord": planar}))
    return projected


def core_type(atom_name: str) -> str:
    if atom_name.startswith("C1"):
        return "C1"
    if atom_name.startswith("C2"):
        return "C2"
    if atom_name.startswith("C3"):
        return "C3"
    if atom_name.startswith("C4"):
        return "C4"
    if atom_name.startswith("CH"):
        return "CH"
    if atom_name.startswith("N"):
        return "N"
    return atom_name


def substituent_anchor(atom_name: str) -> str | None:
    if not atom_name or atom_name in CORE_NAME_SET or atom_name == "FE":
        return atom_name if atom_name in CORE_NAME_SET else ("FE" if atom_name == "FE" else None)
    ring = atom_name[-1]
    if ring not in "ABCD":
        return None
    if atom_name.startswith("CM"):
        return f"C3{ring}"
    if atom_name.startswith("CA") or atom_name.startswith("CB"):
        return f"C2{ring}"
    return None


def relax_in_plane_with_bonds(
    core_names: list[str],
    uv_target: np.ndarray,
    fe_name: str,
    bond_name_pairs: list[tuple[str, str]],
    d0_by_pair: Dict[tuple[str, str], float],
    n_iter: int = 2500,
    lr: float = 0.015,
) -> np.ndarray:
    idx = {name: i for i, name in enumerate(core_names)}
    uv = uv_target.copy()
    uv_fe_ref = uv[idx[fe_name]].copy()

    used_pairs = [pair for pair in bond_name_pairs if pair[0] in idx and pair[1] in idx and pair in d0_by_pair]
    if not used_pairs:
        return uv

    w_bond = 2.0
    w_target = 0.18
    w_fe = 5.0

    for _ in range(n_iter):
        grad = np.zeros_like(uv)

        for a, b in used_pairs:
            ia, ib = idx[a], idx[b]
            diff = uv[ia] - uv[ib]
            dist = float(np.linalg.norm(diff)) + 1e-12
            err = dist - d0_by_pair[(a, b)]
            g = 2.0 * w_bond * err * (diff / dist)
            grad[ia] += g
            grad[ib] -= g

        grad += 2.0 * w_target * (uv - uv_target)
        grad[idx[fe_name]] += 2.0 * w_fe * (uv[idx[fe_name]] - uv_fe_ref)

        uv -= lr * grad

    return uv


def regularize_planar_d4h(
    atoms_original: list[AtomRecord],
    atoms_planar: list[AtomRecord],
    core_idx: list[int],
    fe_idx: int,
    center: np.ndarray,
    e1: np.ndarray,
    e2: np.ndarray,
    core_bonds: list[tuple[str, str]],
) -> tuple[list[AtomRecord], Dict[str, float], Dict[str, float]]:
    atoms = [AtomRecord(**{**a.__dict__, "coord": a.coord.copy()}) for a in atoms_planar]
    idx_by_name = atom_lookup(atoms)

    fe = atoms[fe_idx].coord.copy()
    core_non_fe = [i for i in core_idx if i != fe_idx]
    rel = np.array([atoms[i].coord - fe for i in core_non_fe])
    u = rel @ e1
    v = rel @ e2
    phi = np.arctan2(v, u)
    order = np.argsort(phi)

    radii_by_type: Dict[str, List[float]] = {}
    for i in core_non_fe:
        key = core_type(atoms[i].name)
        r = np.linalg.norm(atoms[i].coord - fe)
        radii_by_type.setdefault(key, []).append(r)
    mean_radius = {k: float(np.mean(vals)) for k, vals in radii_by_type.items()}

    step = 2.0 * math.pi / len(core_non_fe)
    phi_sorted = phi[order]
    phi0 = float(np.mean(phi_sorted - np.arange(len(core_non_fe)) * step))

    original_core_coords = {atoms[i].name: atoms[i].coord.copy() for i in core_non_fe}
    core_names = [atoms[i].name for i in core_idx]

    target_uv = np.zeros((len(core_names), 2), dtype=float)
    name_to_target_xyz: Dict[str, np.ndarray] = {}
    for rank, ord_idx in enumerate(order):
        atom_i = core_non_fe[ord_idx]
        atom_name = atoms[atom_i].name
        key = core_type(atom_name)
        this_phi = phi0 + rank * step
        this_r = mean_radius[key]
        xyz = fe + this_r * math.cos(this_phi) * e1 + this_r * math.sin(this_phi) * e2
        name_to_target_xyz[atom_name] = xyz
    name_to_target_xyz["FE"] = fe.copy()

    for i_name, name in enumerate(core_names):
        xyz = name_to_target_xyz[name]
        rel_xyz = xyz - fe
        target_uv[i_name, 0] = float(np.dot(rel_xyz, e1))
        target_uv[i_name, 1] = float(np.dot(rel_xyz, e2))

    d0_by_pair: Dict[tuple[str, str], float] = {}
    orig_idx = atom_lookup(atoms_original)
    for pair in core_bonds:
        if pair[0] not in orig_idx or pair[1] not in orig_idx:
            continue
        d0_by_pair[pair] = float(np.linalg.norm(atoms_original[orig_idx[pair[0]]].coord - atoms_original[orig_idx[pair[1]]].coord))

    relaxed_uv = relax_in_plane_with_bonds(core_names, target_uv, "FE", core_bonds, d0_by_pair)

    phi_map: Dict[str, float] = {}
    r_map: Dict[str, float] = {}
    for i_name, name in enumerate(core_names):
        u_i, v_i = relaxed_uv[i_name]
        new_coord = fe + u_i * e1 + v_i * e2
        atoms[idx_by_name[name]] = AtomRecord(**{**atoms[idx_by_name[name]].__dict__, "coord": new_coord})
        phi_map[name] = float(math.atan2(v_i, u_i)) if name != "FE" else 0.0
        r_map[name] = float(math.hypot(u_i, v_i)) if name != "FE" else 0.0

    # Shift side chains with their anchor to keep attachment geometry consistent after in-plane regularization.
    for i, atom in enumerate(atoms):
        if i in core_idx:
            continue
        anchor = substituent_anchor(atom.name)
        if not anchor or anchor not in idx_by_name or anchor not in original_core_coords:
            continue
        delta = atoms[idx_by_name[anchor]].coord - original_core_coords[anchor]
        atoms[i] = AtomRecord(**{**atom.__dict__, "coord": atom.coord + delta})

    return atoms, phi_map, r_map


def rotate_around_axis(point: np.ndarray, axis_point: np.ndarray, axis_dir: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = axis_dir / np.linalg.norm(axis_dir)
    v = point - axis_point
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    return axis_point + (v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1.0 - c))


def apply_mode(
    atoms_planar: list[AtomRecord],
    core_idx: list[int],
    fe_idx: int,
    normal: np.ndarray,
    phi_map: Dict[str, float],
    r_map: Dict[str, float],
    mode: str,
) -> list[AtomRecord]:
    amp_saddling = 1.25  # stronger saddling
    doming_fe_shift = 0.35
    doming_ring_amp = 0.90
    ruffling_angle = math.radians(27.0)

    atoms = [AtomRecord(**{**a.__dict__, "coord": a.coord.copy()}) for a in atoms_planar]
    idx_by_name = atom_lookup(atoms)

    core_names = {atoms[i].name for i in core_idx}
    planar_core = {name: atoms[idx_by_name[name]].coord.copy() for name in core_names}
    distorted_core = {name: coord.copy() for name, coord in planar_core.items()}
    fe = planar_core["FE"]
    r_max = max(r_map.get(name, 0.0) for name in core_names if name != "FE")

    if mode == "planar":
        pass
    elif mode == "saddling":
        for name in core_names:
            phi = phi_map.get(name, 0.0)
            dz = 0.0 if name == "FE" else amp_saddling * math.sin(2.0 * phi)
            distorted_core[name] = planar_core[name] + dz * normal
    elif mode == "doming":
        n_dz_values: list[float] = []
        for name in core_names:
            if name == "FE":
                dz = doming_fe_shift
            else:
                r = r_map.get(name, 0.0)
                rn = r / r_max if r_max > 0 else 0.0
                dz = doming_ring_amp * (0.45 + 0.55 * (1.0 - rn * rn))
                if name in {"NA", "NB", "NC", "ND"}:
                    n_dz_values.append(dz)
            distorted_core[name] = planar_core[name] + dz * normal

        # Move Fe together with the inner coordinating N-ring.
        if n_dz_values:
            fe_dz = float(sum(n_dz_values) / len(n_dz_values))
            distorted_core["FE"] = planar_core["FE"] + fe_dz * normal
            doming_fe_shift = fe_dz
    elif mode == "ruffling":
        ring_sets = {
            "A": ["NA", "C1A", "C2A", "C3A", "C4A"],
            "B": ["NB", "C1B", "C2B", "C3B", "C4B"],
            "C": ["NC", "C1C", "C2C", "C3C", "C4C"],
            "D": ["ND", "C1D", "C2D", "C3D", "C4D"],
        }
        for ring, names in ring_sets.items():
            n_name = f"N{ring}"
            if n_name not in idx_by_name:
                continue
            axis_point = planar_core[n_name]
            axis_dir = fe - axis_point
            for name in names:
                if name not in distorted_core:
                    continue
                distorted_core[name] = rotate_around_axis(planar_core[name], axis_point, axis_dir, ruffling_angle)
        # Keep bridge methine atoms close to the original porphyrin plane.
        for name in ("CHA", "CHB", "CHC", "CHD"):
            if name in distorted_core:
                distorted_core[name] = planar_core[name].copy()
    else:
        raise ValueError(f"Unknown mode: {mode}")

    for name in core_names:
        i = idx_by_name[name]
        atoms[i] = AtomRecord(**{**atoms[i].__dict__, "coord": distorted_core[name]})

    # Move side chains with their anchor atom in full 3D.
    for i, atom in enumerate(atoms):
        if i in core_idx:
            continue
        anchor = substituent_anchor(atom.name)
        if not anchor or anchor not in distorted_core:
            continue
        anchor_delta = distorted_core[anchor] - planar_core[anchor]
        new_coord = atom.coord + anchor_delta
        if mode == "doming":
            # Outward substituent direction opposite to Fe displacement, with similar magnitude.
            current_shift = float(np.dot(new_coord - atom.coord, normal))
            correction = (-doming_fe_shift - current_shift) * normal
            new_coord = new_coord + correction
        atoms[i] = AtomRecord(**{**atom.__dict__, "coord": new_coord})

    return atoms


def generate_from_input(
    input_path: Path,
    out_dir: Path,
    system_pdb: Path,
    system_psf: Path,
) -> list[Path]:
    atoms = parse_pdb(input_path)
    core_idx, fe_idx = core_indices(atoms)
    core_coords = np.array([atoms[i].coord for i in core_idx if i != fe_idx])

    core_bonds = extract_core_bonds_from_system(system_pdb, system_psf)
    if not core_bonds:
        raise ValueError("No core bond topology extracted from system PSF.")

    center, normal = fit_plane(core_coords)
    fe = atoms[fe_idx].coord
    if np.dot(normal, fe - center) < 0:
        normal = -normal

    ref_idx = atom_lookup(atoms).get("NA", core_idx[0])
    e1, e2 = plane_basis(center, normal, atoms[ref_idx].coord)

    planar_projected = project_to_plane(atoms, center, normal)
    planar_d4h, phi_map, r_map = regularize_planar_d4h(
        atoms,
        planar_projected,
        core_idx,
        fe_idx,
        center,
        e1,
        e2,
        core_bonds,
    )

    outputs = {
        "heme_planar_D4h.pdb": ("planar", "HEME PLANAR D4H"),
        "heme_ruffled.pdb": ("ruffling", "HEME RUFFLING"),
        "heme_saddled.pdb": ("saddling", "HEME SADDLING"),
        "heme_domed.pdb": ("doming", "HEME DOMING"),
    }

    written: list[Path] = []
    for filename, (mode, title) in outputs.items():
        mode_atoms = apply_mode(planar_d4h, core_idx, fe_idx, normal, phi_map, r_map, mode)
        out_path = out_dir / filename
        write_pdb(out_path, mode_atoms, title)
        written.append(out_path)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate planar, ruffled, saddled, and domed heme structures from a crystal heme PDB."
    )
    parser.add_argument(
        "--input",
        default="PDB/1a2f/1a2f_HEME.pdb",
        help="Input heme PDB file (default: PDB/1a2f/1a2f_HEME.pdb).",
    )
    parser.add_argument(
        "--system-pdb",
        default=None,
        help="System PDB used with PSF for covalent connectivity (default: inferred from input).",
    )
    parser.add_argument(
        "--system-psf",
        default=None,
        help="System PSF used for covalent connectivity (default: inferred from input).",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Output directory for generated PDB files (default: current directory).",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.system_pdb is None or args.system_psf is None:
        inferred_pdb, inferred_psf = infer_system_paths(input_path)
        system_pdb = Path(args.system_pdb) if args.system_pdb else inferred_pdb
        system_psf = Path(args.system_psf) if args.system_psf else inferred_psf
    else:
        system_pdb = Path(args.system_pdb)
        system_psf = Path(args.system_psf)

    written = generate_from_input(input_path, out_dir, system_pdb, system_psf)
    print(f"Input: {input_path}")
    print(f"System PDB: {system_pdb}")
    print(f"System PSF: {system_psf}")
    print("Wrote:")
    for path in written:
        print(f"  - {path}")


if __name__ == "__main__":
    main()
