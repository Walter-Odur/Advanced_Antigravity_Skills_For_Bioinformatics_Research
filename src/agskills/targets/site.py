"""Binding-site definition for docking.

This module fills a hole in the original skills: ``dock`` documented a
``--site binding_site.json`` input and the SLURM generator defaulted to a
``binding_site.json`` path, but **no subcommand anywhere produced that
file**. A researcher following the documented workflow reached the docking
step with no way to define the search box, and the only alternative,
``--center``, crashed (the parser defined ``--size`` while the function read
``args.box_size``).

Three ways to define a site are provided, in descending order of
reliability:

1. :func:`site_from_hetatm` - the centroid of a co-crystallised ligand.
   This is the gold standard when a holo structure is available.
2. :func:`site_from_residues` - the centroid of named pocket residues, for
   when the site is known from mutagenesis or homology.
3. :func:`site_from_structure` - a box enclosing the whole protein, for
   blind docking. Honest but expensive and much less accurate; it is
   labelled as such in the output.
"""

from __future__ import annotations

import statistics
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

from ..errors import InvalidInputError
from ..io_utils import require_file
from .clean import BIOLOGICAL_IONS, COMMON_ADDITIVES, SOLVENT_RESNAMES

__all__ = [
    "BindingSite",
    "Atom",
    "COFACTORS",
    "parse_atoms",
    "list_ligands",
    "site_from_hetatm",
    "site_from_residues",
    "site_from_structure",
]

#: Cofactors, prosthetic groups, sugars and lipids that are part of the
#: protein's own machinery rather than a candidate drug.
#:
#: Distinguishing these matters more than it might appear. Ranking
#: heterogens purely by size picks the cofactor in any flavoenzyme,
#: kinase or dehydrogenase: in PDB 6HEZ, *Mycobacterium tuberculosis*
#: DprE1, the FAD cofactor has 53 atoms while the actual inhibitor (ligand
#: 0SK) has 28. A box centred on FAD points at the flavin site, not the
#: inhibitor pocket, so every subsequent docking result would be against
#: the wrong cavity.
COFACTORS = frozenset({
    # Flavins and nicotinamides
    "FAD", "FMN", "FMA", "RBF",
    "NAD", "NAI", "NAP", "NDP", "NAX", "NDE",
    # Nucleotides and analogues
    "ATP", "ADP", "AMP", "ANP", "AGS", "ACP", "APC",
    "GTP", "GDP", "GMP", "GNP", "GSP",
    "CTP", "CDP", "CMP", "UTP", "UDP", "UMP", "TTP", "TMP",
    "FAD", "COA", "ACO", "CAA", "SAM", "SAH", "TPP", "TDP",
    # Porphyrins and metal clusters
    "HEM", "HEC", "HEA", "HEB", "HAS", "SRM", "CLA", "BCL", "PHO",
    "SF4", "FES", "F3S", "FE2", "MOS", "ICS",
    # Vitamins and derivatives
    "PLP", "PMP", "BTN", "B12", "COB", "CBY", "MQ7", "UQ1", "UQ2",
    "ASC", "THM", "FOL", "MTX",
    # Glutathione and redox
    "GSH", "GDS", "GTT", "TRS",
    # Sugars and glycans
    "NAG", "NDG", "BMA", "MAN", "GAL", "GLC", "BGC", "FUC", "XYS",
    "SIA", "SLB", "A2G", "GLA", "RAM", "ADA",
    # Lipids, detergents and membrane mimetics
    "PLM", "MYR", "STE", "OLA", "OLB", "OLC", "PEE", "PCW", "PC1",
    "LMT", "LDA", "DDQ", "C8E", "BOG", "HEZ", "D10", "D12", "CPS",
    "CHD", "CHS", "CLR", "Y01",
    # Polyamines and buffers that sometimes slip through
    "SPD", "SPM", "PUT", "BTB", "POP", "PPV",
})


def is_cofactor(resname: str) -> bool:
    """Whether a residue name is a known cofactor or prosthetic group.

    >>> is_cofactor("FAD")
    True
    >>> is_cofactor("0SK")
    False
    """
    return resname.strip().upper() in COFACTORS

#: Vina's search box is padded around the reference ligand so that poses
#: larger than the reference still fit. 4 A either side is the common
#: default in the AutoDock literature.
DEFAULT_PADDING_A = 4.0
#: Vina will run with a very small box, but a box under this size cannot
#: accommodate a drug-sized ligand in any orientation.
MIN_BOX_A = 12.0
#: Beyond roughly this size, exhaustiveness has to rise steeply to keep the
#: search converged, so we warn.
LARGE_BOX_A = 30.0


@dataclass(frozen=True)
class Atom:
    """One parsed PDB atom record."""

    record: str
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    element: str

    @property
    def is_hetatm(self) -> bool:
        return self.record == "HETATM"


@dataclass
class BindingSite:
    """A docking search box."""

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    method: str
    reference: str = ""
    confidence: str = ""
    n_reference_atoms: int = 0
    warnings: list[str] = field(default_factory=list)
    selected_ligand: dict[str, Any] | None = None
    alternative_ligands: list[dict[str, Any]] = field(default_factory=list)

    @property
    def center(self) -> tuple[float, float, float]:
        return (self.center_x, self.center_y, self.center_z)

    @property
    def size(self) -> tuple[float, float, float]:
        return (self.size_x, self.size_y, self.size_z)

    @property
    def volume(self) -> float:
        return self.size_x * self.size_y * self.size_z

    def as_dict(self) -> dict[str, Any]:
        """Serialise with the flat key names ``dock --site`` consumes."""
        data = asdict(self)
        for key in ("center_x", "center_y", "center_z",
                    "size_x", "size_y", "size_z"):
            data[key] = round(float(data[key]), 3)
        data["volume_A3"] = round(self.volume, 1)
        return data


def parse_atoms(pdb_text: str, *, first_model_only: bool = True) -> list[Atom]:
    """Parse ATOM and HETATM records from PDB text.

    A deliberately small, dependency-free parser: binding-site geometry
    needs coordinates and residue identity only, and avoiding BioPython here
    keeps the site tools usable in a minimal install.
    """
    atoms: list[Atom] = []
    for line in pdb_text.splitlines():
        record = line[:6]
        if first_model_only and record == "ENDMDL":
            break
        if record not in ("ATOM  ", "HETATM"):
            continue
        try:
            atoms.append(Atom(
                record=record.strip(),
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resname=line[17:20].strip().upper(),
                chain=(line[21].strip() or "_"),
                resseq=int(line[22:26]),
                icode=line[26].strip(),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                element=(line[76:78].strip().upper() if len(line) >= 78 else ""),
            ))
        except (ValueError, IndexError):
            continue
    return atoms


def _is_candidate_ligand(resname: str) -> bool:
    name = resname.strip().upper()
    return (
        name not in SOLVENT_RESNAMES
        and name not in COMMON_ADDITIVES
        and name not in BIOLOGICAL_IONS
        and name not in {"UNK", "UNX", "UNL"}
    )


def list_ligands(pdb_path: str | Path) -> list[dict[str, Any]]:
    """List candidate co-crystallised ligands, best drug candidate first.

    Solvent, buffer additives and lone ions are excluded, since none marks a
    druggable pocket. The remainder are ranked so that a **drug-like**
    ligand outranks a cofactor, and only then by heavy-atom count:

    1. Not a known cofactor, 8-80 heavy atoms (a plausible small molecule).
    2. Not a known cofactor, any size.
    3. A known cofactor (reported, but never preferred).

    Ranking by size alone selects the cofactor in any flavoenzyme or
    kinase; see :data:`COFACTORS`.
    """
    text = require_file(pdb_path).read_text(encoding="utf-8", errors="replace")
    groups: dict[tuple[str, str, int, str], list[Atom]] = {}
    for atom in parse_atoms(text):
        if not atom.is_hetatm or not _is_candidate_ligand(atom.resname):
            continue
        if atom.element == "H":
            continue
        groups.setdefault(
            (atom.resname, atom.chain, atom.resseq, atom.icode), []
        ).append(atom)

    ligands: list[dict[str, Any]] = []
    for key, atoms in groups.items():
        resname, chain, resseq, icode = key
        cofactor = is_cofactor(resname)
        n_atoms = len(atoms)
        drug_sized = 8 <= n_atoms <= 80
        if cofactor:
            tier, role = 2, "cofactor or prosthetic group"
        elif drug_sized:
            tier, role = 0, "drug-like ligand"
        else:
            tier, role = 1, ("fragment or ion-sized heterogen"
                             if n_atoms < 8 else
                             "large heterogen (peptide, polymer or glycan)")
        ligands.append({
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "n_atoms": n_atoms,
            "label": f"{resname}:{chain}:{resseq}{icode}".rstrip(),
            "is_cofactor": cofactor,
            "role": role,
            "_tier": tier,
        })

    ligands.sort(key=lambda e: (e["_tier"], -e["n_atoms"], e["resname"],
                                e["chain"], e["resseq"]))
    for entry in ligands:
        entry.pop("_tier", None)
    return ligands


def _box_from_atoms(atoms: list[Atom], padding: float,
                    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    xs = [a.x for a in atoms]
    ys = [a.y for a in atoms]
    zs = [a.z for a in atoms]
    center = (statistics.fmean(xs), statistics.fmean(ys), statistics.fmean(zs))
    size = tuple(
        max(MIN_BOX_A, (max(v) - min(v)) + 2 * padding) for v in (xs, ys, zs)
    )
    return center, size  # type: ignore[return-value]


def site_from_hetatm(pdb_path: str | Path, *, resname: str | None = None,
                     chain: str | None = None, resseq: int | None = None,
                     padding: float = DEFAULT_PADDING_A) -> BindingSite:
    """Define the box from a co-crystallised ligand.

    With no selector, the largest candidate heterogen is used.

    Raises:
        InvalidInputError: If the structure has no candidate ligand, or the
            requested selector matches nothing.
    """
    path = require_file(pdb_path)
    text = path.read_text(encoding="utf-8", errors="replace")
    atoms = [a for a in parse_atoms(text)
             if a.is_hetatm and _is_candidate_ligand(a.resname)
             and a.element != "H"]
    if not atoms:
        available = sorted({
            a.resname for a in parse_atoms(text) if a.is_hetatm
        })
        raise InvalidInputError(
            f"{path.name} has no co-crystallised ligand to centre a box on.",
            hint=(
                f"Heterogens present: {', '.join(available) or 'none'}. "
                "Define the site from pocket residues instead "
                "(--residues), or use a homologous holo structure."
            ),
        )

    selected = atoms
    if resname:
        selected = [a for a in selected if a.resname == resname.strip().upper()]
    if chain:
        selected = [a for a in selected if a.chain == chain.strip()]
    if resseq is not None:
        selected = [a for a in selected if a.resseq == resseq]

    if not selected:
        options = ", ".join(entry["label"] for entry in list_ligands(path)[:10])
        raise InvalidInputError(
            "No heterogen matched the requested ligand selector.",
            hint=f"Candidates in this structure: {options}",
        )

    chose_automatically = resname is None and chain is None and resseq is None
    alternatives: list[dict[str, Any]] = []
    picked: dict[str, Any] | None = None

    if chose_automatically:
        # Narrow to one ligand rather than averaging over every heterogen,
        # which would centre the box on empty space between two sites.
        candidates = list_ligands(path)
        picked = candidates[0]
        alternatives = candidates[1:]
        selected = [
            a for a in selected
            if a.resname == picked["resname"] and a.chain == picked["chain"]
            and a.resseq == picked["resseq"]
        ]
        reference = picked["label"]
    else:
        first = selected[0]
        reference = f"{first.resname}:{first.chain}:{first.resseq}{first.icode}".rstrip()

    center, size = _box_from_atoms(selected, padding)
    confidence = ("high - the box is centred on an experimentally observed "
                  "bound ligand")
    if picked is not None and picked.get("is_cofactor"):
        confidence = ("medium - the only heterogens present are cofactors, so "
                      "the box is centred on one of those rather than on a "
                      "drug-like ligand")
    site = BindingSite(
        center_x=center[0], center_y=center[1], center_z=center[2],
        size_x=size[0], size_y=size[1], size_z=size[2],
        method="co-crystallised ligand centroid",
        reference=reference,
        confidence=confidence,
        n_reference_atoms=len(selected),
    )
    if picked is not None:
        site.selected_ligand = picked
        site.alternative_ligands = alternatives
        if picked.get("is_cofactor"):
            site.warnings.append(
                f"{picked['label']} is a {picked['role']}, not a drug-like "
                "ligand. Its site is where the protein binds its own "
                "cofactor, which is usually not the pocket you want to dock "
                "into. Select the intended ligand explicitly with "
                "--ligand-resname."
            )
        if alternatives:
            others = ", ".join(
                f"{a['label']} ({a['n_atoms']} atoms, {a['role']})"
                for a in alternatives[:6]
            )
            site.warnings.append(
                f"{len(alternatives)} other heterogen(s) are present: "
                f"{others}. If one of those is your target ligand, re-run "
                "with --ligand-resname."
            )
    _add_box_warnings(site)
    return site


def site_from_residues(pdb_path: str | Path, residues: Iterable[str], *,
                       padding: float = DEFAULT_PADDING_A) -> BindingSite:
    """Define the box from named pocket residues.

    Each selector is ``CHAIN:RESSEQ`` (for example ``A:230``) or a bare
    residue number, which then matches in any chain.

    Raises:
        InvalidInputError: If no selector matches, or none is parseable.
    """
    path = require_file(pdb_path)
    atoms = parse_atoms(path.read_text(encoding="utf-8", errors="replace"))

    wanted: set[tuple[str | None, int]] = set()
    bad: list[str] = []
    for selector in residues:
        token = str(selector).strip()
        if not token:
            continue
        try:
            if ":" in token:
                chain_part, num_part = token.split(":", 1)
                wanted.add((chain_part.strip() or None, int(num_part)))
            else:
                wanted.add((None, int(token)))
        except ValueError:
            bad.append(token)
    if bad:
        raise InvalidInputError(
            f"Could not parse residue selector(s): {', '.join(bad)}",
            hint="Use CHAIN:RESSEQ, for example A:230, or a bare number.",
        )
    if not wanted:
        raise InvalidInputError("No residues were supplied.")

    selected = [
        a for a in atoms
        if any((chain is None or a.chain == chain) and a.resseq == num
               for chain, num in wanted)
    ]
    if not selected:
        raise InvalidInputError(
            "None of the requested residues were found in the structure.",
            hint="Check the numbering matches this file; PDB author "
                 "numbering and sequence position often differ.",
        )

    matched = sorted({f"{a.chain}:{a.resseq}" for a in selected})
    center, size = _box_from_atoms(selected, padding)
    site = BindingSite(
        center_x=center[0], center_y=center[1], center_z=center[2],
        size_x=size[0], size_y=size[1], size_z=size[2],
        method="pocket residue centroid",
        reference=", ".join(matched),
        confidence="medium - depends on the residue selection being correct",
        n_reference_atoms=len(selected),
    )
    if len(matched) < len(wanted):
        site.warnings.append(
            f"Only {len(matched)} of {len(wanted)} requested residues were "
            "found; the box is centred on the ones that matched."
        )
    _add_box_warnings(site)
    return site


def site_from_structure(pdb_path: str | Path, *, padding: float = 2.0) -> BindingSite:
    """Enclose the entire protein, for blind docking.

    Included for completeness and labelled low-confidence: a whole-protein
    box needs far higher exhaustiveness and still samples poorly compared
    with a focused box.
    """
    path = require_file(pdb_path)
    atoms = [a for a in parse_atoms(path.read_text(encoding="utf-8", errors="replace"))
             if not a.is_hetatm]
    if not atoms:
        raise InvalidInputError(f"{path.name} contains no protein atoms")
    center, size = _box_from_atoms(atoms, padding)
    site = BindingSite(
        center_x=center[0], center_y=center[1], center_z=center[2],
        size_x=size[0], size_y=size[1], size_z=size[2],
        method="whole-structure bounding box (blind docking)",
        reference=path.name,
        confidence="low - no pocket information was used",
        n_reference_atoms=len(atoms),
        warnings=[
            "Blind docking over a whole protein samples the search space "
            "sparsely. Raise --exhaustiveness substantially, and prefer a "
            "ligand- or residue-defined box whenever one is available.",
        ],
    )
    _add_box_warnings(site)
    return site


def _add_box_warnings(site: BindingSite) -> None:
    """Attach advisory warnings about box dimensions."""
    if max(site.size) > LARGE_BOX_A:
        site.warnings.append(
            f"The longest box edge is {max(site.size):.1f} A. Boxes above "
            f"{LARGE_BOX_A:.0f} A need a higher --exhaustiveness for the "
            "search to converge."
        )
    if min(site.size) <= MIN_BOX_A:
        site.warnings.append(
            f"A box edge was clamped up to the {MIN_BOX_A:.0f} A minimum; the "
            "reference ligand is small relative to a typical drug."
        )
