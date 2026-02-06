from __future__ import annotations

import base64
import logging
from dataclasses import dataclass
from io import BytesIO
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, Lipinski, rdMolDescriptors, Scaffolds

from .admet import predict_admet

logger = logging.getLogger(__name__)

# Import QED calculator
try:
    from rdkit.Chem import QED
    HAS_QED = True
except ImportError:
    HAS_QED = False
    logger.warning("QED module not available")


@dataclass
class AlertPattern:
    name: str
    smarts: str
    description: str


ALERT_PATTERNS: List[AlertPattern] = [
    # Expanded PAINS and liability patterns
    AlertPattern(
        name="Catechol",
        smarts="c1cc(O)c(O)cc1",
        description="Catechols are redox-active and may chelate metals, causing assay interference.",
    ),
    AlertPattern(
        name="Quinone",
        smarts="O=C1C=CC(=O)C=C1",
        description="Quinones are Michael acceptors and redox-active, potential toxicophores.",
    ),
    AlertPattern(
        name="Michael acceptor (α,β-unsaturated carbonyl)",
        smarts="C=CC=O",
        description="Michael acceptors undergo covalent addition with nucleophiles (proteins, DNA).",
    ),
    AlertPattern(
        name="Rhodanine",
        smarts="S=C1SC(=O)NC1",
        description="Rhodanines are frequent hitters in HTS, known PAINS motif.",
    ),
    AlertPattern(
        name="Nitro aromatic",
        smarts="[NX3](=O)=O-[c]",
        description="Aromatic nitro groups are often toxic and metabolically labile.",
    ),
    AlertPattern(
        name="Aniline (unsubstituted)",
        smarts="Nc1ccccc1",
        description="Anilines can be oxidized to quinone imines, causing toxicity.",
    ),
    AlertPattern(
        name="Aldehyde",
        smarts="[#6][CX3H](=O)",
        description="Aldehydes are reactive electrophiles, may form Schiff bases.",
    ),
    AlertPattern(
        name="Thiocarbonyl",
        smarts="C=S",
        description="Thiocarbonyls are reactive and metabolically unstable.",
    ),
    AlertPattern(
        name="Hydrazine",
        smarts="[NX3][NX3]",
        description="Hydrazines are known mutagens and carcinogens.",
    ),
    AlertPattern(
        name="Acyl hydrazide",
        smarts="C(=O)N[NX3]",
        description="Acyl hydrazides can interfere with metal-dependent assays.",
    ),
]


def _mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    Chem.SanitizeMol(mol)
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def _flip_stereocenters(mol: Chem.Mol, flip_centers: List[int]) -> Chem.Mol:
    mol = Chem.Mol(mol)
    for idx in flip_centers:
        atom = mol.GetAtomWithIdx(idx)
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def _compute_basic_descriptors(mol: Chem.Mol) -> Dict:
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    aromatic_ring_count = sum(
        1 for ring in ring_info.BondRings() if all(mol.GetBondWithIdx(b).GetIsAromatic() for b in ring)
    )

    sp3_carbons = 0
    total_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            total_carbons += 1
            if atom.GetHybridization() == Chem.HybridizationType.SP3:
                sp3_carbons += 1
    fraction_sp3 = float(sp3_carbons) / float(total_carbons) if total_carbons > 0 else 0.0

    return {
        "molecular_weight": mw,
        "clogp": clogp,
        "tpsa": tpsa,
        "hbd": int(hbd),
        "hba": int(hba),
        "rotatable_bonds": int(rot_bonds),
        "fraction_sp3": fraction_sp3,
        "ring_count": int(ring_count),
        "aromatic_ring_count": int(aromatic_ring_count),
    }


def _evaluate_lipinski(props: Dict) -> Tuple[List[Dict], str]:
    mw = props["molecular_weight"]
    clogp = props["clogp"]
    hbd = props["hbd"]
    hba = props["hba"]

    results: List[Dict] = []

    rules = [
        ("Molecular weight ≤ 500", mw <= 500, f"MW = {mw:.1f}"),
        ("cLogP ≤ 5", clogp <= 5, f"cLogP = {clogp:.2f}"),
        ("HBD ≤ 5", hbd <= 5, f"HBD = {hbd}"),
        ("HBA ≤ 10", hba <= 10, f"HBA = {hba}"),
    ]

    passed_count = 0
    for name, passed, detail in rules:
        if passed:
            passed_count += 1
        results.append({"name": name, "passed": passed, "detail": detail})

    if passed_count == 4:
        summary = "Pass"
    elif passed_count == 3:
        summary = "Borderline"
    else:
        summary = "Fail"

    return results, summary


def _find_alerts(mol: Chem.Mol, props: Dict) -> List[Dict]:
    alerts: List[Dict] = []
    for pattern in ALERT_PATTERNS:
        patt = Chem.MolFromSmarts(pattern.smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            alerts.append({"name": pattern.name, "description": pattern.description})

    # Heuristic for excessively lipophilic aromatic systems
    if props["clogp"] > 5 and props["aromatic_ring_count"] >= 3:
        alerts.append(
            {
                "name": "Highly lipophilic aromatic system",
                "description": "High cLogP with multiple aromatic rings suggests permeability and safety risk.",
            }
        )

    return alerts


def _synthetic_accessibility_heuristic(mol: Chem.Mol, props: Dict) -> Tuple[float, str]:
    """
    Lightweight, transparent heuristic inspired by SAS but not predictive.
    Returns a score on 1 (easy) – 10 (difficult) and a class label.
    """
    score = 3.0

    # More rings and stereocenters increase difficulty
    ring_count = props["ring_count"]
    rot_bonds = props["rotatable_bonds"]

    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    n_chiral = len(chiral_centers)

    score += 0.5 * max(0, ring_count - 2)
    score += 0.5 * max(0, n_chiral - 1)
    score += 0.2 * max(0, rot_bonds - 5)

    # Penalize unusual sizes
    mw = props["molecular_weight"]
    if mw > 600:
        score += 1.0
    if mw > 800:
        score += 1.0

    score = max(1.0, min(10.0, score))

    if score <= 4:
        classification = "Easy"
    elif score <= 7:
        classification = "Moderate"
    else:
        classification = "Difficult"

    return float(score), classification


def _flexibility_complexity_warnings(mol: Chem.Mol, props: Dict) -> Tuple[List[str], List[str]]:
    flex_warnings: List[str] = []
    comp_warnings: List[str] = []

    rot_bonds = props["rotatable_bonds"]
    if rot_bonds >= 10:
        flex_warnings.append("High flexibility (≥10 rotatable bonds) → potential entropic penalty on binding.")
    elif rot_bonds >= 7:
        flex_warnings.append("Moderate flexibility (7–9 rotatable bonds) → consider rigidifying the scaffold.")

    fsp3 = props["fraction_sp3"]
    if fsp3 < 0.2:
        comp_warnings.append("Very low Fsp³ (<0.2) → flat, aromatic-heavy scaffold.")
    elif fsp3 < 0.3:
        comp_warnings.append("Low Fsp³ (0.2–0.3) → consider adding 3D character.")

    ring_count = props["ring_count"]
    aromatic_ring_count = props["aromatic_ring_count"]
    if aromatic_ring_count >= 3:
        comp_warnings.append("Multiple aromatic rings (≥3) → watch for solubility and safety liabilities.")
    if ring_count == 0:
        comp_warnings.append("No rings → may be too flexible and metabolically labile.")

    return flex_warnings, comp_warnings


def _stereocenters(mol: Chem.Mol) -> List[Dict]:
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    stereocenters: List[Dict] = []
    for idx, label in centers:
        stereocenters.append(
            {
                "atom_index": idx,
                "cip_label": label if label in ("R", "S") else "?",
            }
        )
    return stereocenters


def _fail_fast_score(
    lipinski_summary: str,
    alerts: List[Dict],
    sas_score: float,
    flex_warnings: List[str],
    comp_warnings: List[str],
) -> Dict:
    """
    Compute transparent heuristic fail-fast score (0–100).
    35% Lipinski, 25% alerts, 20% SAS, 10% flexibility, 10% complexity.
    """
    # Lipinski (35)
    if lipinski_summary == "Pass":
        lip_contrib = 35.0
    elif lipinski_summary == "Borderline":
        lip_contrib = 25.0
    else:
        lip_contrib = 10.0

    # Alerts (25) – simple penalty for each alert
    n_alerts = len(alerts)
    if n_alerts == 0:
        alert_contrib = 25.0
    elif n_alerts == 1:
        alert_contrib = 18.0
    elif n_alerts == 2:
        alert_contrib = 10.0
    else:
        alert_contrib = 5.0

    # SAS (20) – map 1–10 to 20→0 linearly
    sas_contrib = max(0.0, 20.0 * (10.0 - (sas_score - 1.0)) / 9.0)

    # Flexibility (10)
    if not flex_warnings:
        flex_contrib = 10.0
    elif len(flex_warnings) == 1:
        flex_contrib = 6.0
    else:
        flex_contrib = 3.0

    # Complexity (10)
    if not comp_warnings:
        comp_contrib = 10.0
    elif len(comp_warnings) == 1:
        comp_contrib = 6.0
    else:
        comp_contrib = 3.0

    total = lip_contrib + alert_contrib + sas_contrib + flex_contrib + comp_contrib

    if total >= 70:
        decision = "Progress"
    elif total >= 40:
        decision = "Optimize"
    else:
        decision = "Kill"

    rationale_bits: List[str] = [
        f"Lipinski contribution: {lip_contrib:.1f}/35 ({lipinski_summary}).",
        f"Alerts contribution: {alert_contrib:.1f}/25 (n={n_alerts}).",
        f"Synthetic accessibility contribution: {sas_contrib:.1f}/20 (score {sas_score:.1f}).",
        f"Flexibility contribution: {flex_contrib:.1f}/10 (warnings: {len(flex_warnings)}).",
        f"Complexity contribution: {comp_contrib:.1f}/10 (warnings: {len(comp_warnings)}).",
    ]

    return {
        "fail_fast_score": float(total),
        "decision": decision,
        "decision_rationale": " ".join(rationale_bits),
    }


def _stereochemistry_considerations(mol: Chem.Mol) -> Dict:
    """
    Provide educational context about stereochemistry impact.

    2D topological properties (MW, cLogP, TPSA) are IDENTICAL between enantiomers.
    This function explains what DOES differ and why stereochemistry matters.
    """
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)

    # Count defined stereocenters (R or S)
    defined_centers = [c for c in centers if c[1] in ("R", "S")]
    num_stereocenters = len(defined_centers)

    # Maximum possible stereoisomers = 2^n (ignoring meso compounds)
    max_stereoisomers = 2 ** num_stereocenters if num_stereocenters > 0 else 1

    # Scientific considerations that differ between enantiomers
    considerations = []

    if num_stereocenters > 0:
        considerations = [
            "Biological activity: Enantiomers often show 10-1000x activity differences",
            "Receptor binding: Each enantiomer may bind different targets or with different affinity",
            "Metabolism: CYP450 enzymes (e.g., CYP2D6, CYP3A4) are stereospecific",
            "Pharmacokinetics: Half-life can vary 2-10x between enantiomers",
            "Toxicity: One enantiomer may be toxic while the other is therapeutic (e.g., thalidomide)",
            "Drug-drug interactions: Stereoselective inhibition/induction of metabolic enzymes",
        ]

    # Important note about 2D properties
    note_2d = (
        "2D topological properties (MW, cLogP, TPSA, HBD/HBA) are IDENTICAL between "
        "stereoisomers because they depend only on atom connectivity, not 3D arrangement."
    )

    # Generate recommendation based on stereocenter count
    if num_stereocenters == 0:
        recommendation = "No stereocenters present. Achiral molecule."
    elif num_stereocenters == 1:
        recommendation = (
            "Single stereocenter: If SAR data shows activity, strongly consider preparing "
            "and testing pure enantiomers early in development."
        )
    elif num_stereocenters == 2:
        recommendation = (
            f"Two stereocenters ({max_stereoisomers} possible stereoisomers): "
            "Stereoselective synthesis or chiral resolution will likely be needed. "
            "Test individual stereoisomers to identify the eutomer (active form)."
        )
    else:
        recommendation = (
            f"Multiple stereocenters ({num_stereocenters} centers, up to {max_stereoisomers} stereoisomers): "
            "Complex stereochemistry may present significant synthetic challenges. "
            "Consider whether all stereocenters are required for activity."
        )

    # Historical examples for context
    examples = [
        {
            "drug": "Thalidomide",
            "impact": "(R)-enantiomer is sedative; (S)-enantiomer causes teratogenicity"
        },
        {
            "drug": "Omeprazole/Esomeprazole",
            "impact": "(S)-omeprazole (esomeprazole) shows improved pharmacokinetics"
        },
        {
            "drug": "Ibuprofen",
            "impact": "(S)-ibuprofen is the active anti-inflammatory; (R) is largely inactive"
        },
    ]

    return {
        "num_stereocenters": num_stereocenters,
        "max_stereoisomers": max_stereoisomers,
        "considerations": considerations,
        "note_2d_properties": note_2d,
        "recommendation": recommendation,
        "historical_examples": examples if num_stereocenters > 0 else [],
    }


def _rs_comparisons(mol: Chem.Mol) -> List[Dict]:
    """
    For each stereocenter, flip only that center and compare basic descriptors.

    NOTE: 2D properties will be identical - this is scientifically correct!
    The comparison is kept for visual structure comparison purposes.
    """
    base_props = _compute_basic_descriptors(mol)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    entries: List[Dict] = []

    for idx, label in centers:
        if label not in ("R", "S"):
            continue
        flipped = _flip_stereocenters(mol, [idx])
        flipped_props = _compute_basic_descriptors(flipped)

        for metric, pretty in [
            ("clogp", "cLogP"),
            ("tpsa", "TPSA"),
            ("molecular_weight", "Molecular weight"),
            ("rotatable_bonds", "Rotatable bonds"),
        ]:
            entries.append(
                {
                    "description": f"Atom {idx} {label}→"
                    f"{'R' if label == 'S' else 'S'}: {pretty}",
                    "original_value": f"{base_props[metric]:.2f}" if isinstance(base_props[metric], float) else str(
                        base_props[metric]
                    ),
                    "flipped_value": f"{flipped_props[metric]:.2f}"
                    if isinstance(flipped_props[metric], float)
                    else str(flipped_props[metric]),
                }
            )

    return entries


def _generate_molecule_image(mol: Chem.Mol, width: int = 400, height: int = 300) -> str:
    """
    Generate a PNG image of the molecule and return as base64 string.
    """
    # Generate 2D coordinates if not present
    AllChem.Compute2DCoords(mol)

    # Draw the molecule
    drawer = Draw.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.addAtomIndices = False
    opts.bondLineWidth = 2
    opts.setBackgroundColour((0.06, 0.09, 0.16, 1.0))  # Dark background #0f172a

    # Set atom colors for dark theme using updateAtomPalette
    # Dictionary maps atomic number to (R, G, B) tuple
    atom_palette = {
        6: (0.9, 0.9, 0.9),      # C - light gray
        7: (0.23, 0.51, 0.96),   # N - blue
        8: (0.94, 0.27, 0.27),   # O - red
        16: (0.92, 0.70, 0.03),  # S - yellow
        9: (0.13, 0.77, 0.37),   # F - green
        17: (0.13, 0.77, 0.37),  # Cl - green
        35: (0.66, 0.33, 0.97),  # Br - purple
        53: (0.66, 0.33, 0.97),  # I - purple
        1: (0.9, 0.9, 0.9),      # H - light gray
        15: (1.0, 0.5, 0.0),     # P - orange
    }
    opts.updateAtomPalette(atom_palette)

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Get PNG bytes and encode as base64
    png_bytes = drawer.GetDrawingText()
    b64_str = base64.b64encode(png_bytes).decode('utf-8')

    return f"data:image/png;base64,{b64_str}"


def _optimization_suggestions(
    props: Dict,
    alerts: List[Dict],
    sas_class: str,
    flex_warnings: List[str],
    comp_warnings: List[str],
) -> List[Dict]:
    suggestions: List[str] = []

    clogp = props["clogp"]
    if clogp > 4.5:
        suggestions.append(
            "High cLogP → consider introducing heteroatoms or heterocycles to reduce lipophilicity."
        )

    if flex_warnings:
        suggestions.append(
            "High flexibility → explore ring closure or conformational constraints to reduce entropy penalty."
        )

    if any(a["name"] == "Michael acceptor" for a in alerts):
        suggestions.append(
            "Michael acceptor alert → consider replacing the electrophilic motif with a less reactive bioisostere."
        )

    if sas_class == "Difficult":
        suggestions.append(
            "High synthetic complexity → look for ways to reduce stereocenters or simplify fused ring systems."
        )

    if any("low Fsp³" in w.lower() or "very low fsp³" in w.lower() for w in comp_warnings):
        suggestions.append(
            "Flat scaffold → consider adding sp³ centers or small rings to increase 3D character."
        )

    # De-duplicate while preserving order
    seen = set()
    deduped: List[Dict] = []
    for text in suggestions:
        if text not in seen:
            seen.add(text)
            deduped.append({"text": text})

    return deduped


def _compute_qed(mol: Chem.Mol) -> Dict:
    """
    Compute Quantitative Estimate of Drug-likeness (QED).

    QED is a weighted score (0-1) based on 8 molecular properties.
    Reference: Bickerton et al., Nature Chemistry 2012
    """
    if not HAS_QED:
        return {"qed_score": None, "qed_class": "N/A"}

    try:
        qed_score = QED.qed(mol)

        # Classification based on QED score
        if qed_score >= 0.67:
            qed_class = "Highly drug-like"
        elif qed_score >= 0.49:
            qed_class = "Moderately drug-like"
        elif qed_score >= 0.34:
            qed_class = "Poorly drug-like"
        else:
            qed_class = "Not drug-like"

        return {
            "qed_score": float(qed_score),
            "qed_class": qed_class,
        }
    except Exception as e:
        logger.warning(f"QED calculation failed: {e}")
        return {"qed_score": None, "qed_class": "N/A"}


def _evaluate_lead_likeness(props: Dict) -> Dict:
    """
    Evaluate lead-likeness using Rule of 3 (for fragment/lead compounds).

    Reference: Congreve et al., Drug Discovery Today 2003
    Criteria: MW 150-350, cLogP ≤ 3, HBD ≤ 3, HBA ≤ 3, Rotatable bonds ≤ 3
    """
    mw = props["molecular_weight"]
    clogp = props["clogp"]
    hbd = props["hbd"]
    hba = props["hba"]
    rot_bonds = props["rotatable_bonds"]

    violations = []

    if mw < 150 or mw > 350:
        violations.append(f"MW {mw:.1f} (should be 150-350)")
    if clogp > 3:
        violations.append(f"cLogP {clogp:.2f} (should be ≤ 3)")
    if hbd > 3:
        violations.append(f"HBD {hbd} (should be ≤ 3)")
    if hba > 3:
        violations.append(f"HBA {hba} (should be ≤ 3)")
    if rot_bonds > 3:
        violations.append(f"Rotatable bonds {rot_bonds} (should be ≤ 3)")

    num_violations = len(violations)

    if num_violations == 0:
        category = "Lead-like"
    elif num_violations <= 2:
        category = "Borderline lead-like"
    else:
        category = "Not lead-like"

    return {
        "violations": violations,
        "num_violations": num_violations,
        "category": category,
    }


def _analyze_scaffold(mol: Chem.Mol) -> Dict:
    """
    Extract Murcko scaffold and generic framework.

    Reference: Bemis & Murcko, J. Med. Chem. 1996
    """
    try:
        # Get Murcko scaffold (keeps heteroatoms)
        scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold else None

        # Get generic scaffold (all atoms as C, all bonds as single)
        generic_scaffold = Scaffolds.MurckoScaffold.MakeScaffoldGeneric(scaffold) if scaffold else None
        generic_smiles = Chem.MolToSmiles(generic_scaffold) if generic_scaffold else None

        # Count rings in scaffold
        scaffold_rings = scaffold.GetRingInfo().NumRings() if scaffold else 0

        return {
            "murcko_scaffold": scaffold_smiles,
            "generic_scaffold": generic_smiles,
            "scaffold_rings": scaffold_rings,
        }
    except Exception as e:
        logger.warning(f"Scaffold analysis failed: {e}")
        return {
            "murcko_scaffold": None,
            "generic_scaffold": None,
            "scaffold_rings": 0,
        }


def analyze_molecule(smiles: str, flip_centers: List[int]) -> Dict:
    mol = _mol_from_smiles(smiles)
    if flip_centers:
        mol = _flip_stereocenters(mol, flip_centers)

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

    props = _compute_basic_descriptors(mol)
    lipinski_rules, lip_summary = _evaluate_lipinski(props)
    alerts = _find_alerts(mol, props)
    sas_score, sas_class = _synthetic_accessibility_heuristic(mol, props)
    flex_warnings, comp_warnings = _flexibility_complexity_warnings(mol, props)
    stereocenters = _stereocenters(mol)
    rs_comparisons = _rs_comparisons(mol)
    stereochemistry_context = _stereochemistry_considerations(mol)
    ff = _fail_fast_score(lip_summary, alerts, sas_score, flex_warnings, comp_warnings)
    suggestions = _optimization_suggestions(props, alerts, sas_class, flex_warnings, comp_warnings)

    # New properties
    qed = _compute_qed(mol)
    lead_like = _evaluate_lead_likeness(props)
    scaffold = _analyze_scaffold(mol)

    # Generate molecule image
    molecule_image = _generate_molecule_image(mol, width=500, height=350)

    # Get ADMET predictions
    try:
        admet_results = predict_admet(canonical_smiles)
    except Exception as e:
        logger.warning(f"ADMET prediction failed: {e}")
        admet_results = {
            "models_available": False,
            "model_error": str(e),
            "solubility": None,
            "lipophilicity": None,
            "bbb_penetration": None,
            "cns_mpo": None,
            "gi_absorption": None,
            "toxicity": None,
            "cyp_metabolism": None,
            "clinical_toxicity": None,
        }

    return {
        "canonical_smiles": canonical_smiles,
        "molecule_image": molecule_image,
        **props,
        "lipinski_rules": lipinski_rules,
        "lipinski_summary": lip_summary,
        "alerts": alerts,
        "synthetic_accessibility_score": sas_score,
        "synthetic_accessibility_class": sas_class,
        "flexibility_warnings": flex_warnings,
        "complexity_warnings": comp_warnings,
        "stereocenters": stereocenters,
        "rs_comparisons": rs_comparisons,
        "stereochemistry_context": stereochemistry_context,
        "admet": admet_results,
        "qed": qed,
        "lead_likeness": lead_like,
        "scaffold": scaffold,
        "fail_fast_score": ff["fail_fast_score"],
        "decision": ff["decision"],
        "decision_rationale": ff["decision_rationale"],
        "optimization_suggestions": suggestions,
    }

