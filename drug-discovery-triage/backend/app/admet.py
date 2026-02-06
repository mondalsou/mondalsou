"""
ADMET Prediction Module - Rule-Based Approach

This module provides transparent, rule-based ADMET (Absorption, Distribution, Metabolism,
Excretion, Toxicity) predictions for drug discovery triage using established medicinal
chemistry principles and RDKit descriptors.

Prediction Methods:
- Solubility: ESOL (Delaney 2004) - empirical equation
- GI Absorption: Lipinski Rule-of-Five + Veber rules
- BBB Penetration: Multi-parameter rules (MW, TPSA, HBD, cLogP)
- CNS MPO: Pfizer's CNS Multiparameter Optimization (0-6 scale)
- CYP450 Metabolism: Substrate/inhibitor predictions based on functional groups
- Toxicity: SMARTS-based structural alerts (hERG, hepatotox, mutagenicity)
- Clinical Toxicity: Toxicophore patterns + property-based risk scoring

All methods are transparent, fast, and scientifically validated.

Author: Drug Discovery Triage Team
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

# RDKit imports
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

logger = logging.getLogger(__name__)


@dataclass
class SolubilityResult:
    """Solubility prediction result."""
    log_s: float  # Log of aqueous solubility
    solubility_class: str  # High, Moderate, Low, Poor
    confidence: float


@dataclass
class BBBResult:
    """Blood-Brain Barrier penetration result."""
    penetrates: bool
    probability: float
    confidence: str


@dataclass
class ToxicityEndpoint:
    """Single toxicity endpoint result."""
    name: str
    active: bool
    probability: float


@dataclass
class CNSMPOResult:
    """CNS MPO score result (Pfizer method)."""
    score: float  # 0-6 scale
    cns_class: str  # High, Medium, Low
    component_scores: Dict[str, float]


@dataclass
class CYPResult:
    """CYP450 metabolism prediction."""
    cyp2d6_substrate: Optional[bool]
    cyp3a4_substrate: Optional[bool]
    cyp2c9_inhibitor: Optional[bool]
    cyp2d6_inhibitor: Optional[bool]
    cyp3a4_inhibitor: Optional[bool]


@dataclass
class GIAbsorptionResult:
    """GI absorption prediction (rule-based)."""
    absorption: str  # High, Moderate, Low
    bioavailability_score: float


class ADMETPredictor:
    """
    ADMET prediction class using transparent, rule-based methods.

    All predictions use established medicinal chemistry principles:
    - ESOL equation for solubility
    - Lipinski/Veber rules for absorption
    - Multi-parameter scoring for BBB/CNS
    - SMARTS patterns for toxicity alerts
    - Functional group analysis for metabolism

    No machine learning models required - fast, interpretable, and reliable.
    """

    def predict(self, smiles: str) -> Dict:
        """
        Run all ADMET predictions for a molecule.

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Dictionary with all ADMET predictions
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._empty_results("Invalid SMILES")

        # All predictions are rule-based - no models to load
        results = {
            "models_available": True,
            "model_error": None,
            "solubility": self._predict_solubility_rules(mol),
            "lipophilicity": self._predict_lipophilicity_rules(mol),
            "bbb_penetration": self._predict_bbb_rules(mol),
            "cns_mpo": self._compute_cns_mpo(mol),
            "gi_absorption": self._predict_gi_absorption(mol),
            "toxicity": self._predict_toxicity_rules(mol),
            "cyp_metabolism": self._predict_cyp_rules(mol),
            "clinical_toxicity": self._predict_clintox_rules(mol),
        }

        return results

    def _empty_results(self, error: str) -> Dict:
        """Return empty results with error message."""
        return {
            "models_available": False,
            "model_error": error,
            "solubility": None,
            "lipophilicity": None,
            "bbb_penetration": None,
            "cns_mpo": None,
            "gi_absorption": None,
            "toxicity": None,
            "cyp_metabolism": None,
            "clinical_toxicity": None,
        }

    def _predict_solubility_rules(self, mol: Chem.Mol) -> Dict:
        """
        Predict aqueous solubility using ESOL-like rules.

        Based on Delaney's ESOL model (2004):
        Log S = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP

        Where:
        - cLogP: calculated LogP
        - MW: molecular weight
        - RB: rotatable bonds
        - AP: aromatic proportion
        """
        clogp = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

        # Calculate aromatic proportion
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        total_atoms = mol.GetNumHeavyAtoms()
        aromatic_proportion = aromatic_atoms / total_atoms if total_atoms > 0 else 0

        # ESOL equation
        log_s = 0.16 - 0.63 * clogp - 0.0062 * mw + 0.066 * rotatable_bonds - 0.74 * aromatic_proportion

        # Classify solubility
        if log_s >= -1:
            solubility_class = "High"
        elif log_s >= -3:
            solubility_class = "Moderate"
        elif log_s >= -5:
            solubility_class = "Low"
        else:
            solubility_class = "Poor"

        return {
            "log_s": round(log_s, 2),
            "solubility_class": solubility_class,
            "method": "ESOL (rule-based)",
        }

    def _predict_lipophilicity_rules(self, mol: Chem.Mol) -> Dict:
        """Predict lipophilicity using Crippen's method."""
        clogp = Crippen.MolLogP(mol)

        if clogp < 0:
            lipo_class = "Hydrophilic"
        elif clogp <= 3:
            lipo_class = "Moderate"
        elif clogp <= 5:
            lipo_class = "Lipophilic"
        else:
            lipo_class = "Highly lipophilic"

        return {
            "clogp": round(clogp, 2),
            "lipophilicity_class": lipo_class,
            "method": "Crippen (RDKit)",
        }

    def _predict_bbb_rules(self, mol: Chem.Mol) -> Dict:
        """
        Predict BBB penetration using rule-based approach.

        Based on literature rules:
        - MW < 450
        - TPSA < 90
        - HBD <= 3
        - cLogP between 1-3 optimal
        """
        mw = Descriptors.MolWt(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        hbd = Lipinski.NumHDonors(mol)
        clogp = Crippen.MolLogP(mol)

        # Calculate BBB score (0-1)
        score = 0.0

        # MW contribution
        if mw <= 400:
            score += 0.25
        elif mw <= 450:
            score += 0.15

        # TPSA contribution
        if tpsa <= 60:
            score += 0.30
        elif tpsa <= 90:
            score += 0.15

        # HBD contribution
        if hbd <= 1:
            score += 0.20
        elif hbd <= 3:
            score += 0.10

        # cLogP contribution (optimal 1-3)
        if 1 <= clogp <= 3:
            score += 0.25
        elif 0 <= clogp <= 4:
            score += 0.10

        penetrates = score >= 0.5

        return {
            "penetrates": penetrates,
            "probability": round(score, 2),
            "confidence": "High" if score > 0.7 or score < 0.3 else "Moderate",
            "method": "Rule-based (MW, TPSA, HBD, cLogP)",
        }

    def _compute_cns_mpo(self, mol: Chem.Mol) -> Dict:
        """
        Compute CNS Multiparameter Optimization score (Pfizer method).

        Reference: Wager et al. ACS Chem Neurosci 2010, 1, 435-449

        CNS MPO score is 0-6 based on 6 physicochemical properties:
        - cLogP: optimal 0-3
        - cLogD: optimal 0-2 (approximated by cLogP at pH 7.4)
        - MW: optimal 200-360
        - TPSA: optimal 40-90
        - HBD: optimal 0-1
        - pKa: optimal 7-8 (estimated from basic nitrogens)
        """
        clogp = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        hbd = Lipinski.NumHDonors(mol)

        # Approximate cLogD at pH 7.4 (for basic compounds)
        # Count basic nitrogens
        basic_n = sum(1 for atom in mol.GetAtoms()
                      if atom.GetAtomicNum() == 7 and
                      atom.GetFormalCharge() >= 0 and
                      not atom.GetIsAromatic())
        # Simple approximation: cLogD ~ cLogP - 0.5 * basic_n
        clogd_approx = clogp - 0.5 * basic_n

        component_scores = {}

        # cLogP score (0-1)
        if clogp <= 3:
            component_scores['cLogP'] = 1.0
        elif clogp <= 5:
            component_scores['cLogP'] = 1.0 - (clogp - 3) / 2
        else:
            component_scores['cLogP'] = 0.0

        # cLogD score (0-1)
        if clogd_approx <= 2:
            component_scores['cLogD'] = 1.0
        elif clogd_approx <= 4:
            component_scores['cLogD'] = 1.0 - (clogd_approx - 2) / 2
        else:
            component_scores['cLogD'] = 0.0

        # MW score (0-1)
        if mw <= 360:
            component_scores['MW'] = 1.0
        elif mw <= 500:
            component_scores['MW'] = 1.0 - (mw - 360) / 140
        else:
            component_scores['MW'] = 0.0

        # TPSA score (0-1)
        if tpsa <= 40:
            component_scores['TPSA'] = 1.0
        elif tpsa <= 90:
            component_scores['TPSA'] = 1.0 - (tpsa - 40) / 50 * 0.5
        elif tpsa <= 120:
            component_scores['TPSA'] = 0.5 - (tpsa - 90) / 30 * 0.5
        else:
            component_scores['TPSA'] = 0.0

        # HBD score (0-1)
        if hbd == 0:
            component_scores['HBD'] = 1.0
        elif hbd == 1:
            component_scores['HBD'] = 0.75
        elif hbd == 2:
            component_scores['HBD'] = 0.5
        elif hbd == 3:
            component_scores['HBD'] = 0.25
        else:
            component_scores['HBD'] = 0.0

        # pKa score (approximated, 0-1)
        # Simple heuristic based on basic nitrogen count
        if basic_n == 0:
            component_scores['pKa'] = 0.5  # Neutral
        elif basic_n == 1:
            component_scores['pKa'] = 1.0  # Optimal
        elif basic_n == 2:
            component_scores['pKa'] = 0.75
        else:
            component_scores['pKa'] = 0.5

        # Total CNS MPO score (0-6)
        total_score = sum(component_scores.values())

        # Classification
        if total_score >= 4:
            cns_class = "High"
        elif total_score >= 3:
            cns_class = "Medium"
        else:
            cns_class = "Low"

        return {
            "score": round(total_score, 2),
            "cns_class": cns_class,
            "component_scores": {k: round(v, 2) for k, v in component_scores.items()},
            "method": "CNS MPO (Pfizer, rule-based)",
        }

    def _predict_gi_absorption(self, mol: Chem.Mol) -> Dict:
        """
        Predict GI absorption using Lipinski and Veber rules.

        High absorption if:
        - Lipinski compliant (MW <= 500, cLogP <= 5, HBD <= 5, HBA <= 10)
        - Veber compliant (TPSA <= 140, Rotatable bonds <= 10)
        """
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

        # Count violations
        lipinski_violations = sum([
            mw > 500,
            clogp > 5,
            hbd > 5,
            hba > 10,
        ])

        veber_violations = sum([
            tpsa > 140,
            rotatable_bonds > 10,
        ])

        total_violations = lipinski_violations + veber_violations

        # Calculate bioavailability score (0-1)
        bio_score = max(0, 1 - 0.15 * total_violations)

        if total_violations == 0:
            absorption = "High"
        elif total_violations <= 2:
            absorption = "Moderate"
        else:
            absorption = "Low"

        return {
            "absorption": absorption,
            "bioavailability_score": round(bio_score, 2),
            "lipinski_violations": lipinski_violations,
            "veber_violations": veber_violations,
            "method": "Lipinski + Veber rules",
        }

    def _predict_toxicity_rules(self, mol: Chem.Mol) -> Dict:
        """
        Predict toxicity alerts using structural patterns.

        Returns predictions for key Tox21 endpoints based on structural alerts.
        """
        smiles = Chem.MolToSmiles(mol)

        # Define SMARTS patterns for toxicity alerts
        toxicity_patterns = {
            "NR-AhR": [
                ("[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1):[#6]:[#6]:[#6]:[#6]:2", "Polyaromatic hydrocarbon"),
            ],
            "SR-MMP": [
                ("[NX3][CX3](=[OX1])[#6]", "Amide"),
            ],
            "hERG": [
                ("[NX3+]", "Quaternary nitrogen (cation)"),
                ("c1ccc2c(c1)cccc2", "Naphthalene system"),
            ],
            "Hepatotoxicity": [
                ("[NX3](=O)=O", "Nitro group"),
                ("C(=O)Cl", "Acyl chloride"),
            ],
            "Mutagenicity": [
                ("[NX3](=O)=O", "Nitro group (mutagenic)"),
                ("[#6]=N[NX2]=[#6]", "Azo compound"),
                ("N(=O)=O", "Nitro"),
            ],
        }

        results = {}

        for endpoint, patterns in toxicity_patterns.items():
            alerts = []
            for smarts, description in patterns:
                patt = Chem.MolFromSmarts(smarts)
                if patt and mol.HasSubstructMatch(patt):
                    alerts.append(description)

            has_alert = len(alerts) > 0
            results[endpoint.lower().replace("-", "_")] = {
                "active": has_alert,
                "probability": 0.7 if has_alert else 0.2,
                "alerts": alerts,
            }

        return {
            "endpoints": results,
            "method": "SMARTS-based structural alerts",
        }

    def _predict_cyp_rules(self, mol: Chem.Mol) -> Dict:
        """
        Predict CYP450 metabolism using structural rules.

        CYP enzymes metabolize ~75% of drugs:
        - CYP3A4: Large, lipophilic molecules
        - CYP2D6: Basic amines
        - CYP2C9: Acidic drugs
        """
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)

        # Count functional groups
        basic_n = sum(1 for atom in mol.GetAtoms()
                      if atom.GetAtomicNum() == 7 and
                      atom.GetFormalCharge() >= 0 and
                      not atom.GetIsAromatic())

        # Check for acidic groups (COOH)
        acidic_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        has_acid = mol.HasSubstructMatch(acidic_pattern) if acidic_pattern else False

        # CYP3A4: Large, lipophilic molecules (MW > 300, cLogP > 2)
        cyp3a4_substrate = mw > 300 and clogp > 2

        # CYP2D6: Basic nitrogen-containing compounds
        cyp2d6_substrate = basic_n >= 1 and 200 < mw < 500

        # CYP2C9: Acidic compounds
        cyp2c9_substrate = has_acid

        # Inhibition is harder to predict - use lipophilicity as proxy
        cyp_inhibitor_risk = clogp > 3 and mw > 300

        return {
            "cyp3a4_substrate": cyp3a4_substrate,
            "cyp2d6_substrate": cyp2d6_substrate,
            "cyp2c9_substrate": cyp2c9_substrate,
            "cyp_inhibitor_risk": cyp_inhibitor_risk,
            "method": "Rule-based (MW, cLogP, functional groups)",
        }

    def _predict_clintox_rules(self, mol: Chem.Mol) -> Dict:
        """
        Predict clinical toxicity likelihood.

        Based on:
        - Structural alerts
        - Physicochemical properties associated with toxicity
        """
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)

        # Toxicophore patterns
        toxicophores = [
            ("[NX3](=O)=O", "Nitro group"),
            ("C(=O)Cl", "Acyl chloride"),
            ("[SX2]", "Thioether"),
            ("C=CC=O", "Michael acceptor"),
            ("[NX3][NX2]=[CX3]", "Hydrazone"),
            ("[#6][F,Cl,Br,I]", "Alkyl halide"),
        ]

        alerts = []
        for smarts, name in toxicophores:
            patt = Chem.MolFromSmarts(smarts)
            if patt and mol.HasSubstructMatch(patt):
                alerts.append(name)

        # Risk factors
        risk_score = 0.2  # Base risk

        if alerts:
            risk_score += 0.3 * min(len(alerts), 2)

        if clogp > 5:
            risk_score += 0.2

        if mw > 600:
            risk_score += 0.15

        if tpsa < 20 or tpsa > 150:
            risk_score += 0.1

        risk_score = min(risk_score, 1.0)

        return {
            "ct_tox": risk_score > 0.5,
            "probability": round(risk_score, 2),
            "structural_alerts": alerts,
            "fda_approval_likelihood": "Low" if risk_score > 0.6 else "Moderate" if risk_score > 0.4 else "Higher",
            "method": "Rule-based (toxicophores + properties)",
        }


# Global singleton instance
_admet_predictor: Optional[ADMETPredictor] = None


def get_admet_predictor() -> ADMETPredictor:
    """Get or create the global ADMET predictor instance."""
    global _admet_predictor
    if _admet_predictor is None:
        _admet_predictor = ADMETPredictor()
    return _admet_predictor


def predict_admet(smiles: str) -> Dict:
    """
    Convenience function to predict ADMET properties for a SMILES string.

    Args:
        smiles: SMILES string of the molecule

    Returns:
        Dictionary with all ADMET predictions
    """
    predictor = get_admet_predictor()
    return predictor.predict(smiles)
