from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Any, Dict, List, Literal, Optional

from .chemistry import analyze_molecule


class Stereocenter(BaseModel):
    atom_index: int
    cip_label: Literal["R", "S", "?"]


class LipinskiRuleResult(BaseModel):
    name: str
    passed: bool
    detail: str


class Alert(BaseModel):
    name: str
    description: str


class RSComparisonEntry(BaseModel):
    description: str
    original_value: str
    flipped_value: str


class OptimizationSuggestion(BaseModel):
    text: str


class HistoricalExample(BaseModel):
    drug: str
    impact: str


class StereochemistryContext(BaseModel):
    num_stereocenters: int
    max_stereoisomers: int
    considerations: List[str]
    note_2d_properties: str
    recommendation: str
    historical_examples: List[HistoricalExample]


# ADMET Response Models
class SolubilityPrediction(BaseModel):
    log_s: float
    solubility_class: str
    method: str


class LipophilicityPrediction(BaseModel):
    clogp: float
    lipophilicity_class: str
    method: str


class BBBPrediction(BaseModel):
    penetrates: bool
    probability: float
    confidence: str
    method: str


class CNSMPOPrediction(BaseModel):
    score: float
    cns_class: str
    component_scores: Dict[str, float]
    method: str


class GIAbsorptionPrediction(BaseModel):
    absorption: str
    bioavailability_score: float
    lipinski_violations: int
    veber_violations: int
    method: str


class ToxicityEndpointResult(BaseModel):
    active: bool
    probability: float
    alerts: List[str]


class ToxicityPrediction(BaseModel):
    endpoints: Dict[str, ToxicityEndpointResult]
    method: str


class CYPPrediction(BaseModel):
    cyp3a4_substrate: bool
    cyp2d6_substrate: bool
    cyp2c9_substrate: bool
    cyp_inhibitor_risk: bool
    method: str


class ClinToxPrediction(BaseModel):
    ct_tox: bool
    probability: float
    structural_alerts: List[str]
    fda_approval_likelihood: str
    method: str


class ADMETPredictions(BaseModel):
    models_available: bool
    model_error: Optional[str] = None
    solubility: Optional[SolubilityPrediction] = None
    lipophilicity: Optional[LipophilicityPrediction] = None
    bbb_penetration: Optional[BBBPrediction] = None
    cns_mpo: Optional[CNSMPOPrediction] = None
    gi_absorption: Optional[GIAbsorptionPrediction] = None
    toxicity: Optional[ToxicityPrediction] = None
    cyp_metabolism: Optional[CYPPrediction] = None
    clinical_toxicity: Optional[ClinToxPrediction] = None


class QEDResult(BaseModel):
    qed_score: Optional[float]
    qed_class: str


class LeadLikeness(BaseModel):
    violations: List[str]
    num_violations: int
    category: str


class ScaffoldAnalysis(BaseModel):
    murcko_scaffold: Optional[str]
    generic_scaffold: Optional[str]
    scaffold_rings: int


class MoleculeAnalysisRequest(BaseModel):
    smiles: str
    # Optional: list of atom indices where R/S should be flipped relative to the input
    flip_centers: Optional[List[int]] = None


class MoleculeAnalysisResponse(BaseModel):
    canonical_smiles: str
    molecule_image: str  # Base64 PNG image data URL
    molecular_weight: float
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    fraction_sp3: float
    ring_count: int
    aromatic_ring_count: int
    lipinski_rules: List[LipinskiRuleResult]
    lipinski_summary: str
    alerts: List[Alert]
    synthetic_accessibility_score: float
    synthetic_accessibility_class: Literal["Easy", "Moderate", "Difficult"]
    flexibility_warnings: List[str]
    complexity_warnings: List[str]
    stereocenters: List[Stereocenter]
    rs_comparisons: List[RSComparisonEntry]
    stereochemistry_context: StereochemistryContext
    admet: ADMETPredictions
    qed: QEDResult
    lead_likeness: LeadLikeness
    scaffold: ScaffoldAnalysis
    fail_fast_score: float
    decision: Literal["Progress", "Optimize", "Kill"]
    decision_rationale: str
    optimization_suggestions: List[OptimizationSuggestion]


app = FastAPI(title="Early Drug Discovery Screening Dashboard API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/analyze", response_model=MoleculeAnalysisResponse)
def analyze(request: MoleculeAnalysisRequest) -> MoleculeAnalysisResponse:
    if not request.smiles or not request.smiles.strip():
        raise HTTPException(status_code=400, detail="SMILES string must not be empty.")

    try:
        result = analyze_molecule(
            smiles=request.smiles,
            flip_centers=request.flip_centers or [],
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e

    return MoleculeAnalysisResponse(**result)

