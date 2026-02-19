import { useState } from 'react'
import {
  AppBar,
  Box,
  Button,
  Chip,
  CircularProgress,
  Container,
  IconButton,
  LinearProgress,
  Paper,
  Stack,
  Tab,
  Tabs,
  TextField,
  Toolbar,
  Tooltip,
  Typography,
} from '@mui/material'
import Grid from '@mui/material/Grid'
import './App.css'

// Icons as simple SVG components
const ScienceIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M19.8 18.4L14 10.67V6.5l1.35-1.69c.26-.33.03-.81-.39-.81H9.04c-.42 0-.65.48-.39.81L10 6.5v4.17L4.2 18.4c-.49.66-.02 1.6.8 1.6h14c.82 0 1.29-.94.8-1.6z" />
  </svg>
)

const MoleculeIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <circle cx="7" cy="7" r="2.5" />
    <circle cx="17" cy="7" r="2.5" />
    <circle cx="12" cy="17" r="2.5" />
    <line x1="9" y1="8" x2="15" y2="8" stroke="currentColor" strokeWidth="1.5" />
    <line x1="8" y1="9" x2="11" y2="15" stroke="currentColor" strokeWidth="1.5" />
    <line x1="16" y1="9" x2="13" y2="15" stroke="currentColor" strokeWidth="1.5" />
  </svg>
)

const SpeedIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M20.38 8.57l-1.23 1.85a8 8 0 01-.22 7.58H5.07A8 8 0 0115.58 6.85l1.85-1.23A10 10 0 003 17a2 2 0 001.73 3h14.54A2 2 0 0021 17a10 10 0 00-.62-8.43z" />
    <path d="M10.59 15.41a2 2 0 002.83 0l5.66-8.49-8.49 5.66a2 2 0 000 2.83z" />
  </svg>
)

const AlertIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M12 2L1 21h22L12 2zm0 3.99L19.53 19H4.47L12 5.99zM11 10v4h2v-4h-2zm0 6v2h2v-2h-2z" />
  </svg>
)

const ChiralIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8zm-1-4h2v2h-2zm0-2h2V7h-2z" />
  </svg>
)

const LightbulbIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M9 21c0 .55.45 1 1 1h4c.55 0 1-.45 1-1v-1H9v1zm3-19C8.14 2 5 5.14 5 9c0 2.38 1.19 4.47 3 5.74V17c0 .55.45 1 1 1h6c.55 0 1-.45 1-1v-2.26c1.81-1.27 3-3.36 3-5.74 0-3.86-3.14-7-7-7z" />
  </svg>
)

const CopyIcon = () => (
  <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
    <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z" />
  </svg>
)

const CompareIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M10 3H5c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h5v2h2V1h-2v2zm0 15H5l5-6v6zm9-15h-5v2h5v13l-5-6v9h5c1.1 0 2-.9 2-2V5c0-1.1-.9-2-2-2z" />
  </svg>
)

// ADMET Icon (pill/capsule shape)
const ADMETIcon = () => (
  <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
    <path d="M4.22 11.29l5.49-5.49A5.24 5.24 0 0117.5 2.5a5.24 5.24 0 013.79 8.79l-5.49 5.49A5.24 5.24 0 018 20a5.24 5.24 0 01-3.78-8.71zm1.41 1.41A3.24 3.24 0 008 18a3.24 3.24 0 004.29-1.3l2.92-2.92-5.66-5.66-2.92 2.92A3.24 3.24 0 005.63 12.7z" />
  </svg>
)

// Brain icon for CNS (currently unused)
// const BrainIcon = () => (
//   <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
//     <path d="M12 2a9 9 0 00-9 9c0 4.17 2.84 7.67 6.69 8.69L12 22l2.31-2.31C18.16 18.67 21 15.17 21 11a9 9 0 00-9-9zm0 16c-3.87 0-7-3.13-7-7s3.13-7 7-7 7 3.13 7 7-3.13 7-7 7zm-1-11h2v6h-2zm0 8h2v2h-2z" />
//   </svg>
// )

// Liver icon for metabolism (currently unused)
// const LiverIcon = () => (
//   <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
//     <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-2 15l-5-5 1.41-1.41L10 14.17l7.59-7.59L19 8l-9 9z" />
//   </svg>
// )

// Toxicity icon (skull, currently unused)
// const ToxicityIcon = () => (
//   <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
//     <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm1 15h-2v-2h2v2zm0-4h-2V7h2v6z" />
//   </svg>
// )

type Stereocenter = {
  atom_index: number
  cip_label: 'R' | 'S' | '?'
}

type LipinskiRuleResult = {
  name: string
  passed: boolean
  detail: string
}

type Alert = {
  name: string
  description: string
}

type RSComparisonEntry = {
  description: string
  original_value: string
  flipped_value: string
}

type OptimizationSuggestion = {
  text: string
}

type HistoricalExample = {
  drug: string
  impact: string
}

type StereochemistryContext = {
  num_stereocenters: number
  max_stereoisomers: number
  considerations: string[]
  note_2d_properties: string
  recommendation: string
  historical_examples: HistoricalExample[]
}

// ADMET Types
type SolubilityPrediction = {
  log_s: number
  solubility_class: string
  method: string
}

type LipophilicityPrediction = {
  clogp: number
  lipophilicity_class: string
  method: string
}

type BBBPrediction = {
  penetrates: boolean
  probability: number
  confidence: string
  method: string
}

type CNSMPOPrediction = {
  score: number
  cns_class: string
  component_scores: Record<string, number>
  method: string
}

type GIAbsorptionPrediction = {
  absorption: string
  bioavailability_score: number
  lipinski_violations: number
  veber_violations: number
  method: string
}

type ToxicityEndpointResult = {
  active: boolean
  probability: number
  alerts: string[]
}

type ToxicityPrediction = {
  endpoints: Record<string, ToxicityEndpointResult>
  method: string
}

type CYPPrediction = {
  cyp3a4_substrate: boolean
  cyp2d6_substrate: boolean
  cyp2c9_substrate: boolean
  cyp_inhibitor_risk: boolean
  method: string
}

type ClinToxPrediction = {
  ct_tox: boolean
  probability: number
  structural_alerts: string[]
  fda_approval_likelihood: string
  method: string
}

type ADMETPredictions = {
  models_available: boolean
  model_error: string | null
  solubility: SolubilityPrediction | null
  lipophilicity: LipophilicityPrediction | null
  bbb_penetration: BBBPrediction | null
  cns_mpo: CNSMPOPrediction | null
  gi_absorption: GIAbsorptionPrediction | null
  toxicity: ToxicityPrediction | null
  cyp_metabolism: CYPPrediction | null
  clinical_toxicity: ClinToxPrediction | null
}

type QEDResult = {
  qed_score: number | null
  qed_class: string
}

type LeadLikeness = {
  violations: string[]
  num_violations: number
  category: string
}

type ScaffoldAnalysis = {
  murcko_scaffold: string | null
  generic_scaffold: string | null
  scaffold_rings: number
}

type MoleculeAnalysisResponse = {
  canonical_smiles: string
  molecule_image: string  // Base64 PNG image
  molecular_weight: number
  clogp: number
  tpsa: number
  hbd: number
  hba: number
  rotatable_bonds: number
  fraction_sp3: number
  ring_count: number
  aromatic_ring_count: number
  lipinski_rules: LipinskiRuleResult[]
  lipinski_summary: 'Pass' | 'Borderline' | 'Fail'
  alerts: Alert[]
  synthetic_accessibility_score: number
  synthetic_accessibility_class: 'Easy' | 'Moderate' | 'Difficult'
  flexibility_warnings: string[]
  complexity_warnings: string[]
  stereocenters: Stereocenter[]
  rs_comparisons: RSComparisonEntry[]
  stereochemistry_context: StereochemistryContext
  admet: ADMETPredictions
  qed: QEDResult
  lead_likeness: LeadLikeness
  scaffold: ScaffoldAnalysis
  fail_fast_score: number
  decision: 'Progress' | 'Optimize' | 'Kill'
  decision_rationale: string
  optimization_suggestions: OptimizationSuggestion[]
}

// Use environment variable for API URL, fallback to localhost for development
const API_BASE = import.meta.env.VITE_API_URL || 'https://drug-discovery-triage-1.onrender.com'

// Example molecules with defined stereochemistry
const EXAMPLE_MOLECULES = [
  { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O' },
  { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
  { name: '(S)-Ibuprofen', smiles: 'CC(C)Cc1ccc([C@H](C)C(=O)O)cc1' },
  { name: 'Penicillin V', smiles: 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)COC3=CC=CC=C3)C(=O)O)C' },
]

// Property thresholds for color coding
const THRESHOLDS = {
  mw: { good: 500, warning: 600 },
  clogp: { good: 5, warning: 6 },
  tpsa: { good: 140, warning: 160 },
  hbd: { good: 5, warning: 7 },
  hba: { good: 10, warning: 12 },
  rotatable: { good: 10, warning: 15 },
  fsp3: { good: 0.25, warning: 0.15 },
  sas: { good: 4, warning: 7 },
}

function App() {
  const [smiles, setSmiles] = useState('')
  const [selectedFlips, setSelectedFlips] = useState<number[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [data, setData] = useState<MoleculeAnalysisResponse | null>(null)
  const [flippedData, setFlippedData] = useState<MoleculeAnalysisResponse | null>(null)
  const [flippedLoading, setFlippedLoading] = useState(false)
  const [activeTab, setActiveTab] = useState(0)
  const [copied, setCopied] = useState(false)
  const [showComparison, setShowComparison] = useState(false)

  const handleAnalyze = async () => {
    setLoading(true)
    setError(null)
    setData(null)
    setFlippedData(null)
    setShowComparison(false)

    try {
      const res = await fetch(`${API_BASE}/analyze`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          smiles,
          flip_centers: selectedFlips,
        }),
      })

      if (!res.ok) {
        const body = await res.json().catch(() => null)
        throw new Error(body?.detail ?? 'Analysis failed')
      }

      const json = (await res.json()) as MoleculeAnalysisResponse
      setData(json)
      setSelectedFlips([])
    } catch (e) {
      setError(e instanceof Error ? e.message : 'Unexpected error')
    } finally {
      setLoading(false)
    }
  }

  // Fetch flipped molecule comparison
  const fetchFlippedComparison = async (atomIndex: number) => {
    if (!smiles) return

    setFlippedLoading(true)
    setShowComparison(true)

    try {
      const res = await fetch(`${API_BASE}/analyze`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          smiles,
          flip_centers: [atomIndex],
        }),
      })

      if (!res.ok) {
        throw new Error('Failed to analyze flipped molecule')
      }

      const json = (await res.json()) as MoleculeAnalysisResponse
      setFlippedData(json)
    } catch (e) {
      console.error('Error fetching flipped comparison:', e)
      setFlippedData(null)
    } finally {
      setFlippedLoading(false)
    }
  }

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && smiles.trim() && !loading) {
      handleAnalyze()
    }
  }

  const handleExampleClick = (exampleSmiles: string) => {
    setSmiles(exampleSmiles)
    setSelectedFlips([])
    setFlippedData(null)
    setShowComparison(false)
  }

  const copyToClipboard = async () => {
    if (data?.canonical_smiles) {
      await navigator.clipboard.writeText(data.canonical_smiles)
      setCopied(true)
      setTimeout(() => setCopied(false), 2000)
    }
  }

  const toggleFlip = (idx: number) => {
    setSelectedFlips((prev) =>
      prev.includes(idx) ? prev.filter((i) => i !== idx) : [...prev, idx]
    )
    // Fetch comparison when toggling
    fetchFlippedComparison(idx)
  }

  const closeComparison = () => {
    setShowComparison(false)
    setFlippedData(null)
  }

  // Color helpers
  const getDecisionColor = (decision: string) => {
    if (decision === 'Progress') return '#22c55e'
    if (decision === 'Optimize') return '#eab308'
    return '#ef4444'
  }

  const getDecisionBg = (decision: string) => {
    if (decision === 'Progress') return 'rgba(34, 197, 94, 0.15)'
    if (decision === 'Optimize') return 'rgba(234, 179, 8, 0.15)'
    return 'rgba(239, 68, 68, 0.15)'
  }

  const getScoreColor = (score: number) => {
    if (score >= 70) return '#22c55e'
    if (score >= 40) return '#eab308'
    return '#ef4444'
  }

  const getPropertyColor = (
    value: number,
    threshold: { good: number; warning: number },
    inverted = false
  ) => {
    if (inverted) {
      if (value >= threshold.good) return '#22c55e'
      if (value >= threshold.warning) return '#eab308'
      return '#ef4444'
    }
    if (value <= threshold.good) return '#22c55e'
    if (value <= threshold.warning) return '#eab308'
    return '#ef4444'
  }

  const aromaticPercent =
    data && data.ring_count > 0
      ? (100 * data.aromatic_ring_count) / data.ring_count
      : 0

  // Property card component
  const PropertyCard = ({
    label,
    value,
    unit = '',
    tooltip,
    color,
    progress,
    comparisonValue,
  }: {
    label: string
    value: string | number
    unit?: string
    tooltip?: string
    color?: string
    progress?: number
    comparisonValue?: string | number
  }) => (
    <Tooltip title={tooltip || ''} arrow placement="top">
      <Box
        sx={{
          p: 1.5,
          borderRadius: 2,
          backgroundColor: 'rgba(15, 23, 42, 0.6)',
          border: '1px solid rgba(148, 163, 184, 0.2)',
          transition: 'all 0.2s',
          '&:hover': {
            backgroundColor: 'rgba(15, 23, 42, 0.8)',
            borderColor: 'rgba(148, 163, 184, 0.4)',
          },
        }}
      >
        <Typography
          variant="caption"
          sx={{ color: 'rgb(148, 163, 184)', display: 'block', mb: 0.5 }}
        >
          {label}
        </Typography>
        <Stack direction="row" alignItems="baseline" spacing={1}>
          <Typography
            variant="h6"
            sx={{ color: color || '#e5e7eb', fontWeight: 600, lineHeight: 1.2 }}
          >
            {value}
            {unit && (
              <Typography
                component="span"
                variant="caption"
                sx={{ ml: 0.5, color: 'rgb(148, 163, 184)' }}
              >
                {unit}
              </Typography>
            )}
          </Typography>
          {comparisonValue !== undefined && (
            <Typography variant="caption" sx={{ color: '#38bdf8' }}>
              → {comparisonValue}{unit}
            </Typography>
          )}
        </Stack>
        {progress !== undefined && (
          <LinearProgress
            variant="determinate"
            value={Math.min(progress, 100)}
            sx={{
              mt: 1,
              height: 4,
              borderRadius: 2,
              backgroundColor: 'rgba(148, 163, 184, 0.2)',
              '& .MuiLinearProgress-bar': {
                backgroundColor: color || '#38bdf8',
                borderRadius: 2,
              },
            }}
          />
        )}
      </Box>
    </Tooltip>
  )

  return (
    <Box
      sx={{
        minHeight: '100vh',
        position: 'relative',
        overflow: 'hidden',
        '&::before': {
          content: '""',
          position: 'absolute',
          inset: '-20%',
          background:
            'radial-gradient(circle at 10% 0%, rgba(56,189,248,0.18) 0, transparent 50%), ' +
            'radial-gradient(circle at 90% 0%, rgba(168,85,247,0.2) 0, transparent 55%), ' +
            'radial-gradient(circle at 50% 110%, rgba(56,189,248,0.18) 0, transparent 55%)',
          filter: 'blur(40px)',
          opacity: 0.9,
          zIndex: 0,
        },
      }}
    >
      {/* Header */}
      <AppBar
        position="sticky"
        elevation={0}
        sx={{
          backgroundColor: 'rgba(15,23,42,0.85)',
          backdropFilter: 'blur(20px)',
          borderBottom: '1px solid rgba(148,163,184,0.15)',
        }}
      >
        <Toolbar sx={{ py: 1 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1.5, flexGrow: 1 }}>
            <Box sx={{ color: '#38bdf8' }}>
              <ScienceIcon />
            </Box>
            <Typography variant="h6" sx={{ fontWeight: 600 }}>
              Molecular Triage
            </Typography>
          </Box>
          <Chip
            size="small"
            variant="outlined"
            label="Chemistry-first · Transparent · Fail-fast"
            sx={{
              borderColor: 'rgba(148,163,184,0.4)',
              color: 'rgb(148,163,184)',
              fontSize: '0.7rem',
              display: { xs: 'none', sm: 'flex' },
            }}
          />
        </Toolbar>
      </AppBar>

      <Container sx={{ py: 3, position: 'relative', zIndex: 1 }} maxWidth="xl">
        {/* Input Section */}
        <Paper
          sx={{
            p: 2.5,
            mb: 3,
            borderRadius: 3,
            border: '1px solid rgba(148,163,184,0.25)',
            background:
              'linear-gradient(135deg, rgba(15,23,42,0.92), rgba(15,23,42,0.88)),' +
              'radial-gradient(circle at 0% 0%, rgba(56,189,248,0.12) 0, transparent 55%)',
            boxShadow: '0 8px 32px rgba(0,0,0,0.4)',
          }}
          elevation={0}
        >
          <Stack spacing={2}>
            <Stack
              direction={{ xs: 'column', md: 'row' }}
              spacing={2}
              alignItems={{ xs: 'stretch', md: 'center' }}
            >
              <TextField
                fullWidth
                label="SMILES"
                placeholder="Enter molecular SMILES string (e.g., CCO for ethanol)"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                onKeyDown={handleKeyDown}
                variant="outlined"
                size="small"
                InputProps={{
                  sx: {
                    borderRadius: 2,
                    backgroundColor: 'rgba(15,23,42,0.7)',
                  },
                }}
              />
              <Button
                variant="contained"
                onClick={handleAnalyze}
                disabled={!smiles.trim() || loading}
                sx={{
                  borderRadius: 2,
                  px: 4,
                  py: 1,
                  textTransform: 'none',
                  fontWeight: 600,
                  minWidth: 120,
                  boxShadow: '0 4px 16px rgba(56,189,248,0.3)',
                  '&:hover': {
                    boxShadow: '0 6px 20px rgba(56,189,248,0.4)',
                  },
                }}
              >
                {loading ? <CircularProgress size={20} color="inherit" /> : 'Analyze'}
              </Button>
            </Stack>

            {/* Example molecules */}
            <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
              <Typography
                variant="caption"
                sx={{ color: 'rgb(148,163,184)', alignSelf: 'center' }}
              >
                Try:
              </Typography>
              {EXAMPLE_MOLECULES.map((ex) => (
                <Chip
                  key={ex.name}
                  label={ex.name}
                  size="small"
                  variant="outlined"
                  onClick={() => handleExampleClick(ex.smiles)}
                  sx={{
                    borderColor: 'rgba(148,163,184,0.3)',
                    color: 'rgb(148,163,184)',
                    fontSize: '0.75rem',
                    cursor: 'pointer',
                    '&:hover': {
                      borderColor: '#38bdf8',
                      color: '#38bdf8',
                      backgroundColor: 'rgba(56,189,248,0.1)',
                    },
                  }}
                />
              ))}
            </Stack>

            {error && (
              <Typography color="error" variant="body2">
                {error}
              </Typography>
            )}
          </Stack>
        </Paper>

        {/* Results Dashboard */}
        {data && (
          <Stack spacing={2.5}>
            {/* Hero Row: Decision Score + Molecule Viewer */}
            <Grid container spacing={2.5}>
              {/* Decision Score Card */}
              <Grid size={{ xs: 12, md: 4 }}>
                <Paper
                  sx={{
                    p: 3,
                    height: '100%',
                    borderRadius: 3,
                    border: `2px solid ${getDecisionColor(data.decision)}40`,
                    background: `linear-gradient(135deg, ${getDecisionBg(data.decision)}, rgba(15,23,42,0.95))`,
                    boxShadow: `0 8px 32px ${getDecisionColor(data.decision)}20`,
                    display: 'flex',
                    flexDirection: 'column',
                  }}
                  elevation={0}
                >
                  <Stack direction="row" alignItems="center" spacing={1} sx={{ mb: 2 }}>
                    <Box sx={{ color: getDecisionColor(data.decision) }}>
                      <SpeedIcon />
                    </Box>
                    <Typography variant="subtitle2" sx={{ color: 'rgb(148,163,184)' }}>
                      Triage Decision
                    </Typography>
                  </Stack>

                  <Box sx={{ textAlign: 'center', flex: 1, display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
                    <Typography
                      variant="h2"
                      sx={{
                        fontWeight: 700,
                        color: getScoreColor(data.fail_fast_score),
                        lineHeight: 1,
                      }}
                    >
                      {data.fail_fast_score.toFixed(0)}
                    </Typography>
                    <Typography variant="caption" sx={{ color: 'rgb(148,163,184)', mb: 2 }}>
                      Drug-Likeness Score
                    </Typography>

                    <Chip
                      label={data.decision}
                      sx={{
                        alignSelf: 'center',
                        backgroundColor: getDecisionColor(data.decision),
                        color: '#000',
                        fontWeight: 700,
                        fontSize: '1rem',
                        py: 2.5,
                        px: 2,
                      }}
                    />
                  </Box>

                  <Typography
                    variant="body2"
                    sx={{ mt: 2, color: 'rgb(148,163,184)', textAlign: 'center', fontSize: '0.8rem' }}
                  >
                    {data.decision_rationale}
                  </Typography>

                  {/* Lipinski Summary */}
                  <Box sx={{ mt: 2, pt: 2, borderTop: '1px solid rgba(148,163,184,0.2)' }}>
                    <Stack direction="row" justifyContent="space-between" alignItems="center">
                      <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>
                        Lipinski RO5
                      </Typography>
                      <Chip
                        size="small"
                        label={data.lipinski_summary}
                        sx={{
                          backgroundColor:
                            data.lipinski_summary === 'Pass'
                              ? 'rgba(34, 197, 94, 0.2)'
                              : data.lipinski_summary === 'Borderline'
                                ? 'rgba(234, 179, 8, 0.2)'
                                : 'rgba(239, 68, 68, 0.2)',
                          color:
                            data.lipinski_summary === 'Pass'
                              ? '#22c55e'
                              : data.lipinski_summary === 'Borderline'
                                ? '#eab308'
                                : '#ef4444',
                          fontWeight: 600,
                        }}
                      />
                    </Stack>
                    <Stack direction="row" spacing={0.5} sx={{ mt: 1 }}>
                      {data.lipinski_rules.map((rule) => (
                        <Tooltip key={rule.name} title={`${rule.name}: ${rule.detail}`} arrow>
                          <Box
                            sx={{
                              flex: 1,
                              height: 6,
                              borderRadius: 1,
                              backgroundColor: rule.passed ? '#22c55e' : '#ef4444',
                            }}
                          />
                        </Tooltip>
                      ))}
                    </Stack>
                  </Box>
                </Paper>
              </Grid>

              {/* Molecule Viewer */}
              <Grid size={{ xs: 12, md: 8 }}>
                <Paper
                  sx={{
                    p: 2.5,
                    height: '100%',
                    minHeight: 420,
                    borderRadius: 3,
                    border: '1px solid rgba(148,163,184,0.25)',
                    background: 'radial-gradient(circle at 50% 0%, rgba(56,189,248,0.08) 0, transparent 60%), rgba(15,23,42,0.95)',
                    boxShadow: '0 8px 32px rgba(0,0,0,0.3)',
                  }}
                  elevation={0}
                >
                  <Stack direction="row" justifyContent="space-between" alignItems="center">
                    <Stack direction="row" alignItems="center" spacing={1}>
                      <Box sx={{ color: '#38bdf8' }}>
                        <MoleculeIcon />
                      </Box>
                      <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
                        {showComparison ? '2D Structure Comparison' : '2D Structure'}
                      </Typography>
                      {showComparison && (
                        <Chip
                          size="small"
                          label="Comparing R↔S"
                          icon={<CompareIcon />}
                          onDelete={closeComparison}
                          sx={{
                            ml: 1,
                            backgroundColor: 'rgba(56,189,248,0.15)',
                            color: '#38bdf8',
                            '& .MuiChip-deleteIcon': { color: '#38bdf8' },
                          }}
                        />
                      )}
                    </Stack>
                    <Stack direction="row" alignItems="center" spacing={1}>
                      <Tooltip title={copied ? 'Copied!' : 'Copy canonical SMILES'}>
                        <IconButton
                          size="small"
                          onClick={copyToClipboard}
                          sx={{ color: copied ? '#22c55e' : 'rgb(148,163,184)' }}
                        >
                          <CopyIcon />
                        </IconButton>
                      </Tooltip>
                      <Typography
                        variant="caption"
                        sx={{
                          color: 'rgb(148,163,184)',
                          fontFamily: 'monospace',
                          maxWidth: 200,
                          overflow: 'hidden',
                          textOverflow: 'ellipsis',
                          whiteSpace: 'nowrap',
                        }}
                      >
                        {data.canonical_smiles}
                      </Typography>
                    </Stack>
                  </Stack>

                  {/* Molecule Images */}
                  <Grid container spacing={2} sx={{ mt: 1 }}>
                    <Grid size={{ xs: 12, md: showComparison ? 6 : 12 }}>
                      <Box
                        sx={{
                          display: 'flex',
                          justifyContent: 'center',
                          alignItems: 'center',
                          flexDirection: 'column',
                          borderRadius: 2,
                          border: '1px solid rgba(148,163,184,0.15)',
                          backgroundColor: '#0f172a',
                          minHeight: 320,
                          overflow: 'hidden',
                          p: 1,
                        }}
                      >
                        <img
                          src={data.molecule_image}
                          alt="Molecule structure"
                          style={{
                            maxWidth: '100%',
                            maxHeight: 300,
                            objectFit: 'contain',
                          }}
                        />
                        <Typography variant="caption" sx={{ color: 'rgb(148,163,184)', mt: 1 }}>
                          Original
                        </Typography>
                      </Box>
                    </Grid>

                    {showComparison && (
                      <Grid size={{ xs: 12, md: 6 }}>
                        <Box
                          sx={{
                            display: 'flex',
                            justifyContent: 'center',
                            alignItems: 'center',
                            flexDirection: 'column',
                            borderRadius: 2,
                            border: '1px solid rgba(56,189,248,0.3)',
                            backgroundColor: '#0f172a',
                            minHeight: 320,
                            overflow: 'hidden',
                            p: 1,
                          }}
                        >
                          {flippedLoading ? (
                            <CircularProgress size={40} sx={{ color: '#38bdf8' }} />
                          ) : flippedData?.molecule_image ? (
                            <>
                              <img
                                src={flippedData.molecule_image}
                                alt="Flipped molecule structure"
                                style={{
                                  maxWidth: '100%',
                                  maxHeight: 300,
                                  objectFit: 'contain',
                                }}
                              />
                              <Typography variant="caption" sx={{ color: '#38bdf8', mt: 1 }}>
                                Flipped (R↔S)
                              </Typography>
                            </>
                          ) : (
                            <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>
                              Loading flipped structure...
                            </Typography>
                          )}
                        </Box>
                      </Grid>
                    )}
                  </Grid>

                  {/* Atom color legend */}
                  <Stack direction="row" spacing={2} sx={{ mt: 1.5, justifyContent: 'center', flexWrap: 'wrap' }}>
                    {[
                      { atom: 'C', color: '#e5e5e5' },
                      { atom: 'O', color: '#f04545' },
                      { atom: 'N', color: '#3b82f6' },
                      { atom: 'S', color: '#eab308' },
                      { atom: 'F/Cl', color: '#22c55e' },
                      { atom: 'Br/I', color: '#a855f7' },
                    ].map(({ atom, color }) => (
                      <Stack key={atom} direction="row" spacing={0.5} alignItems="center">
                        <Box
                          sx={{
                            width: 10,
                            height: 10,
                            borderRadius: '50%',
                            backgroundColor: color,
                          }}
                        />
                        <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>
                          {atom}
                        </Typography>
                      </Stack>
                    ))}
                  </Stack>
                </Paper>
              </Grid>
            </Grid>

            {/* Property Comparison (when R/S is toggled) */}
            {showComparison && flippedData && (
              <Paper
                sx={{
                  p: 2.5,
                  borderRadius: 3,
                  border: '1px solid rgba(56,189,248,0.3)',
                  backgroundColor: 'rgba(15,23,42,0.92)',
                  boxShadow: '0 8px 32px rgba(56,189,248,0.1)',
                }}
                elevation={0}
              >
                <Stack direction="row" alignItems="center" spacing={1} sx={{ mb: 2 }}>
                  <CompareIcon />
                  <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
                    Property Comparison (Original → Flipped)
                  </Typography>
                </Stack>
                <Grid container spacing={1.5}>
                  <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                    <PropertyCard
                      label="Mol. Weight"
                      value={data.molecular_weight.toFixed(1)}
                      comparisonValue={flippedData.molecular_weight.toFixed(1)}
                      unit="Da"
                      tooltip="Molecular weight"
                    />
                  </Grid>
                  <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                    <PropertyCard
                      label="cLogP"
                      value={data.clogp.toFixed(2)}
                      comparisonValue={flippedData.clogp.toFixed(2)}
                      tooltip="Calculated partition coefficient"
                    />
                  </Grid>
                  <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                    <PropertyCard
                      label="TPSA"
                      value={data.tpsa.toFixed(1)}
                      comparisonValue={flippedData.tpsa.toFixed(1)}
                      unit="Ų"
                      tooltip="Topological polar surface area"
                    />
                  </Grid>
                  <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                    <PropertyCard
                      label="HBD / HBA"
                      value={`${data.hbd}/${data.hba}`}
                      comparisonValue={`${flippedData.hbd}/${flippedData.hba}`}
                      tooltip="Hydrogen bond donors/acceptors"
                    />
                  </Grid>
                  <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                    <PropertyCard
                      label="Drug Score"
                      value={data.fail_fast_score.toFixed(0)}
                      comparisonValue={flippedData.fail_fast_score.toFixed(0)}
                      tooltip="Fail-fast drug-likeness score"
                      color={getScoreColor(data.fail_fast_score)}
                    />
                  </Grid>
                  <Grid size={{ xs: 6, sm: 4, md: 2 }}>
                    <Box
                      sx={{
                        p: 1.5,
                        borderRadius: 2,
                        backgroundColor: 'rgba(15, 23, 42, 0.6)',
                        border: '1px solid rgba(148, 163, 184, 0.2)',
                      }}
                    >
                      <Typography variant="caption" sx={{ color: 'rgb(148, 163, 184)', display: 'block', mb: 0.5 }}>
                        Decision
                      </Typography>
                      <Stack direction="row" alignItems="center" spacing={1}>
                        <Chip
                          size="small"
                          label={data.decision}
                          sx={{ backgroundColor: getDecisionColor(data.decision), color: '#000', fontWeight: 600 }}
                        />
                        <Typography sx={{ color: 'rgb(148,163,184)' }}>→</Typography>
                        <Chip
                          size="small"
                          label={flippedData.decision}
                          sx={{ backgroundColor: getDecisionColor(flippedData.decision), color: '#000', fontWeight: 600 }}
                        />
                      </Stack>
                    </Box>
                  </Grid>
                </Grid>
              </Paper>
            )}

            {/* Properties Grid */}
            <Paper
              sx={{
                p: 2.5,
                borderRadius: 3,
                border: '1px solid rgba(148,163,184,0.25)',
                backgroundColor: 'rgba(15,23,42,0.92)',
                boxShadow: '0 8px 32px rgba(0,0,0,0.3)',
              }}
              elevation={0}
            >
              <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2 }}>
                Physicochemical Properties
              </Typography>
              <Grid container spacing={1.5}>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="Mol. Weight"
                    value={data.molecular_weight.toFixed(1)}
                    unit="Da"
                    tooltip="Molecular weight. Ideal < 500 Da for oral drugs"
                    color={getPropertyColor(data.molecular_weight, THRESHOLDS.mw)}
                    progress={(data.molecular_weight / 600) * 100}
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="cLogP"
                    value={data.clogp.toFixed(2)}
                    tooltip="Calculated partition coefficient. Ideal < 5 for drug-likeness"
                    color={getPropertyColor(data.clogp, THRESHOLDS.clogp)}
                    progress={(Math.max(0, data.clogp + 2) / 8) * 100}
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="TPSA"
                    value={data.tpsa.toFixed(1)}
                    unit="Ų"
                    tooltip="Topological polar surface area. Ideal < 140 Ų for oral bioavailability"
                    color={getPropertyColor(data.tpsa, THRESHOLDS.tpsa)}
                    progress={(data.tpsa / 200) * 100}
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="HBD / HBA"
                    value={`${data.hbd} / ${data.hba}`}
                    tooltip="Hydrogen bond donors / acceptors. Ideal: HBD ≤ 5, HBA ≤ 10"
                    color={
                      data.hbd <= 5 && data.hba <= 10
                        ? '#22c55e'
                        : data.hbd <= 7 && data.hba <= 12
                          ? '#eab308'
                          : '#ef4444'
                    }
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="Rotatable"
                    value={data.rotatable_bonds}
                    tooltip="Rotatable bonds. Ideal ≤ 10 for oral bioavailability"
                    color={getPropertyColor(data.rotatable_bonds, THRESHOLDS.rotatable)}
                    progress={(data.rotatable_bonds / 15) * 100}
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="Fsp³"
                    value={data.fraction_sp3.toFixed(2)}
                    tooltip="Fraction of sp³ carbons. Higher values (> 0.25) often indicate better developability"
                    color={getPropertyColor(data.fraction_sp3, THRESHOLDS.fsp3, true)}
                    progress={data.fraction_sp3 * 100}
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="Rings"
                    value={data.ring_count}
                    tooltip="Total ring count in the molecule"
                  />
                </Grid>
                <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                  <PropertyCard
                    label="Aromatic %"
                    value={`${aromaticPercent.toFixed(0)}%`}
                    tooltip="Percentage of aromatic rings. High aromatic content may indicate poor solubility"
                    color={aromaticPercent > 80 ? '#ef4444' : aromaticPercent > 60 ? '#eab308' : '#22c55e'}
                    progress={aromaticPercent}
                  />
                </Grid>

                {/* New properties */}
                {data.qed && data.qed.qed_score !== null && (
                  <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                    <PropertyCard
                      label="QED"
                      value={data.qed.qed_score.toFixed(2)}
                      tooltip={`Quantitative Estimate of Drug-likeness (Bickerton 2012). ${data.qed.qed_class}`}
                      color={data.qed.qed_score >= 0.67 ? '#22c55e' : data.qed.qed_score >= 0.49 ? '#eab308' : '#ef4444'}
                      progress={data.qed.qed_score * 100}
                    />
                  </Grid>
                )}

                {data.lead_likeness && (
                  <Grid size={{ xs: 6, sm: 4, md: 3, lg: 1.5 }}>
                    <PropertyCard
                      label="Lead-like"
                      value={data.lead_likeness.category.replace('lead-like', 'LL').replace('Lead-like', 'Yes').replace('Borderline ', 'Border').replace('Not ', 'No')}
                      tooltip={`Rule of 3 (Congreve 2003) for fragment/lead screening:\n• MW 150-350 Da\n• cLogP ≤ 3\n• HBD ≤ 3, HBA ≤ 3\n• Rotatable bonds ≤ 3\n\n${data.lead_likeness.num_violations === 0 ? '✓ All criteria met!' : `Violations (${data.lead_likeness.num_violations}):\n${data.lead_likeness.violations.map(v => '• ' + v).join('\n')}`}`}
                      color={data.lead_likeness.num_violations === 0 ? '#22c55e' : data.lead_likeness.num_violations <= 2 ? '#eab308' : '#ef4444'}
                    />
                  </Grid>
                )}
              </Grid>
            </Paper>

            {/* Tabbed Analysis Section */}
            <Paper
              sx={{
                borderRadius: 3,
                border: '1px solid rgba(148,163,184,0.25)',
                backgroundColor: 'rgba(15,23,42,0.92)',
                boxShadow: '0 8px 32px rgba(0,0,0,0.3)',
                overflow: 'hidden',
              }}
              elevation={0}
            >
              <Tabs
                value={activeTab}
                onChange={(_, v) => setActiveTab(v)}
                variant="scrollable"
                scrollButtons="auto"
                sx={{
                  borderBottom: '1px solid rgba(148,163,184,0.2)',
                  '& .MuiTab-root': {
                    textTransform: 'none',
                    fontWeight: 500,
                    minHeight: 48,
                  },
                }}
              >
                <Tab icon={<ADMETIcon />} iconPosition="start" label="ADMET" />
                <Tab icon={<AlertIcon />} iconPosition="start" label={`Alerts (${data.alerts.length})`} />
                <Tab icon={<ChiralIcon />} iconPosition="start" label={`Stereochemistry (${data.stereocenters.length})`} />
                <Tab icon={<LightbulbIcon />} iconPosition="start" label={`Suggestions (${data.optimization_suggestions.length})`} />
              </Tabs>

              <Box sx={{ p: 2.5 }}>
                {/* ADMET Tab */}
                {activeTab === 0 && data.admet && (
                  <Stack spacing={2.5}>
                    {/* Row 1: Absorption & Distribution */}
                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1.5, color: '#38bdf8' }}>
                        Absorption & Distribution
                      </Typography>
                      <Grid container spacing={2}>
                        {/* Solubility */}
                        <Grid size={{ xs: 12, sm: 6, md: 3 }}>
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(15,23,42,0.6)', border: '1px solid rgba(148,163,184,0.2)' }}>
                            <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>Aqueous Solubility</Typography>
                            {data.admet.solubility ? (
                              <>
                                <Typography variant="h6" sx={{ color: data.admet.solubility.solubility_class === 'High' ? '#22c55e' : data.admet.solubility.solubility_class === 'Moderate' ? '#eab308' : '#ef4444' }}>
                                  {data.admet.solubility.solubility_class}
                                </Typography>
                                <Typography variant="caption" sx={{ color: 'rgb(100,116,139)' }}>
                                  Log S: {data.admet.solubility.log_s}
                                </Typography>
                              </>
                            ) : (
                              <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>N/A</Typography>
                            )}
                          </Box>
                        </Grid>

                        {/* GI Absorption */}
                        <Grid size={{ xs: 12, sm: 6, md: 3 }}>
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(15,23,42,0.6)', border: '1px solid rgba(148,163,184,0.2)' }}>
                            <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>GI Absorption</Typography>
                            {data.admet.gi_absorption ? (
                              <>
                                <Typography variant="h6" sx={{ color: data.admet.gi_absorption.absorption === 'High' ? '#22c55e' : data.admet.gi_absorption.absorption === 'Moderate' ? '#eab308' : '#ef4444' }}>
                                  {data.admet.gi_absorption.absorption}
                                </Typography>
                                <Typography variant="caption" sx={{ color: 'rgb(100,116,139)' }}>
                                  Bioavailability: {(data.admet.gi_absorption.bioavailability_score * 100).toFixed(0)}%
                                </Typography>
                              </>
                            ) : (
                              <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>N/A</Typography>
                            )}
                          </Box>
                        </Grid>

                        {/* BBB Penetration */}
                        <Grid size={{ xs: 12, sm: 6, md: 3 }}>
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(15,23,42,0.6)', border: '1px solid rgba(148,163,184,0.2)' }}>
                            <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>BBB Penetration</Typography>
                            {data.admet.bbb_penetration ? (
                              <>
                                <Stack direction="row" alignItems="center" spacing={1}>
                                  <Typography variant="h6" sx={{ color: data.admet.bbb_penetration.penetrates ? '#a855f7' : '#eab308' }}>
                                    {data.admet.bbb_penetration.penetrates ? 'BBB+' : 'BBB-'}
                                  </Typography>
                                  <Chip size="small" label={data.admet.bbb_penetration.confidence} sx={{ backgroundColor: 'rgba(148,163,184,0.2)', color: 'rgb(148,163,184)', fontSize: '0.65rem' }} />
                                </Stack>
                                <Typography variant="caption" sx={{ color: 'rgb(100,116,139)' }}>
                                  Probability: {(data.admet.bbb_penetration.probability * 100).toFixed(0)}%
                                </Typography>
                              </>
                            ) : (
                              <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>N/A</Typography>
                            )}
                          </Box>
                        </Grid>

                        {/* CNS MPO */}
                        <Grid size={{ xs: 12, sm: 6, md: 3 }}>
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(15,23,42,0.6)', border: '1px solid rgba(148,163,184,0.2)' }}>
                            <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>CNS MPO Score</Typography>
                            {data.admet.cns_mpo ? (
                              <>
                                <Stack direction="row" alignItems="baseline" spacing={0.5}>
                                  <Typography variant="h6" sx={{ color: data.admet.cns_mpo.cns_class === 'High' ? '#22c55e' : data.admet.cns_mpo.cns_class === 'Medium' ? '#eab308' : '#ef4444' }}>
                                    {data.admet.cns_mpo.score.toFixed(1)}
                                  </Typography>
                                  <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>/6</Typography>
                                </Stack>
                                <Chip size="small" label={`${data.admet.cns_mpo.cns_class} CNS`} sx={{ backgroundColor: data.admet.cns_mpo.cns_class === 'High' ? 'rgba(34,197,94,0.2)' : data.admet.cns_mpo.cns_class === 'Medium' ? 'rgba(234,179,8,0.2)' : 'rgba(239,68,68,0.2)', color: data.admet.cns_mpo.cns_class === 'High' ? '#22c55e' : data.admet.cns_mpo.cns_class === 'Medium' ? '#eab308' : '#ef4444', fontSize: '0.7rem', mt: 0.5 }} />
                              </>
                            ) : (
                              <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>N/A</Typography>
                            )}
                          </Box>
                        </Grid>
                      </Grid>
                    </Box>

                    {/* Row 2: Metabolism (CYP450) */}
                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1.5, color: '#a855f7' }}>
                        Metabolism (CYP450)
                      </Typography>
                      {data.admet.cyp_metabolism ? (
                        <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap>
                          <Chip
                            label="CYP3A4 Substrate"
                            size="small"
                            sx={{
                              backgroundColor: data.admet.cyp_metabolism.cyp3a4_substrate ? 'rgba(168,85,247,0.2)' : 'rgba(148,163,184,0.1)',
                              color: data.admet.cyp_metabolism.cyp3a4_substrate ? '#a855f7' : 'rgb(148,163,184)',
                              border: data.admet.cyp_metabolism.cyp3a4_substrate ? '1px solid rgba(168,85,247,0.4)' : '1px solid transparent',
                            }}
                          />
                          <Chip
                            label="CYP2D6 Substrate"
                            size="small"
                            sx={{
                              backgroundColor: data.admet.cyp_metabolism.cyp2d6_substrate ? 'rgba(56,189,248,0.2)' : 'rgba(148,163,184,0.1)',
                              color: data.admet.cyp_metabolism.cyp2d6_substrate ? '#38bdf8' : 'rgb(148,163,184)',
                              border: data.admet.cyp_metabolism.cyp2d6_substrate ? '1px solid rgba(56,189,248,0.4)' : '1px solid transparent',
                            }}
                          />
                          <Chip
                            label="CYP2C9 Substrate"
                            size="small"
                            sx={{
                              backgroundColor: data.admet.cyp_metabolism.cyp2c9_substrate ? 'rgba(234,179,8,0.2)' : 'rgba(148,163,184,0.1)',
                              color: data.admet.cyp_metabolism.cyp2c9_substrate ? '#eab308' : 'rgb(148,163,184)',
                              border: data.admet.cyp_metabolism.cyp2c9_substrate ? '1px solid rgba(234,179,8,0.4)' : '1px solid transparent',
                            }}
                          />
                          {data.admet.cyp_metabolism.cyp_inhibitor_risk && (
                            <Chip
                              label="CYP Inhibitor Risk"
                              size="small"
                              sx={{
                                backgroundColor: 'rgba(239,68,68,0.2)',
                                color: '#ef4444',
                                border: '1px solid rgba(239,68,68,0.4)',
                              }}
                            />
                          )}
                        </Stack>
                      ) : (
                        <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>Metabolism data not available</Typography>
                      )}
                    </Box>

                    {/* Row 3: Toxicity */}
                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1.5, color: '#ef4444' }}>
                        Toxicity Predictions
                      </Typography>
                      <Grid container spacing={2}>
                        {/* Clinical Toxicity */}
                        <Grid size={{ xs: 12, md: 6 }}>
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(15,23,42,0.6)', border: '1px solid rgba(148,163,184,0.2)' }}>
                            <Stack direction="row" justifyContent="space-between" alignItems="center" sx={{ mb: 1 }}>
                              <Typography variant="body2" sx={{ fontWeight: 600 }}>Clinical Toxicity Risk</Typography>
                              {data.admet.clinical_toxicity && (
                                <Chip
                                  size="small"
                                  label={data.admet.clinical_toxicity.ct_tox ? 'High Risk' : 'Lower Risk'}
                                  sx={{
                                    backgroundColor: data.admet.clinical_toxicity.ct_tox ? 'rgba(239,68,68,0.2)' : 'rgba(34,197,94,0.2)',
                                    color: data.admet.clinical_toxicity.ct_tox ? '#ef4444' : '#22c55e',
                                  }}
                                />
                              )}
                            </Stack>
                            {data.admet.clinical_toxicity ? (
                              <>
                                <LinearProgress
                                  variant="determinate"
                                  value={data.admet.clinical_toxicity.probability * 100}
                                  sx={{
                                    height: 6,
                                    borderRadius: 3,
                                    backgroundColor: 'rgba(148,163,184,0.2)',
                                    '& .MuiLinearProgress-bar': {
                                      backgroundColor: data.admet.clinical_toxicity.probability > 0.5 ? '#ef4444' : data.admet.clinical_toxicity.probability > 0.3 ? '#eab308' : '#22c55e',
                                    },
                                  }}
                                />
                                <Stack direction="row" justifyContent="space-between" sx={{ mt: 1 }}>
                                  <Typography variant="caption" sx={{ color: 'rgb(100,116,139)' }}>
                                    FDA Approval: {data.admet.clinical_toxicity.fda_approval_likelihood}
                                  </Typography>
                                  <Typography variant="caption" sx={{ color: 'rgb(100,116,139)' }}>
                                    {(data.admet.clinical_toxicity.probability * 100).toFixed(0)}% risk
                                  </Typography>
                                </Stack>
                                {data.admet.clinical_toxicity.structural_alerts.length > 0 && (
                                  <Box sx={{ mt: 1 }}>
                                    <Typography variant="caption" sx={{ color: '#ef4444' }}>
                                      Alerts: {data.admet.clinical_toxicity.structural_alerts.join(', ')}
                                    </Typography>
                                  </Box>
                                )}
                              </>
                            ) : (
                              <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>Data not available</Typography>
                            )}
                          </Box>
                        </Grid>

                        {/* Toxicity Endpoints */}
                        <Grid size={{ xs: 12, md: 6 }}>
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(15,23,42,0.6)', border: '1px solid rgba(148,163,184,0.2)' }}>
                            <Typography variant="body2" sx={{ fontWeight: 600, mb: 1 }}>Toxicity Endpoints (Tox21)</Typography>
                            {data.admet.toxicity?.endpoints ? (
                              <Stack spacing={0.5}>
                                {Object.entries(data.admet.toxicity.endpoints).map(([key, endpoint]) => (
                                  <Stack key={key} direction="row" justifyContent="space-between" alignItems="center">
                                    <Typography variant="caption" sx={{ color: 'rgb(148,163,184)', textTransform: 'uppercase', fontSize: '0.65rem' }}>
                                      {key.replace(/_/g, ' ')}
                                    </Typography>
                                    <Chip
                                      size="small"
                                      label={endpoint.active ? 'Active' : 'Inactive'}
                                      sx={{
                                        height: 18,
                                        fontSize: '0.6rem',
                                        backgroundColor: endpoint.active ? 'rgba(239,68,68,0.2)' : 'rgba(34,197,94,0.2)',
                                        color: endpoint.active ? '#ef4444' : '#22c55e',
                                      }}
                                    />
                                  </Stack>
                                ))}
                              </Stack>
                            ) : (
                              <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>Data not available</Typography>
                            )}
                          </Box>
                        </Grid>
                      </Grid>
                    </Box>

                    {/* Method Attribution */}
                    <Box sx={{ pt: 1, borderTop: '1px solid rgba(148,163,184,0.1)' }}>
                      <Typography variant="caption" sx={{ color: 'rgb(100,116,139)' }}>
                        <strong>Rule-based predictions:</strong> ESOL (Delaney 2004) for solubility, Lipinski/Veber rules for absorption,
                        CNS MPO (Pfizer), SMARTS-based structural alerts for toxicity. Fast, transparent, and scientifically validated.
                      </Typography>
                    </Box>
                  </Stack>
                )}

                {/* Alerts Tab */}
                {activeTab === 1 && (
                  <Stack spacing={2}>
                    <Box>
                      <Stack direction="row" justifyContent="space-between" alignItems="center" sx={{ mb: 1 }}>
                        <Typography variant="subtitle2">Synthetic Accessibility</Typography>
                        <Chip
                          size="small"
                          label={`${data.synthetic_accessibility_score.toFixed(1)} - ${data.synthetic_accessibility_class}`}
                          sx={{
                            backgroundColor:
                              data.synthetic_accessibility_class === 'Easy'
                                ? 'rgba(34, 197, 94, 0.2)'
                                : data.synthetic_accessibility_class === 'Moderate'
                                  ? 'rgba(234, 179, 8, 0.2)'
                                  : 'rgba(239, 68, 68, 0.2)',
                            color:
                              data.synthetic_accessibility_class === 'Easy'
                                ? '#22c55e'
                                : data.synthetic_accessibility_class === 'Moderate'
                                  ? '#eab308'
                                  : '#ef4444',
                          }}
                        />
                      </Stack>
                      <LinearProgress
                        variant="determinate"
                        value={(data.synthetic_accessibility_score / 10) * 100}
                        sx={{
                          height: 8,
                          borderRadius: 4,
                          backgroundColor: 'rgba(148,163,184,0.2)',
                          '& .MuiLinearProgress-bar': {
                            backgroundColor: getPropertyColor(data.synthetic_accessibility_score, THRESHOLDS.sas),
                          },
                        }}
                      />
                      <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>
                        Scale: 1 (easy) to 10 (difficult)
                      </Typography>
                    </Box>

                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1 }}>Medicinal Chemistry Alerts</Typography>
                      {data.alerts.length === 0 ? (
                        <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(34, 197, 94, 0.1)', border: '1px solid rgba(34, 197, 94, 0.3)' }}>
                          <Typography variant="body2" sx={{ color: '#22c55e' }}>No structural alerts detected</Typography>
                        </Box>
                      ) : (
                        <Stack spacing={1}>
                          {data.alerts.map((alert) => (
                            <Box key={alert.name} sx={{ p: 1.5, borderRadius: 2, backgroundColor: 'rgba(239, 68, 68, 0.1)', border: '1px solid rgba(239, 68, 68, 0.3)' }}>
                              <Typography variant="body2" sx={{ fontWeight: 600, color: '#ef4444' }}>{alert.name}</Typography>
                              <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>{alert.description}</Typography>
                            </Box>
                          ))}
                        </Stack>
                      )}
                    </Box>

                    <Grid container spacing={2}>
                      <Grid size={{ xs: 12, md: 6 }}>
                        <Typography variant="subtitle2" sx={{ mb: 1 }}>Flexibility Warnings</Typography>
                        {data.flexibility_warnings.length === 0 ? (
                          <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>No flexibility concerns</Typography>
                        ) : (
                          <Stack spacing={0.5}>
                            {data.flexibility_warnings.map((w, i) => (
                              <Typography key={i} variant="body2" sx={{ color: '#eab308' }}>• {w}</Typography>
                            ))}
                          </Stack>
                        )}
                      </Grid>
                      <Grid size={{ xs: 12, md: 6 }}>
                        <Typography variant="subtitle2" sx={{ mb: 1 }}>Complexity Warnings</Typography>
                        {data.complexity_warnings.length === 0 ? (
                          <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>No complexity concerns</Typography>
                        ) : (
                          <Stack spacing={0.5}>
                            {data.complexity_warnings.map((w, i) => (
                              <Typography key={i} variant="body2" sx={{ color: '#eab308' }}>• {w}</Typography>
                            ))}
                          </Stack>
                        )}
                      </Grid>
                    </Grid>
                  </Stack>
                )}

                {/* Stereochemistry Tab */}
                {activeTab === 2 && (
                  <Grid container spacing={2.5}>
                    {/* Left Column: Chiral Centers */}
                    <Grid size={{ xs: 12, md: 5 }}>
                      <Typography variant="subtitle2" sx={{ mb: 1.5 }}>Chiral Centers</Typography>
                      {data.stereocenters.length === 0 ? (
                        <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(148,163,184,0.1)' }}>
                          <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>No chiral centers detected - achiral molecule</Typography>
                        </Box>
                      ) : (
                        <Stack spacing={1}>
                          {/* Stereoisomer count summary */}
                          <Box sx={{ p: 1.5, borderRadius: 2, backgroundColor: 'rgba(56,189,248,0.1)', border: '1px solid rgba(56,189,248,0.3)', mb: 1 }}>
                            <Stack direction="row" justifyContent="space-between" alignItems="center">
                              <Typography variant="body2" sx={{ color: '#38bdf8' }}>
                                {data.stereochemistry_context.num_stereocenters} chiral center{data.stereochemistry_context.num_stereocenters > 1 ? 's' : ''}
                              </Typography>
                              <Chip
                                size="small"
                                label={`Up to ${data.stereochemistry_context.max_stereoisomers} stereoisomers`}
                                sx={{ backgroundColor: 'rgba(56,189,248,0.2)', color: '#38bdf8', fontWeight: 500 }}
                              />
                            </Stack>
                          </Box>

                          {data.stereocenters.map((c) => (
                            <Stack
                              key={c.atom_index}
                              direction="row"
                              spacing={2}
                              alignItems="center"
                              sx={{
                                p: 1.5,
                                borderRadius: 2,
                                backgroundColor: selectedFlips.includes(c.atom_index) ? 'rgba(56,189,248,0.15)' : 'rgba(148,163,184,0.1)',
                                border: selectedFlips.includes(c.atom_index) ? '1px solid rgba(56,189,248,0.4)' : '1px solid transparent',
                              }}
                            >
                              <Box sx={{ flex: 1 }}>
                                <Typography variant="body2">Atom #{c.atom_index}</Typography>
                                <Chip
                                  size="small"
                                  label={c.cip_label}
                                  sx={{
                                    mt: 0.5,
                                    backgroundColor:
                                      c.cip_label === 'R' ? 'rgba(56,189,248,0.2)' : c.cip_label === 'S' ? 'rgba(168,85,247,0.2)' : 'rgba(148,163,184,0.2)',
                                    color: c.cip_label === 'R' ? '#38bdf8' : c.cip_label === 'S' ? '#a855f7' : 'rgb(148,163,184)',
                                    fontWeight: 600,
                                  }}
                                />
                              </Box>
                              <Button
                                size="small"
                                variant={selectedFlips.includes(c.atom_index) ? 'contained' : 'outlined'}
                                onClick={() => toggleFlip(c.atom_index)}
                                sx={{ textTransform: 'none' }}
                              >
                                {selectedFlips.includes(c.atom_index) ? 'Comparing...' : 'Flip R↔S'}
                              </Button>
                            </Stack>
                          ))}
                          <Typography variant="caption" sx={{ color: 'rgb(148,163,184)', mt: 1 }}>
                            Click "Flip R↔S" to visualize the enantiomer structure above
                          </Typography>
                        </Stack>
                      )}
                    </Grid>

                    {/* Right Column: Educational Content */}
                    <Grid size={{ xs: 12, md: 7 }}>
                      <Typography variant="subtitle2" sx={{ mb: 1.5 }}>Stereochemistry Impact</Typography>

                      {data.stereocenters.length === 0 ? (
                        <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(148,163,184,0.1)' }}>
                          <Typography variant="body2" sx={{ color: 'rgb(148,163,184)' }}>
                            No stereochemistry considerations for achiral molecules.
                          </Typography>
                        </Box>
                      ) : (
                        <Stack spacing={2}>
                          {/* Important Note about 2D properties */}
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(234,179,8,0.1)', border: '1px solid rgba(234,179,8,0.3)' }}>
                            <Typography variant="body2" sx={{ color: '#eab308', fontWeight: 600, mb: 0.5 }}>
                              Why 2D Properties Are Identical
                            </Typography>
                            <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>
                              {data.stereochemistry_context.note_2d_properties}
                            </Typography>
                          </Box>

                          {/* What differs between enantiomers */}
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(168,85,247,0.1)', border: '1px solid rgba(168,85,247,0.3)' }}>
                            <Typography variant="body2" sx={{ color: '#a855f7', fontWeight: 600, mb: 1 }}>
                              What DIFFERS Between Enantiomers (requires experimental data)
                            </Typography>
                            <Stack spacing={0.5}>
                              {data.stereochemistry_context.considerations.map((consideration, i) => (
                                <Typography key={i} variant="caption" sx={{ color: 'rgb(200,200,200)' }}>
                                  {consideration}
                                </Typography>
                              ))}
                            </Stack>
                          </Box>

                          {/* Recommendation */}
                          <Box sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(34,197,94,0.1)', border: '1px solid rgba(34,197,94,0.3)' }}>
                            <Typography variant="body2" sx={{ color: '#22c55e', fontWeight: 600, mb: 0.5 }}>
                              Recommendation
                            </Typography>
                            <Typography variant="caption" sx={{ color: 'rgb(148,163,184)' }}>
                              {data.stereochemistry_context.recommendation}
                            </Typography>
                          </Box>

                          {/* Historical Examples */}
                          {data.stereochemistry_context.historical_examples.length > 0 && (
                            <Box>
                              <Typography variant="caption" sx={{ color: 'rgb(148,163,184)', fontWeight: 600, mb: 1, display: 'block' }}>
                                Historical Examples
                              </Typography>
                              <Stack spacing={0.5}>
                                {data.stereochemistry_context.historical_examples.map((ex, i) => (
                                  <Typography key={i} variant="caption" sx={{ color: 'rgb(120,140,160)' }}>
                                    <strong style={{ color: '#38bdf8' }}>{ex.drug}:</strong> {ex.impact}
                                  </Typography>
                                ))}
                              </Stack>
                            </Box>
                          )}
                        </Stack>
                      )}
                    </Grid>
                  </Grid>
                )}

                {/* Suggestions Tab */}
                {activeTab === 3 && (
                  <Box>
                    {data.optimization_suggestions.length === 0 ? (
                      <Box sx={{ p: 3, borderRadius: 2, backgroundColor: 'rgba(34, 197, 94, 0.1)', border: '1px solid rgba(34, 197, 94, 0.3)', textAlign: 'center' }}>
                        <Typography variant="body1" sx={{ color: '#22c55e' }}>No optimization suggestions - molecule looks good!</Typography>
                      </Box>
                    ) : (
                      <Stack spacing={1.5}>
                        {data.optimization_suggestions.map((s, i) => (
                          <Box key={i} sx={{ p: 2, borderRadius: 2, backgroundColor: 'rgba(56,189,248,0.08)', border: '1px solid rgba(56,189,248,0.2)', display: 'flex', alignItems: 'flex-start', gap: 1.5 }}>
                            <Box sx={{ width: 24, height: 24, borderRadius: '50%', backgroundColor: 'rgba(56,189,248,0.2)', display: 'flex', alignItems: 'center', justifyContent: 'center', flexShrink: 0, color: '#38bdf8', fontSize: '0.75rem', fontWeight: 600 }}>
                              {i + 1}
                            </Box>
                            <Typography variant="body2" sx={{ color: '#e5e7eb' }}>{s.text}</Typography>
                          </Box>
                        ))}
                      </Stack>
                    )}
                  </Box>
                )}
              </Box>
            </Paper>

            {/* Lipinski Details */}
            <Paper
              sx={{
                p: 2.5,
                borderRadius: 3,
                border: '1px solid rgba(148,163,184,0.25)',
                backgroundColor: 'rgba(15,23,42,0.92)',
                boxShadow: '0 8px 32px rgba(0,0,0,0.3)',
              }}
              elevation={0}
            >
              <Typography variant="subtitle1" sx={{ fontWeight: 600, mb: 2 }}>Lipinski Rule-of-Five Details</Typography>
              <Grid container spacing={1.5}>
                {data.lipinski_rules.map((rule) => (
                  <Grid key={rule.name} size={{ xs: 12, sm: 6, md: 3 }}>
                    <Box
                      sx={{
                        p: 1.5,
                        borderRadius: 2,
                        backgroundColor: rule.passed ? 'rgba(34, 197, 94, 0.1)' : 'rgba(239, 68, 68, 0.1)',
                        border: `1px solid ${rule.passed ? 'rgba(34, 197, 94, 0.3)' : 'rgba(239, 68, 68, 0.3)'}`,
                      }}
                    >
                      <Stack direction="row" justifyContent="space-between" alignItems="center">
                        <Typography variant="body2" sx={{ fontWeight: 600 }}>{rule.name}</Typography>
                        <Chip size="small" label={rule.passed ? 'Pass' : 'Fail'} sx={{ backgroundColor: rule.passed ? '#22c55e' : '#ef4444', color: '#000', fontWeight: 600, fontSize: '0.65rem' }} />
                      </Stack>
                      <Typography variant="caption" sx={{ color: 'rgb(148,163,184)', mt: 0.5, display: 'block' }}>{rule.detail}</Typography>
                    </Box>
                  </Grid>
                ))}
              </Grid>
            </Paper>
          </Stack>
        )}

        {/* Empty State */}
        {!data && !loading && (
          <Paper
            sx={{
              p: 6,
              textAlign: 'center',
              borderRadius: 3,
              border: '1px solid rgba(56,189,248,0.2)',
              background: 'linear-gradient(135deg, rgba(15,23,42,0.95), rgba(15,23,42,0.9))',
              boxShadow: '0 8px 32px rgba(56,189,248,0.1)',
            }}
            elevation={0}
          >
            <Box
              sx={{
                width: 80,
                height: 80,
                mx: 'auto',
                mb: 3,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                borderRadius: '50%',
                background: 'radial-gradient(circle, rgba(56,189,248,0.2) 0%, transparent 70%)',
              }}
            >
              <Box sx={{ color: '#38bdf8', transform: 'scale(2.5)' }}>
                <MoleculeIcon />
              </Box>
            </Box>
            <Typography variant="h5" sx={{ color: '#e5e7eb', fontWeight: 600, mb: 2 }}>
              Ready to Analyze Molecules
            </Typography>
            <Typography variant="body1" sx={{ color: 'rgb(148,163,184)', mb: 4, maxWidth: 600, mx: 'auto' }}>
              Enter a SMILES string to get comprehensive drug-likeness predictions including ADMET properties, toxicity alerts, and stereochemistry analysis.
            </Typography>
            <Stack direction="row" spacing={2} justifyContent="center" flexWrap="wrap" useFlexGap>
              <Box sx={{ textAlign: 'left' }}>
                <Typography variant="caption" sx={{ color: 'rgb(100,116,139)', display: 'block', mb: 0.5 }}>
                  What you'll get:
                </Typography>
                <Stack spacing={0.5}>
                  <Typography variant="caption" sx={{ color: '#38bdf8' }}>• QED drug-likeness score</Typography>
                  <Typography variant="caption" sx={{ color: '#a855f7' }}>• ADMET predictions (BBB, CNS, GI)</Typography>
                  <Typography variant="caption" sx={{ color: '#22c55e' }}>• Toxicity & CYP450 alerts</Typography>
                  <Typography variant="caption" sx={{ color: '#eab308' }}>• Stereochemistry context</Typography>
                </Stack>
              </Box>
            </Stack>
          </Paper>
        )}
      </Container>
    </Box>
  )
}

export default App
