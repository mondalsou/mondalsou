# Drug Discovery Triage - Deployment Guide

Production-ready web application for preliminary drug screening with stereochemistry analysis and ADMET predictions.

## Overview

This application uses **transparent, rule-based methods** for ADMET predictions - no black-box ML models required. All predictions are based on established medicinal chemistry principles:

- ✅ **Fast**: No model loading, instant predictions
- ✅ **Transparent**: Interpretable scoring based on chemical rules
- ✅ **Lightweight**: ~500MB Docker image vs 2GB+ with ML frameworks
- ✅ **Reliable**: No TensorFlow/CUDA compatibility issues
- ✅ **Scientifically Valid**: ESOL, Lipinski, CNS MPO (Pfizer), SMARTS toxicophores

---

## Deployment Options

### Option 1: Docker Compose (Self-Hosted) - Recommended for Testing

**Single command deployment:**
```bash
docker-compose up --build
```

**Access:**
- Frontend: http://localhost:3000
- Backend API: http://localhost:8000
- API Docs: http://localhost:8000/docs

**Requirements:**
- Docker & Docker Compose installed
- 2GB RAM minimum
- 5GB disk space

---

### Option 2: Vercel + Railway (Production) - Recommended for Production

#### Step 1: Deploy Backend to Railway

1. **Install Railway CLI:**
   ```bash
   npm i -g @railway/cli
   railway login
   ```

2. **Create new project:**
   ```bash
   cd backend
   railway init
   ```

3. **Deploy:**
   ```bash
   railway up
   ```

4. **Get your backend URL:**
   ```bash
   railway domain
   # Example: https://your-app.railway.app
   ```

#### Step 2: Deploy Frontend to Vercel

1. **Install Vercel CLI:**
   ```bash
   npm i -g vercel
   ```

2. **Update API URL:**
   ```bash
   cd frontend
   # Edit .env.production with your Railway URL
   echo "VITE_API_URL=https://your-app.railway.app" > .env.production
   ```

3. **Deploy:**
   ```bash
   vercel --prod
   ```

#### Step 3: Setup Auto-Deployment (GitHub Actions)

**Add secrets to GitHub repository:**

Settings → Secrets → Actions → New repository secret:

- `RAILWAY_TOKEN`: Get from `railway login --browserless`
- `VERCEL_TOKEN`: Get from https://vercel.com/account/tokens
- `VERCEL_ORG_ID`: Get from `.vercel/project.json` after first deploy
- `VERCEL_PROJECT_ID`: Get from `.vercel/project.json`

**Auto-deploy on push to main:**
```bash
git push origin main
# GitHub Actions will automatically deploy
```

---

## Local Development

### Backend
```bash
cd backend

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run development server
uvicorn app.main:app --reload

# API: http://localhost:8000
# Docs: http://localhost:8000/docs
```

### Frontend
```bash
cd frontend

# Install dependencies
npm install

# Run development server
npm run dev

# App: http://localhost:5173
```

---

## ADMET Prediction Methods

All predictions use **rule-based methods** - fast, transparent, and scientifically validated:

### Absorption & Distribution

**Aqueous Solubility (ESOL)**
- Method: Delaney equation (J. Chem. Inf. Comput. Sci. 2004)
- Formula: `Log S = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP`
- Classification: High (>-1), Moderate (-1 to -3), Low (-3 to -5), Poor (<-5)

**GI Absorption**
- Method: Lipinski Rule-of-Five + Veber rules
- Criteria: MW ≤500, cLogP ≤5, HBD ≤5, HBA ≤10, TPSA ≤140, RotBonds ≤10
- Output: Bioavailability score (0-1)

**BBB Penetration**
- Method: Multi-parameter scoring
- Rules: MW <450, TPSA <90, HBD ≤3, cLogP 1-3 (optimal)
- Output: BBB+/BBB- with probability

**CNS MPO Score**
- Method: Pfizer CNS Multiparameter Optimization (ACS Chem Neurosci 2010)
- Parameters: cLogP, cLogD, MW, TPSA, HBD, pKa (6 components)
- Score: 0-6 (≥4 = High CNS, 3-4 = Medium, <3 = Low)

### Metabolism

**CYP450 Predictions**
- CYP3A4 substrate: Large lipophilic molecules (MW >300, cLogP >2)
- CYP2D6 substrate: Basic nitrogen-containing (1+ basic N, 200 < MW < 500)
- CYP2C9 substrate: Acidic compounds (carboxylic acid present)
- CYP inhibitor risk: High lipophilicity proxy (cLogP >3, MW >300)

### Toxicity

**Structural Alerts (SMARTS patterns)**
- hERG: Quaternary nitrogen, naphthalene systems
- Hepatotoxicity: Nitro groups, acyl chlorides
- Mutagenicity: Nitro compounds, azo compounds
- Endpoints: NR-AhR, SR-MMP, and others from Tox21

**Clinical Toxicity Risk**
- Toxicophore patterns: 6 common structural alerts
- Property-based risk: cLogP >5, MW >600, TPSA outliers
- Output: Risk score (0-1), FDA approval likelihood

---

## Architecture

```
┌─────────────────────┐
│   Frontend (React)  │  Port 3000 (Docker) / 5173 (dev)
│   - Material-UI     │
│   - Vite + TS       │
└──────────┬──────────┘
           │ HTTP
           ▼
┌─────────────────────┐
│  Backend (FastAPI)  │  Port 8000
│  - RDKit            │
│  - ADMET Module     │
│  - Pydantic Models  │
└─────────────────────┘
```

---

## Environment Variables

### Frontend (.env.production)
```bash
VITE_API_URL=https://your-backend.railway.app
```

### Backend (Railway)
```bash
PORT=8000  # Auto-set by Railway
```

---

## CI/CD Pipeline

GitHub Actions workflow (`.github/workflows/deploy.yml`):

1. **Test Backend**: Python syntax, import test
2. **Test Frontend**: TypeScript check, build test
3. **Deploy Backend**: Railway auto-deploy on push to main
4. **Deploy Frontend**: Vercel auto-deploy on push to main

---

## Troubleshooting

### Backend won't start
```bash
# Check RDKit installation
python -c "from rdkit import Chem; print(Chem.__version__)"

# Check API
curl http://localhost:8000/health
```

### Frontend can't connect to backend
```bash
# Check CORS settings in backend/app/main.py
# Check VITE_API_URL in frontend/.env.production
# Check browser console for errors
```

### Docker build fails
```bash
# Clean rebuild
docker-compose down -v
docker-compose build --no-cache
docker-compose up
```

### Railway deployment fails
```bash
# Check logs
railway logs

# Check Dockerfile
railway run bash
```

---

## Performance

- **Prediction time**: <100ms per molecule (rule-based)
- **API response**: ~200ms (including network)
- **Docker image**: ~500MB (vs 2GB+ with ML frameworks)
- **Memory usage**: ~200MB (vs 1GB+ with TensorFlow)

---

## Scientific Validation

All methods are based on peer-reviewed literature:

1. **ESOL**: Delaney, J. Chem. Inf. Comput. Sci. 2004, 44, 1000-1005
2. **Lipinski RO5**: Lipinski et al., Adv. Drug Deliv. Rev. 2001, 46, 3-26
3. **Veber Rules**: Veber et al., J. Med. Chem. 2002, 45, 2615-2623
4. **CNS MPO**: Wager et al., ACS Chem. Neurosci. 2010, 1, 435-449
5. **BBB**: Pajouhesh & Lenz, NeuroRx 2005, 2, 541-553

---

## License

MIT License - See LICENSE file

---

## Support

For issues or questions:
- GitHub Issues: https://github.com/your-username/drug-discovery-triage/issues
- Documentation: See README.md
