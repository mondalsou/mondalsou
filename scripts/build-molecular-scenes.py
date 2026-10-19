"""Fresh illustrative geometries, generated locally; no company assets or results."""
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, Descriptors, QED
out=Path(__file__).resolve().parents[1]/'docs/assets/molecular'
def build(smiles,name,n=1):
 m=Chem.AddHs(Chem.MolFromSmiles(smiles));p=AllChem.ETKDGv3();p.randomSeed=2026;p.pruneRmsThresh=.3
 ids=list(AllChem.EmbedMultipleConfs(m,numConfs=max(n,12 if n>1 else 1),params=p))
 AllChem.MMFFOptimizeMoleculeConfs(m,maxIters=500)
 ring=list(m.GetRingInfo().AtomRings()[0]) if m.GetRingInfo().AtomRings() else list(range(min(3,m.GetNumAtoms())))
 if len(ids)>1:rdMolAlign.AlignMolConformers(m,atomIds=ring)
 heavy=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()>1]
 coords={i:np.array(m.GetConformer(i).GetPositions()) for i in ids};chosen=[ids[0]]
 while len(chosen)<min(n,len(ids)):
  chosen.append(max((i for i in ids if i not in chosen),key=lambda i:min(np.sqrt(np.mean((coords[i][heavy]-coords[j][heavy])**2)) for j in chosen)))
 center=coords[chosen[0]][heavy].mean(axis=0)
 record={'name':name,'smiles':smiles,'atoms':[a.GetSymbol() for a in m.GetAtoms()],'bonds':[[b.GetBeginAtomIdx(),b.GetEndAtomIdx(),b.GetBondTypeAsDouble()] for b in m.GetBonds()],'conformers':[(coords[i]-center).round(4).tolist() for i in chosen], 'descriptors':{'logP':round(Descriptors.MolLogP(m),2),'TPSA':round(Descriptors.TPSA(m),1),'QED':round(QED.qed(m),2)}}
 return m,record
items={}
for key,smi,name,n in [('conformer','CCN(CC)CC(=O)Nc1c(C)cccc1C','Lidocaine',3),('acid','CC(=O)Oc1ccccc1C(=O)O','Aspirin',1),('lead','CC(=O)Nc1ccccc1','Acetanilide',1),('candidate','CC(=O)Nc1ccc(O)cc1','Paracetamol',1),('water','O','Water',1)]:
 m,d=build(smi,name,n);items[key]=d
 if key=='acid':
  for b in m.GetBonds():
   a,z=b.GetBeginAtom(),b.GetEndAtom()
   if {a.GetSymbol(),z.GetSymbol()}=={'O','H'}:d['donor']=[a.GetIdx(),z.GetIdx()] if a.GetSymbol()=='O' else [z.GetIdx(),a.GetIdx()]
# Align the shared scaffold for the illustrative para-H to OH edit.
parent=Chem.MolFromSmiles(items['lead']['smiles'])
child=Chem.MolFromSmiles(items['candidate']['smiles'])
match=child.GetSubstructMatch(parent)
a=np.array(items['lead']['conformers'][0]);b=np.array(items['candidate']['conformers'][0])
x=b[list(match)];y=a[:parent.GetNumAtoms()];xc=x.mean(0);yc=y.mean(0)
u,_,vt=np.linalg.svd((x-xc).T@(y-yc));rotation=u@vt
if np.linalg.det(rotation)<0:u[:,-1]*=-1;rotation=u@vt
items['candidate']['conformers'][0]=((b-xc)@rotation+yc).round(4).tolist()
# Validate physical bonds and consistent topology in every generated conformer.
for key,d in items.items():
 for conf in d['conformers']:
  xyz=np.array(conf)
  for a,b,order in d['bonds']:
   dist=np.linalg.norm(xyz[a]-xyz[b]);assert .85<dist<1.9,(key,a,b,dist)
(out/'scenes.json').write_text(json.dumps(items,separators=(',',':')))
print({k:(len(d['atoms']),len(d['conformers'])) for k,d in items.items()})
