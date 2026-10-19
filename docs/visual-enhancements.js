import * as T from './assets/vendor/three.module.min.js';
// Fresh RDKit geometries. Authored motion illustrates a process, not a simulation.
const data = await fetch(new URL('./assets/molecular/scenes.json', import.meta.url)).then(r=>{if(!r.ok)throw Error('Molecular data unavailable');return r.json();});
const reduced=matchMedia('(prefers-reduced-motion: reduce)');
const V=(x=0,y=0,z=0)=>new T.Vector3(x,y,z);
const colors={C:'#91aac5',H:'#f1f5fa',O:'#f17179',N:'#669ff4'};
const ease=t=>{t=Math.max(0,Math.min(1,t));return t*t*(3-2*t);};
const viewers=[];
function createView(host){
 const kind=host.dataset.scene,canvas=host.querySelector('canvas'),stage=host.querySelector('.scene-stage'),readout=host.querySelector('.scene-readout'),footer=host.querySelector('.scene-footer span');
 const renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:false,preserveDrawingBuffer:true,powerPreference:'low-power'});
 renderer.setPixelRatio(Math.min(devicePixelRatio,1.7));renderer.outputColorSpace=T.SRGBColorSpace;renderer.toneMapping=T.ACESFilmicToneMapping;renderer.toneMappingExposure=1.3;
 const scene=new T.Scene();scene.background=new T.Color('#0a1222');
 scene.add(new T.HemisphereLight('#dcecff','#28314b',2));
 for(const [color,intensity,pos] of [['#fff3e6',3,[3,6,8]],['#7ed9ff',2,[-5,2,-2]],['#b392ff',2,[4,-1,-4]]]){const l=new T.DirectionalLight(color,intensity);l.position.set(...pos);scene.add(l);}
 const camera=new T.PerspectiveCamera(35,1,.1,200);camera.position.set(0,2,21);camera.lookAt(0,0,0);
 const root=new T.Group();scene.add(root);
 const sphere=new T.SphereGeometry(1,24,16),cylinder=new T.CylinderGeometry(1,1,1,12),materials=new Map();
 const mat=(c,opacity=1)=>{const key=c+opacity;if(!materials.has(key))materials.set(key,new T.MeshPhysicalMaterial({color:c,roughness:.25,metalness:.07,clearcoat:.5,transparent:opacity<1,opacity,depthWrite:opacity===1}));return materials.get(key);};
 function ball(parent,pos,r,c,opacity=1){const m=new T.Mesh(sphere,mat(c,opacity));m.position.copy(pos);m.scale.setScalar(r);parent.add(m);return m;}
 function stick(parent,a,b,r,c){const m=new T.Mesh(cylinder,mat(c));parent.add(m);track(m,a,b,r);return m;}
 function track(mesh,a,b,r=.085){const delta=b.clone().sub(a);mesh.position.copy(a).add(b).multiplyScalar(.5);mesh.quaternion.setFromUnitVectors(V(0,1,0),delta.clone().normalize());mesh.scale.set(r,Math.max(.0001,delta.length()),r);}
 function molecular(key,{scale=1,parent=root,conf=0,hide=[],tint=null,hydrogen=true}={}){
  const d=data[key],g=new T.Group();parent.add(g);g.scale.setScalar(scale);
  const points=d.conformers[conf].map(p=>V(...p)),balls=[],bonds=[];
  d.atoms.forEach((e,i)=>{if(hide.includes(i)||(!hydrogen&&e==='H'))return;balls[i]=ball(g,points[i],e==='H'?.19:.34,tint||colors[e]||colors.C);});
  d.bonds.forEach(([a,b,order])=>{
   if(!balls[a]||!balls[b])return;
   const middle=points[a].clone().lerp(points[b],.5);
   bonds.push({a,b,half:0,mesh:stick(g,points[a],middle,.092,tint||colors[d.atoms[a]])});
   bonds.push({a,b,half:1,mesh:stick(g,middle,points[b],.092,tint||colors[d.atoms[b]])});
   if(order===2){const offset=points[b].clone().sub(points[a]).cross(V(0,0,1)).normalize().multiplyScalar(.17);bonds.push({a,b,offset,mesh:stick(g,points[a].clone().add(offset),points[b].clone().add(offset),.045,tint||'#c9d2e0')});}
  });
  function positions(coords,alpha=1){coords.forEach((p,i)=>{if(balls[i])balls[i].position.copy(p);});bonds.forEach(b=>{const a=coords[b.a],z=coords[b.b],mid=a.clone().lerp(z,.5);b.mesh.visible=alpha>.02;if(b.offset)track(b.mesh,a.clone().add(b.offset),z.clone().add(b.offset),.045*alpha);else track(b.mesh,b.half?mid:a,b.half?z:mid,.092*alpha);});}
  return {g,d,points,balls,bonds,positions};
 }
 function waters(parent,count,radius,baseY=0){const list=[];for(let i=0;i<count;i++){const a=i*2.399,rad=radius*(.65+(i%3)*.14);const m=molecular('water',{scale:.34,parent});m.g.position.set(Math.cos(a)*rad,baseY+(i%5-2)*.58,Math.sin(a)*rad);m.g.rotation.set(i*.7,i*.5,0);list.push({g:m.g,p:m.g.position.clone(),i});}return t=>list.forEach(({g,p,i})=>g.position.copy(p).add(V(Math.sin(t*.3+i)*.045,Math.cos(t*.4+i)*.06,0)));}
 const extras=[];
 function glass(){
  const g=new T.CylinderGeometry(3.7,3.5,6.4,64,1,true);extras.push(g);
  const material=new T.MeshPhysicalMaterial({color:'#9fcde8',transparent:true,opacity:.09,roughness:.09,side:T.DoubleSide,depthWrite:false});
  const wall=new T.Mesh(g,material);root.add(wall);
  for(const y of [-3.2,3.2]){const geo=new T.TorusGeometry(y<0?3.5:3.7,.04,8,64);extras.push(geo);const edge=new T.Mesh(geo,mat('#a2d9ee',.65));edge.rotation.x=Math.PI/2;edge.position.y=y;root.add(edge);}
  const fillGeo=new T.CylinderGeometry(3.59,3.47,4.7,64);extras.push(fillGeo);const liquid=new T.Mesh(fillGeo,new T.MeshPhysicalMaterial({color:'#5aa5cf',transparent:true,opacity:.07,depthWrite:false,roughness:.15}));liquid.position.y=-.8;root.add(liquid);
  const surfaceGeo=new T.CircleGeometry(3.59,64);extras.push(surfaceGeo);const surface=new T.Mesh(surfaceGeo,new T.MeshPhysicalMaterial({color:'#81bddf',transparent:true,opacity:.14,side:T.DoubleSide,depthWrite:false}));surface.rotation.x=-Math.PI/2;surface.position.y=1.55;root.add(surface);
  for(let i=0;i<4;i++)stick(root,V(-3.45,-1.8+i*.8,.9),V(-3.1,-1.8+i*.8,.9),.025,'#87bad2');
 }
 let update=()=>{};
 if(kind==='conformer'){
  const molecule=molecular('conformer',{scale:.73});
  const forms=data.conformer.conformers.map(c=>c.map(p=>V(...p)));
  const scattered=molecule.points.map((p,i)=>p.clone().multiplyScalar(1.4).add(V(Math.sin(i*5)*2,Math.cos(i*3)*2,Math.sin(i*2)*2)));
  footer.textContent='Lidocaine · RDKit conformers · illustrative transitions';
  update=t=>{
   const phase=t%24;
   root.rotation.y=.25+Math.sin(t*.18)*.35;root.rotation.x=.15;
   if(phase<4){const k=ease(phase/3.4);molecule.positions(scattered.map((p,i)=>p.clone().lerp(forms[0][i],k)),k);stage.textContent='Assembly';}
   else{const idx=Math.min(2,Math.floor((phase-4)/6.67));molecule.positions(forms[idx]);stage.textContent=`Conformer ${idx+1} / 3`;}
   readout.innerHTML='';
  };
 }else if(kind==='crystal'){
  const assembly=new T.Group();root.add(assembly);const molecules=[];
  for(let i=0;i<8;i++){const m=molecular('acid',{scale:.36,parent:assembly,hydrogen:false});const target=V((i%2-.5)*4.6,(Math.floor(i/2)%2-.5)*3.3,(Math.floor(i/4)-.5)*3.5);m.g.rotation.set(.4,Math.PI*(i%2),.2);molecules.push({m,target});}
  const vertices=[[-2.8,-2.2,-2.3],[2.8,-2.2,-2.3],[2.8,2.2,-2.3],[-2.8,2.2,-2.3],[-2.8,-2.2,2.3],[2.8,-2.2,2.3],[2.8,2.2,2.3],[-2.8,2.2,2.3]];
  for(const [a,b] of [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])stick(assembly,V(...vertices[a]),V(...vertices[b]),.017,'#537ba4');
  footer.textContent='Aspirin · illustrative packing, not a predicted polymorph';
  update=t=>{const p=t%16,k=ease(p/5);assembly.rotation.set(.25,t*.09,.08);molecules.forEach(({m,target},i)=>m.g.position.copy(target).multiplyScalar(1+(1-k)*.45));stage.textContent=p<5?'Assembling':'Crystal packing';};
 }else if(kind==='pka'){
  glass();const d=data.acid,[oi,hi]=d.donor;
  const acid=molecular('acid',{scale:.42,hide:[hi],hydrogen:false});
  const direction=V(...d.conformers[0][hi]).sub(V(...d.conformers[0][oi])).normalize();
  acid.g.quaternion.setFromUnitVectors(direction,V(1,0,0));
  const localO=V(...d.conformers[0][oi]).multiplyScalar(.42).applyQuaternion(acid.g.quaternion);
  acid.g.position.copy(V(-.4,-.1,.8).sub(localO));acid.g.updateMatrixWorld(true);
  const o=V(...d.conformers[0][oi]).applyMatrix4(acid.g.matrixWorld),start=V(...d.conformers[0][hi]).applyMatrix4(acid.g.matrixWorld);
  const acceptor=V(1.8,-.1,.8),end=acceptor.clone().add(V(-.52,.17,0));
  ball(root,acceptor,.3,colors.O);for(const delta of [[.38,.4,.05],[.3,-.44,-.12]]){const h=acceptor.clone().add(V(...delta));ball(root,h,.19,colors.H);stick(root,acceptor,h,.07,'#b5c7d9');}
  const proton=ball(root,start,.22,'#ffcf69'),halo=ball(root,start,.34,'#ffce6b',.14);
  const oldBond=stick(root,o,start,.07,'#f6c76d'),newBond=stick(root,end,acceptor,.07,'#f6c76d');const solvent=new T.Group();root.add(solvent);for(const x of [-2.5,2.5]){const water=molecular('water',{scale:.20,parent:solvent});water.g.position.set(x,-2,-1.7);water.g.rotation.z=x*.2;}const drift=()=>{};
  footer.textContent='Aspirin in water · illustrative proton transfer';
  update=t=>{const p=t%18,k=ease((p-3)/4)*(1-ease((p-13)/4));const h=start.clone().lerp(end,k);h.y+=Math.sin(k*Math.PI)*.22;proton.position.copy(h);halo.position.copy(h);oldBond.visible=k<.5;newBond.visible=k>=.5;track(oldBond,o,h,.07);track(newBond,h,acceptor,.07);stage.textContent=k<.1?'Before transfer':k>.9?'After transfer':'Proton transfer';readout.innerHTML='';drift(t);root.rotation.y=0;};
 }else if(kind==='property'){
  // A qualitative shake-flask view: two immiscible phases and a neutral solute.
  const vessel=new T.Group();root.add(vessel);vessel.position.y=.35;
  const wallGeo=new T.CylinderGeometry(3.15,3.15,6.2,64,1,true);extras.push(wallGeo);
  vessel.add(new T.Mesh(wallGeo,new T.MeshPhysicalMaterial({color:'#c5d9ef',transparent:true,opacity:.075,roughness:.12,side:T.DoubleSide,depthWrite:false})));
  for(const y of [-3.1,3.1]){const g=new T.TorusGeometry(3.15,.035,8,64);extras.push(g);const rim=new T.Mesh(g,mat('#bacce0',.65));rim.position.y=y;rim.rotation.x=Math.PI/2;vessel.add(rim);}
  const phaseGeo=new T.CylinderGeometry(3.07,3.07,2.65,64);extras.push(phaseGeo);
  for(const [y,color] of [[1.35,'#e2ae59'],[-1.35,'#438fce']]){const liquid=new T.Mesh(phaseGeo,new T.MeshPhysicalMaterial({color,transparent:true,opacity:.23,roughness:.22,depthWrite:false}));liquid.position.y=y;vessel.add(liquid);}
  const discGeo=new T.CircleGeometry(3.07,64);extras.push(discGeo);const interfacePlane=new T.Mesh(discGeo,new T.MeshBasicMaterial({color:'#d8dccf',transparent:true,opacity:.2,side:T.DoubleSide,depthWrite:false}));interfacePlane.rotation.x=-Math.PI/2;vessel.add(interfacePlane);
  const residents=[];
  for(const [i,pos] of [[0,[-1.35,1.7,-.4]],[1,[1.2,1.35,-.7]],[2,[-.4,2.1,1]],[3,[-1.25,-1.55,-.2]]]){const m=molecular('lead',{scale:.24,parent:vessel,hydrogen:false});m.g.position.set(...pos);m.g.rotation.set(.2,i*.9,.3);residents.push({g:m.g,p:m.g.position.clone(),i});}
  const exchange=molecular('lead',{scale:.35,parent:vessel,hydrogen:false});
  const labels=document.createElement('div');labels.className='logp-phase-labels';labels.innerHTML='<span class="octanol-label">n-Octanol</span><span class="water-label">Water</span>';host.appendChild(labels);
  const membrane=new T.Group();root.add(membrane);membrane.visible=false;
  // Two leaflets: polar heads face water, hydrocarbon tails face inward.
  // Remove only foreground lipids in the centre as a viewing cutaway, not a pore.
  const lipids=[];
  for(let i=0;i<11;i++)for(let j=0;j<4;j++){
   const x=-3.5+i*.7,z=-1.2+j*.8;
   if(Math.abs(x)<1.1&&z>.2)continue;
   for(const side of [-1,1]){
    const lipid=new T.Group();membrane.add(lipid);lipid.position.set(x,side*.95,z);
    ball(lipid,V(),.19,'#73c8cc');
    for(const offset of [-.11,.11]){
     const a=V(offset,-side*.15,0),mid=V(offset+.05,-side*.48,.06),end=V(offset-.03,-side*.79,.02);
     stick(lipid,a,mid,.055,'#c9af7e');stick(lipid,mid,end,.05,'#c9af7e');
    }
    lipids.push({g:lipid,x,z,side});
   }
  }
  const permeant=molecular('lead',{scale:.44,parent:membrane});
  for(let i=0;i<6;i++){
   const water=molecular('water',{scale:.25,parent:membrane});water.g.position.set(-2.8+(i%3)*2.8,i<3?2.7:-2.7,-.65);water.g.rotation.set(i*.3,i,0);
  }
  const tabs=document.createElement('div');tabs.className='scene-tabs';tabs.setAttribute('role','group');tabs.setAttribute('aria-label','ADMET illustration');
  let mode='logp';
  for(const [value,title] of [['logp','logP'],['permeability','Permeability']]){
   const tab=document.createElement('button');tab.type='button';tab.textContent=title;tab.setAttribute('aria-label',`Show ${title} illustration`);tab.setAttribute('aria-pressed',String(mode===value));
   tab.addEventListener('click',()=>{mode=value;time=0;dirty=true;[...tabs.children].forEach(button=>button.setAttribute('aria-pressed',String(button===tab)));update(time);});tabs.appendChild(tab);
  }
  host.appendChild(tabs);
  update=t=>{
   const isLogP=mode==='logp';vessel.visible=isLogP;membrane.visible=!isLogP;labels.hidden=!isLogP;
   if(isLogP){
    stage.textContent='Partitioning';footer.textContent='Neutral solute · qualitative illustration';
    readout.innerHTML='<div class="logp-equation">logP = log₁₀(c<sub>octanol</sub> / c<sub>water</sub>)</div>';
    const p=t%18,k=ease((p-2)/5)*(1-ease((p-11)/5));
    exchange.g.position.set(1.1,-1.6+3.2*k,1.5);exchange.g.rotation.set(.15,.35+Math.sin(t*.35)*.25,.2);
    residents.forEach(({g,p,i})=>g.position.copy(p).add(V(Math.sin(t*.4+i)*.07,Math.cos(t*.3+i)*.06,0)));
   }else{
    const p=t%22,k=ease((p-2)/15);
    permeant.g.position.set(.12,2.6-5.2*k,.5);permeant.g.rotation.set(.12,Math.sin(t*.2)*.25,.9);
    // The return is a gentle rewind of the illustration, not a second mechanism.
    if(p>19)permeant.g.position.y=-2.6+5.2*ease((p-19)/3);
    lipids.forEach(({g,x,z,side},i)=>{g.position.set(x+Math.sin(t*.5+i)*.025,side*.95,z);});
    stage.textContent=k<.3?'Approach':k<.7?'Crossing':'Release';
    footer.textContent='Passive diffusion · schematic, not a predicted rate';
    readout.innerHTML='';
   }
   root.rotation.y=0;
  };
 }else if(kind==='agent'){
  const bot=new T.Group();root.add(bot);bot.position.set(-3.5,2.15,0);
  const headGeo=new T.BoxGeometry(1.65,1.25,.9);extras.push(headGeo);const head=new T.Mesh(headGeo,mat('#7189aa'));bot.add(head);
  const faceGeo=new T.BoxGeometry(1.39,.8,.06);extras.push(faceGeo);const face=new T.Mesh(faceGeo,mat('#0b2036'));face.position.z=.49;bot.add(face);
  for(const x of [-.4,.4])ball(bot,V(x,.08,.55),.14,'#65edee');
  stick(bot,V(0,.65,0),V(0,1.02,0),.045,'#b6c7dc');const antenna=ball(bot,V(0,1.12,0),.15,'#b697ff');
  const thoughts=[];for(let i=0;i<3;i++)thoughts.push(ball(root,V(-1.8+i*.55,3.05,0),.1,'#bca0ef'));
  const a=molecular('lead',{scale:.84}),b=molecular('candidate',{scale:.84});a.g.position.set(1,.2,0);b.g.position.copy(a.g.position);
  const ox=data.candidate.atoms.map((e,i)=>e==='O'?i:-1).filter(i=>i>=0).pop();
  const attachedH=data.candidate.bonds.find(([x,y])=>(x===ox&&data.candidate.atoms[y]==='H')||(y===ox&&data.candidate.atoms[x]==='H'));
  const hi=attachedH[0]===ox?attachedH[1]:attachedH[0];
  const joint=V(-2,1,1),armA=stick(root,V(-3.5,1.55,0),joint,.1,'#8297b6'),armB=stick(root,joint,V(0,0,0),.075,'#a5bbd4');
  ball(root,joint,.18,'#8badcd');const gripper=ball(root,V(),.18,'#78e3e7');
  const changeHalo=ball(b.g,b.points[ox],.5,'#f09dcf',.12);
  const ribbon=document.createElement('div');ribbon.className='agent-steps';ribbon.innerHTML='<span>Plan</span><span>Edit</span><span>Check</span><span>Refine</span>';host.appendChild(ribbon);
  footer.textContent='Illustrative design agent · not a live AI run';
  const messages=['Keep the scaffold. Explore a polar substituent.','Replace a para hydrogen with an OH group.','Check the structure and compare design constraints.','Use tool feedback to choose the next modification.'];
  update=t=>{const p=t%32,step=Math.floor(p/8),within=p%8;const rotation=.15; a.g.rotation.set(.12,rotation,0);b.g.rotation.copy(a.g.rotation);
   a.g.visible=step===0;b.g.visible=step!==0;
   const k=step===1?ease(within/4):1;
   const coords=b.points.map((v,i)=>v.clone().add((i===ox||i===hi)?V(1.8*(1-k),2.4*(1-k),0):V()));b.positions(coords);
   b.g.updateMatrixWorld(true);const target=coords[ox].clone().applyMatrix4(b.g.matrixWorld);
   const idle=V(-1.7,1.9,1);const tip=step===1?target:idle;track(armB,joint,tip,.075);gripper.position.copy(tip);
   changeHalo.position.copy(coords[ox]);changeHalo.visible=step===1||step===2;
   thoughts.forEach((dot,i)=>{dot.visible=step===0||step===3;dot.scale.setScalar(.08+.05*(.5+.5*Math.sin(t*3-i)));});
   antenna.scale.setScalar(.13+.025*Math.sin(t*2));
   stage.textContent=['01 · Agent plans an edit','02 · Agent modifies the molecule','03 · Agent calls chemistry tools','04 · Feedback guides the next design'][step];
   readout.innerHTML=`<strong>${step===0?'Acetanilide → candidate design':step===1?'Proposed OH substitution':step===2?'Validate → evaluate → compare':'Next design cycle'}</strong><br>${messages[step]}`;
   [...ribbon.children].forEach((e,i)=>e.classList.toggle('active',i===step));
  };
 }

 let time=reduced.matches?7:0,paused=reduced.matches,visible=false,dirty=true,dragging=false,lastX=0,turn=0;
 const button=host.querySelector('.visual-pause');
 function syncButton(){button.setAttribute('aria-pressed',String(paused));button.textContent=paused?'▷':'Ⅱ';button.setAttribute('aria-label',`${paused?'Play':'Pause'} ${kind} animation`);}
 button.addEventListener('click',()=>{paused=!paused;syncButton();dirty=true;});syncButton();
 host.querySelector('.scene-replay').addEventListener('click',()=>{time=0;dirty=true;});
 new IntersectionObserver(entries=>{visible=entries[0].isIntersecting;dirty=true;},{rootMargin:'60px'}).observe(host);
 new ResizeObserver(()=>{const w=host.clientWidth,h=host.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.position.z=kind==='agent'?23:Math.max(16,14/camera.aspect);camera.updateProjectionMatrix();dirty=true;}).observe(host);
 canvas.addEventListener('pointerdown',e=>{if(e.pointerType==='touch')return;dragging=true;lastX=e.clientX;canvas.setPointerCapture(e.pointerId);});
 canvas.addEventListener('pointermove',e=>{if(dragging){turn+=(e.clientX-lastX)*.009;lastX=e.clientX;dirty=true;}});
 canvas.addEventListener('pointerup',()=>dragging=false);canvas.addEventListener('pointercancel',()=>dragging=false);
 reduced.addEventListener('change',()=>{paused=reduced.matches;syncButton();dirty=true;});
 canvas.addEventListener('webglcontextlost',event=>{event.preventDefault();paused=true;stage.textContent='3D view paused — reload to restore';syncButton();});
 update(time);host.title=footer.textContent;renderer.render(scene,camera);host.dataset.ready='true';
 return dt=>{if(!visible||document.hidden)return;if(!paused&&!reduced.matches){time+=dt;dirty=true;}if(dirty){update(time);root.rotation.y+=turn;renderer.render(scene,camera);root.rotation.y-=turn;dirty=false;}};
}
for(const host of document.querySelectorAll('.science-visual')){
 try{viewers.push(createView(host));}catch(error){console.error(error);host.querySelector('.scene-stage').textContent='3D rendering unavailable in this browser';}
}
let previous=0;
function frame(now){requestAnimationFrame(frame);if(now-previous<40)return;const dt=previous?Math.min((now-previous)/1000,.08):0;previous=now;viewers.forEach(draw=>draw(dt));}
requestAnimationFrame(frame);
const menuToggle=document.getElementById('menu-toggle'),menu=document.getElementById('mobile-menu');
menuToggle.setAttribute('aria-controls','mobile-menu');menuToggle.setAttribute('aria-expanded','false');
new MutationObserver(()=>menuToggle.setAttribute('aria-expanded',String(menu.classList.contains('open')))).observe(menu,{attributes:true,attributeFilter:['class']});
document.addEventListener('keydown',e=>{if(e.key==='Escape'&&menu.classList.contains('open')){menu.classList.remove('open');menuToggle.focus();}});
