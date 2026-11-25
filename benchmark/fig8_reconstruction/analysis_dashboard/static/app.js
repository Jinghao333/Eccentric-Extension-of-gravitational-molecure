'use strict';
/* Client is exclusively presentation: scientific formulae live in Python core.py.
   Plotly is served locally from the installed Python package (works offline). */
const $ = id => document.getElementById(id);
const state={run:null,overview:null,atlas:null,geometry:null,channel:null,
             angular:null, selected:null, page:'overview'};
const color=['#7bd1dd','#edac67','#bc9ff0','#92d3b1','#ef8494','#8baeff','#d8c086','#b3cbd6'];
const C={paper:'#141f2e',grid:'#2c3e50',ink:'#e9f3fd',muted:'#9db3c8'};
function fmt(x,n=4){if(x===null||x===undefined||!Number.isFinite(+x))return '—';let v=Number(x);return v!==0&&(Math.abs(v)<.001||Math.abs(v)>=1e5)?v.toExponential(3):v.toLocaleString(undefined,{maximumSignificantDigits:n});}
function percent(x){return x==null?'—':(x*100).toFixed(x>=.01?2:4)+'%';}
function escapeHtml(x){return String(x).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;').replaceAll("'",'&#039;');}
function showError(s){$('error').hidden=!s;$('error').textContent=s||'';}
function layout(extra={}){return Object.assign({paper_bgcolor:C.paper,plot_bgcolor:C.paper,font:{family:'Inter, system-ui, sans-serif',color:C.ink,size:11},margin:{l:58,r:22,t:24,b:52},hovermode:'closest',legend:{orientation:'h',y:-.25,x:0,font:{size:10}},xaxis:{gridcolor:C.grid,zerolinecolor:C.grid,automargin:true},yaxis:{gridcolor:C.grid,zerolinecolor:C.grid,automargin:true}},extra);}
const config={responsive:true,displaylogo:false,modeBarButtonsToRemove:['lasso2d','select2d','toggleSpikelines'],toImageButtonOptions:{format:'svg',scale:1}};
function draw(id,data,lay,onClick){let el=$(id);Plotly.react(el,data,layout(lay),config).then(()=>{if(el.removeAllListeners)el.removeAllListeners('plotly_click');if(onClick)el.on('plotly_click',onClick);});}
async function getJSON(path){let res=await fetch(path);let o=await res.json();if(!res.ok)throw new Error(o.error||('HTTP '+res.status));return o;}
function query(url,params={}){return url+'?'+new URLSearchParams(Object.assign({run:state.run},params));}
function card(label,v,caption){return `<div class="kpi"><label>${escapeHtml(label)}</label><div class="v">${escapeHtml(v)}</div>${caption?`<small>${escapeHtml(caption)}</small>`:''}</div>`;}
function kv(label,value){return `<div class="kv"><span>${escapeHtml(label)}</span><b>${escapeHtml(value)}</b></div>`;}
function setPage(page){state.page=page;document.querySelectorAll('.page').forEach(p=>p.hidden=p.id!=='page-'+page);document.querySelectorAll('.sidebar button[data-page]').forEach(b=>b.classList.toggle('active',b.dataset.page===page));if(page==='deep'&&state.selected&&!state.channel)loadChannel(...state.selected);setTimeout(()=>{document.querySelectorAll('.page:not([hidden]) .js-plotly-plot').forEach(p=>Plotly.Plots.resize(p));},100);}
async function loadRuns(){try{let r=await getJSON('/api/runs');let sel=$('runSelector');let before=state.run;sel.innerHTML=r.runs.length?r.runs.map(s=>`<option value="${escapeHtml(s)}">${escapeHtml(s)}</option>`).join(''):'<option value="">No completed runs</option>';if(r.runs.length){sel.value=r.runs.includes(before)?before:r.runs[0];await loadRun(sel.value);}else showError('No completed v4 run in this directory. Supply a runs directory with config.json + physics.h5 for each simulation.');}catch(err){showError('Cannot list runs: '+err.message);}}
async function loadRun(name){state.run=name;state.overview=null;state.atlas=null;state.geometry=null;state.channel=null;state.angular=null;state.selected=null;showError('');$('pointLabel').textContent='Reading '+name;
  try{let [o,a,g]=await Promise.all([getJSON(query('/api/overview')),getJSON(query('/api/atlas')),getJSON(query('/api/geometry',{n:180}))]);state.overview=o;state.atlas=a;state.geometry=g;
    $('pointLabel').textContent=`α=${o.parameters.alpha} · e=${o.parameters.e} · p=${o.parameters.p}`;
    renderOverview();renderAtlas();fillSelector();if(a.top.length)selectChannel(a.top[0].ell,a.top[0].m,a.top[0].g,false);
  }catch(err){showError('Run cannot be interpreted: '+err.message);$('pointLabel').textContent='Incomplete or invalid run';}}
function renderOverview(){let {parameters:p,totals:t,scales:s,warnings,radial}=state.overview;
  $('kpis').innerHTML=card('Orbital period T',fmt(p.period)+' M','Prescribed Newtonian orbit')+
    card('Total orbital ionization power',fmt(t.orbital_power),'Raw signed orbital energy loss')+
    card('Fig. 8 normalized power',fmt(t.normalized_power),'q⁻² (M/Mb) P')+
    card('Fractional depletion rate Γ',fmt(t.fractional_cloud_depletion_rate)+' /M','Leading-order charge convention');
  const pData=[['Primary mass M',fmt(p.M)],['Companion mass mp=qM',fmt(p.m_p)],['Cloud mass Mb',fmt(p.M_b)],['Scalar mass μ',fmt(p.mu)],['α=Mμ',fmt(p.alpha)],['q=mp/M',fmt(p.q)],['Initial state',`|${p.N_i},${p.ell_i},${p.m_i}⟩`],['Binding energy Ei',fmt(p.E_i)],['Semilatus rectum p',fmt(p.p)+' M'],['Eccentricity e',fmt(p.e)],['Semi-major axis a',fmt(p.a)+' M'],['Periapsis rp',fmt(p.periapsis)+' M'],['Apoapsis ra',fmt(p.apoapsis)+' M'],['Ωorb',fmt(p.omega)+' /M'],['Orbital period T',fmt(p.period)+' M'],['Time of periapsis tp',fmt(p.t_p)+' M'],['Periapsis angle ϖ',fmt(p.varpi)+' rad'],['Open computed channels',fmt(t.n_open)],['Sampled-zero channels',fmt(t.n_sampled_zero)],['Charge current ∑ J',fmt(t.outgoing_number_current)],['1/Γ (instantaneous)',fmt(t.approximate_depletion_timescale)+' M']];
  $('paramsTable').innerHTML=pData.map(x=>kv(...x)).join('');
  $('scalesTable').innerHTML=[['Bohr radius aB',fmt(s.bohr_radius)+' M'],['Equatorial density peak',fmt(s.density_peak)+' M'],['Shell-mass peak',fmt(s.shell_peak)+' M'],['Mean cloud radius',fmt(s.mean_radius)+' M'],['50% mass radius',fmt(s.r50)+' M'],['90% mass radius',fmt(s.r90)+' M']].map(x=>kv(...x)).join('');
  $('warnings').innerHTML=warnings.map(x=>`<p>${escapeHtml(x)}</p>`).join('')+`<p>Mass profile analytic integral / configured Mb = ${fmt(s.mass_integral_from_analytic_profile/p.M_b,7)}. A mismatch indicates a normalization or sampling issue.</p>`;
  renderDensity();renderRadial();
}
function circleTrace(r,name,clr,dash){let phi=Array.from({length:160},(_,j)=>j/159*Math.PI*2);return {type:'scatter',mode:'lines',name,x:phi.map(a=>r*Math.cos(a)),y:phi.map(a=>r*Math.sin(a)),line:{color:clr,width:1,dash:dash||'dot'},hoverinfo:'skip'};}
function markerAt(t){let path=state.overview.orbit;let m=2*Math.PI*t;let i=path.M_a.findIndex(v=>v>=m);i=i<1?1:i;let before=path.M_a[i-1],after=path.M_a[i];let frac=(m-before)/Math.max(after-before,1e-16);return {x:path.x[i-1]+frac*(path.x[i]-path.x[i-1]),y:path.y[i-1]+frac*(path.y[i]-path.y[i-1]),r:path.r[i-1]+frac*(path.r[i]-path.r[i-1]),E:path.E[i-1]+frac*(path.E[i]-path.E[i-1])};}
function renderDensity(){if(!state.geometry||!state.overview)return;let s=state.overview.scales,p=state.overview.parameters,g=state.geometry;let linear=$('densityScale').value==='linear';
 let z=g.rho.map(row=>row.map(x=>linear?x:(x>0?Math.log10(x):null)));let scale=linear?'Mass density':'log₁₀ mass density';let phase=markerAt(+$('orbitalPhase').value);
 let traces=[{type:'heatmap',x:g.coords,y:g.coords,z,colorscale:'Viridis',colorbar:{title:{text:scale,side:'right'},thickness:11},hovertemplate:'x=%{x:.2f} M<br>y=%{y:.2f} M<br>density scale=%{z:.4g}<extra></extra>',name:scale},
 {type:'scatter',mode:'lines',name:'Kepler orbit',x:state.overview.orbit.x,y:state.overview.orbit.y,line:{color:'#fae4a0',width:2}},
 {type:'scatter',mode:'markers',name:'Primary',x:[0],y:[0],marker:{color:'#0c1018',size:12,line:{color:'#dce5f4',width:1.3}}},
 {type:'scatter',mode:'markers',name:'Companion',x:[phase.x],y:[phase.y],marker:{color:'#f0a575',size:12,line:{color:'#fdf3df',width:1.5}}},
 circleTrace(s.density_peak,'Density peak',color[1]),circleTrace(s.shell_peak,'Shell peak',color[2]),circleTrace(s.r50,'50% mass radius','#8bcfa8')];
 let extent=g.extent;
 draw('geometryPlot',traces,{margin:{l:44,r:18,t:12,b:40},xaxis:{title:'x / M',range:[-extent,extent],scaleanchor:'y',scaleratio:1,gridcolor:C.grid},yaxis:{title:'y / M',range:[-extent,extent],gridcolor:C.grid},legend:{orientation:'h',x:0,y:-.15,font:{size:9}},uirevision:state.run+'-geometry'});
 updatePhase();}
function updatePhase(){if(!state.overview)return;let t=+$('orbitalPhase').value;let m=markerAt(t);$('orbitalReadout').textContent=`t/T=${t.toFixed(3)} · r=${fmt(m.r)} M · E=${fmt(m.E)} rad`;let el=$('geometryPlot');if(el.data&&el.data.length>3)Plotly.restyle(el,{x:[[m.x]],y:[[m.y]]},[3]);}
function renderRadial(){if(!state.overview)return;let d=state.overview.radial,s=state.overview.scales,p=state.overview.parameters;
 let key=$('radialChoice').value;let labels={rho_equatorial:'Equatorial mass density',dM_dr:'Mass per radial shell',enclosed_fraction:'Enclosed mass fraction'};
 let vertical=[['rp',p.periapsis,'#e2ac70'],['ra',p.apoapsis,'#e2ac70'],['ρ peak',s.density_peak,'#80c6ce'],['shell peak',s.shell_peak,'#bc9fee'],['r50',s.r50,'#8cd3b5']];
 let t=[{type:'scatter',mode:'lines',x:d.radius,y:d[key],name:labels[key],line:{color:color[0],width:2.4}}];
 for(let [name,r,c] of vertical){t.push({type:'scatter',mode:'lines',x:[r,r],y:[Math.min(...d[key]),Math.max(...d[key])],name,visible:'legendonly',line:{color:c,width:1.5,dash:'dash'},hoverinfo:'name'});}
 let isLog=key!=='enclosed_fraction';draw('radialPlot',t,{xaxis:{title:'r / M',gridcolor:C.grid},yaxis:{title:labels[key],type:isLog?'log':'linear',dtick:isLog?1:undefined,gridcolor:C.grid,autorange:true},legend:{orientation:'h',y:-.22,x:0,font:{size:9}},uirevision:state.run+'-radial-'+key});
}
function renderAtlas(){if(!state.atlas)return;renderFig11();renderTop();renderAngular();renderAngularModes();renderHarmonics();renderCumulative();}
function selectAngular(ell,m){state.angular={ell,m};renderHarmonics();$('selectedChannel').textContent=`Angular mode (ℓf,mf)=(${ell},${m}). Select a harmonic to open its response.`;}
function selectChannel(ell,m,g,jump=false){state.selected=[+ell,+m,+g];state.channel=null;state.angular={ell:+ell,m:+m};$('selectedChannel').textContent=`Selected open channel (${ell},${m},${g}). Click Channel physics to explore.`;if(state.atlas)renderHarmonics();let sel=$('channelSelector');sel.value=state.selected.join(',');if(jump){setPage('deep');loadChannel(...state.selected);}}
function renderFig11(){let d=state.atlas;let tr=d.fig11.map((s,i)=>({type:'scatter',mode:'markers',name:s.name,x:s.points.filter(p=>p.power>0).map(p=>p.ell),y:s.points.filter(p=>p.power>0).map(p=>p.power),customdata:s.points.filter(p=>p.power>0).map(p=>[p.ell,p.m]),marker:{symbol:['diamond','square','triangle-up','circle'][i],size:9,color:color[i],line:{width:1,color:'#eff4ff'}}}));
 draw('fig11Plot',tr,{xaxis:{title:'Outgoing degree ℓf',dtick:1,gridcolor:C.grid},yaxis:{title:'Normalized orbital power',type:'log',dtick:1,autorange:true,gridcolor:C.grid},legend:{orientation:'h',x:0,y:-.2,font:{size:9}},uirevision:state.run+'-fig11'},ev=>{let p=ev.points[0];if(p.customdata)selectAngular(...p.customdata);});}
function renderTop(){let d=state.atlas;let top=d.top.slice(0,10).reverse();$('topSummary').textContent=`Top 10 individual channels account for ${percent(d.top_10_fraction)} of full-run normalized power.`;
 draw('topPlot',[{type:'bar',orientation:'h',y:top.map(p=>`(${p.ell},${p.m},${p.g})`),x:top.map(p=>p.power),customdata:top.map(p=>[p.ell,p.m,p.g,p.fraction,p.delta_g]),marker:{color:color[0]},hovertemplate:'(ℓf,mf,g)=(%{customdata[0]},%{customdata[1]},%{customdata[2]})<br>Power=%{x:.5g}<br>Share=%{customdata[3]:.3%}<br>Δg=%{customdata[4]:.5g}<extra></extra>'}],
 {xaxis:{type:'log',title:'Normalized ionization power',gridcolor:C.grid,autorange:true},yaxis:{automargin:true,gridcolor:C.grid},margin:{l:83,r:12,t:20,b:60},showlegend:false,uirevision:state.run+'-top'},ev=>{let d=ev.points[0].customdata;selectChannel(d[0],d[1],d[2],true);});}
function renderAngular(){let d=state.atlas.angular_heatmap;
 let log=d.power.map(row=>row.map(p=>p>0?Math.log10(p):null));let zs=[];for(let j=0;j<d.m.length;j++)for(let i=0;i<d.ell.length;i++)if(d.status[j][i]==='zero_or_closed')zs.push([d.ell[i],d.m[j]]);
 let tr=[{type:'heatmap',x:d.ell,y:d.m,z:log,colorscale:'Viridis',colorbar:{title:{text:'log₁₀ power',side:'right'},thickness:11},hovertemplate:'ℓf=%{x}, mf=%{y}<br>log₁₀ power=%{z:.3f}<extra></extra>'}];
 if(zs.length)tr.push({type:'scatter',mode:'markers',name:'No positive flux',x:zs.map(x=>x[0]),y:zs.map(x=>x[1]),marker:{symbol:'circle-open',color:'#849db2',size:7},hovertemplate:'ℓf=%{x}, mf=%{y}<br>Zero/closed<extra></extra>'});
 draw('angularPlot',tr,{xaxis:{title:'Outgoing degree ℓf',dtick:1,gridcolor:C.grid},yaxis:{title:'Azimuthal mf',dtick:2,gridcolor:C.grid},margin:{l:53,r:22,t:20,b:54},legend:{orientation:'h',x:0,y:-.22},uirevision:state.run+'-angular'},ev=>{let p=ev.points[0];if(Number.isFinite(p.x)&&Number.isFinite(p.y))selectAngular(Math.round(p.x),Math.round(p.y));});}
function renderAngularModes(){let modes=state.atlas.angular_modes||[];
 let labels=modes.map(p=>`(${p.ell},${p.m})`);
 let custom=modes.map(p=>[p.ell,p.m,p.raw_power,p.share,p.harmonic_count]);
 $('angularModeSummary').textContent=`${modes.length} positive angular pairs. Power is summed over all stored open g and shares are relative to the full-run ionization power.`;
 let traces=[{type:'bar',name:'Normalized power',x:labels,y:modes.map(p=>p.power),customdata:custom,marker:{color:modes.map(p=>p.share),colorscale:'Viridis',colorbar:{title:{text:'Share',side:'right'},tickformat:'.1%',thickness:11}},hovertemplate:'(\u2113f,mf)=(%{customdata[0]},%{customdata[1]})<br>\u03a3g normalized power=%{y:.6g}<br>\u03a3g raw power=%{customdata[2]:.6g}<br>Share=%{customdata[3]:.3%}<br>Open g count=%{customdata[4]}<extra></extra>'},
 {type:'scatter',mode:'lines+markers',name:'Share of total',x:labels,y:modes.map(p=>p.share),customdata:custom,yaxis:'y2',line:{color:color[3],width:1.5},marker:{size:5},hovertemplate:'(\u2113f,mf)=(%{customdata[0]},%{customdata[1]})<br>Share=%{y:.3%}<extra></extra>'}];
 draw('angularModePlot',traces,{xaxis:{title:'Angular pair (\u2113f, mf), descending contribution',tickangle:-55,gridcolor:C.grid,automargin:true},yaxis:{title:'Summed normalized ionization power',type:'log',dtick:1,autorange:true,gridcolor:C.grid},yaxis2:{title:'Share of total',overlaying:'y',side:'right',tickformat:'.1%',rangemode:'tozero',showgrid:false},legend:{orientation:'h',x:0,y:-.28,font:{size:9}},margin:{l:70,r:72,t:20,b:100},uirevision:state.run+'-angular-modes'},ev=>{let p=ev.points[0];if(p.customdata)selectAngular(p.customdata[0],p.customdata[1]);});}
function renderHarmonics(){if(!state.atlas)return;let d=state.atlas;let selected=state.angular;
 let filtered=$('harmonicFilter').checked&&selected;
 let rows=filtered?d.harmonics.map(h=>({...h,power:d.channels.filter(c=>c.ell===selected.ell&&c.m===selected.m&&c.g===h.g).reduce((s,c)=>s+c.power,0)})):d.harmonics;
 $('harmonicSub').textContent=filtered?`Only angular mode (ℓf,mf)=(${selected.ell},${selected.m}); click a nonzero g to inspect it.`:'Summed over all outgoing angular channels';
 let active=rows.filter(x=>x.power>0);
 draw('harmonicPlot',[{type:'bar',x:active.map(x=>x.g),y:active.map(x=>x.power),marker:{color:color[2]},customdata:active.map(x=>x.delta_g),hovertemplate:'g=%{x}<br>Power=%{y:.5g}<br>Δg=%{customdata:.5g}<extra></extra>'}],
 {xaxis:{title:'Temporal harmonic g',dtick:1,gridcolor:C.grid},yaxis:{title:'Normalized power',type:'log',dtick:1,autorange:true,gridcolor:C.grid},showlegend:false,uirevision:state.run+'-harmonics-'+(filtered?selected.ell+','+selected.m:'all')},ev=>{let g=Math.round(ev.points[0].x);if(selected&&filtered){let r=d.channels.find(x=>x.ell===selected.ell&&x.m===selected.m&&x.g===g);if(r&&r.status==='open')selectChannel(selected.ell,selected.m,g,true);else $('selectedChannel').textContent='This selected harmonic has no stored nonzero outgoing response.';}});
}
function renderCumulative(){let t=state.atlas.top;draw('cumulativePlot',[{type:'scatter',mode:'lines+markers',name:'Ranked cumulative fraction',x:t.map((_,i)=>i+1),y:t.map(p=>p.cumulative_fraction),line:{color:color[3],width:2},marker:{size:6},hovertemplate:'Top %{x} channels: %{y:.2%}<extra></extra>'}],
 {xaxis:{title:'Number of top-ranked outgoing channels',dtick:2,gridcolor:C.grid},yaxis:{title:'Fraction of all power',range:[0,1.05],tickformat:'.0%',gridcolor:C.grid},showlegend:false,margin:{l:62,r:17,t:15,b:50}});}
function fillSelector(){let options=state.atlas.channels.filter(x=>x.status==='open').sort((a,b)=>b.power-a.power);$('channelSelector').innerHTML=options.map(c=>`<option value="${c.ell},${c.m},${c.g}">(${c.ell}, ${c.m}, ${c.g}) · ${fmt(c.power)}</option>`).join('');}
async function loadChannel(ell,m,g){if(!state.run)return;let requested=[ell,m,g].join(',');$('channelMessage').textContent=`Loading (${requested})…`;try{let d=await getJSON(query('/api/channel',{ell,m,g}));if(!state.selected||state.selected.join(',')!==requested)return;state.channel=d;renderDeep();$('channelMessage').textContent='Reconstructed from saved primitives';}catch(e){$('channelMessage').textContent='Unable to reconstruct';showError('Channel analysis: '+e.message);}}
function deepLineChart(id,partials,total,title,yaxisLabel,phaseKey){let tr=[];let r=state.channel.radial.r;partials.forEach((p,i)=>tr.push({type:'scatter',mode:'lines',x:r,y:p.abs,name:'l='+p.l,line:{color:color[i%color.length],width:1.5}}));
 if(total)tr.push({type:'scatter',mode:'lines',x:r,y:total,name:'Coherent total',line:{color:'#ffffff',width:2.6}});
 if(phaseKey&&$('showPhase').checked)partials.forEach((p,i)=>tr.push({type:'scatter',mode:'lines',x:r,y:p[phaseKey],yaxis:'y2',name:'arg(l='+p.l+')',line:{color:color[i%color.length],width:1,dash:'dot'},opacity:.65}));
 let xx={title:'r / M',gridcolor:C.grid,type:$('radialAxis').value};
 let ly={title:yaxisLabel,type:'log',dtick:1,gridcolor:C.grid,autorange:true};let o={xaxis:xx,yaxis:ly,uirevision:state.run+'-deep-'+id,legend:{orientation:'h',y:-.24,x:0,font:{size:9}}};
 if(phaseKey&&$('showPhase').checked)o.yaxis2={overlaying:'y',side:'right',title:'phase (rad)',range:[-Math.PI,Math.PI],showgrid:false};
 draw(id,tr,o);}
function renderDeep(){let d=state.channel;if(!d)return;let sp=d.spectral,rep=d.response,rad=d.radial;
 $('channelKpis').innerHTML=card('Driven energy Δg',fmt(sp.Delta_g),'Open continuum')+
 card('KG momentum k',fmt(sp.k_kg),'Flux phase space')+
 card('Channel normalized power',fmt(rep.normalized_power),percent(rep.fraction)+' of total')+
 card('Multipole coherence C',fmt(rep.multipole_coherence),'1 aligned; ≪1 cancellation');
 deepLineChart('forcingPlot',rad.forcing,null,'Fourier forcing','|u_lmg(r)|','phase');
 deepLineChart('sourcePlot',rad.source_partials,rad.source_total_abs,'Partial scalar sources','|S(r)|','phase');
 deepLineChart('integrandPlot',rad.overlap_integrands,rad.overlap_total_abs,'Green overlap','|Z₋ S(r)|',null);
 let rr=rad.r,I=rep.I;
 let tr=[{type:'scatter',mode:'lines',name:'Re accumulated I',x:rr,y:rad.cumulative_real,line:{color:color[0],width:2}},
 {type:'scatter',mode:'lines',name:'Im accumulated I',x:rr,y:rad.cumulative_imag,line:{color:color[2],width:2}},
 {type:'scatter',mode:'lines',name:'Re production I',x:[rr[0],rr[rr.length-1]],y:[I.re,I.re],line:{color:color[0],dash:'dash'}},
 {type:'scatter',mode:'lines',name:'Im production I',x:[rr[0],rr[rr.length-1]],y:[I.im,I.im],line:{color:color[2],dash:'dash'}}];
 draw('cumulativeIPlot',tr,{xaxis:{title:'r / M',gridcolor:C.grid,type:$('radialAxis').value},yaxis:{title:'Complex overlap components',gridcolor:C.grid},legend:{orientation:'h',x:0,y:-.25,font:{size:9}},uirevision:state.run+'-cumulI'});
 let pairs=d.partial_overlaps;
 draw('overlapPlot',[{type:'bar',name:'|I_l|',x:pairs.map(p=>'l='+p.l),y:pairs.map(p=>p.abs),marker:{color:color[1]},hovertemplate:'%{x}<br>|I_l|=%{y:.4e}<extra></extra>'},
 {type:'scatter',mode:'markers+lines',name:'arg I_l',x:pairs.map(p=>'l='+p.l),y:pairs.map(p=>p.phase),yaxis:'y2',marker:{color:color[2],size:8}}],
 {xaxis:{title:'Forcing multipole',gridcolor:C.grid},yaxis:{title:'|I_l|',type:'log',dtick:1,gridcolor:C.grid},yaxis2:{overlaying:'y',side:'right',title:'Phase (rad)',range:[-Math.PI,Math.PI],showgrid:false},legend:{orientation:'h',y:-.2,x:0,font:{size:9}}});
 let inter=d.interference;
 draw('interferencePlot',[{type:'bar',x:inter.map(x=>`${x.l1} × ${x.l2}`),y:inter.map(x=>x.value),marker:{color:inter.map(x=>x.value>=0?color[3]:color[4])},name:'2 Re(I_l I_l′*)'}],
 {xaxis:{title:'Forcing-multipole pair',gridcolor:C.grid},yaxis:{title:'Signed interference',gridcolor:C.grid,zeroline:true,zerolinecolor:'#e7edf6'},showlegend:false});
 let table=[['Final labels',`(ℓf,mf,g)=(${d.label.ell},${d.label.m},${d.label.g})`],['Forcing m=mf−mi',d.label.forcing_m],['Forcing multipoles',pairs.map(p=>p.l).join(', ')],['Response energy Δg',fmt(sp.Delta_g)],['Carrier ωg',fmt(sp.omega_g)],['kNR',fmt(sp.k_nr)],['kKG',fmt(sp.k_kg)],['Coulomb η',fmt(sp.eta)],['Coherent overlap I',fmt(I.re)+' + i('+fmt(I.im)+')'],['Wronskian W',fmt(rep.W.re)+' + i('+fmt(rep.W.im)+')'],['Physical outgoing A∞',fmt(rep.A_out.re)+' + i('+fmt(rep.A_out.im)+')'],['Channel orbital power',fmt(rep.orbital_power)],['Share of total',percent(rep.fraction)],['Coherence factor',fmt(rep.multipole_coherence)],['Interference enhancement',fmt(rep.multipole_enhancement)],['Source reassembly relative mismatch',fmt(rep.source_relative_difference)],['Diagnostic grid I − saved I',fmt(rep.grid_minus_production_I.abs)]];
 $('responseTable').innerHTML=table.map(x=>kv(...x)).join('');}
function wire(){document.querySelectorAll('.sidebar button[data-page]').forEach(b=>b.onclick=()=>setPage(b.dataset.page));$('openAtlas').onclick=()=>setPage('channels');$('runSelector').onchange=ev=>loadRun(ev.target.value);$('reload').onclick=loadRuns;$('densityScale').onchange=renderDensity;$('extent').onchange=async()=>{try{let ext=$('extent').value;state.geometry=await getJSON(query('/api/geometry',ext==='auto'?{n:180}:{n:180,extent:ext}));renderDensity();}catch(e){showError(e.message);}};$('radialChoice').onchange=renderRadial;$('orbitalPhase').oninput=updatePhase;$('harmonicFilter').onchange=renderHarmonics;$('channelSelector').onchange=ev=>{if(!ev.target.value)return;selectChannel(...ev.target.value.split(',').map(Number),true);};$('showPhase').onchange=()=>{if(state.channel)renderDeep();};$('radialAxis').onchange=()=>{if(state.channel)renderDeep();};$('resetFilters').onclick=()=>{state.angular=null;renderHarmonics();$('selectedChannel').textContent='Showing all outgoing channels';};}
wire();loadRuns();
