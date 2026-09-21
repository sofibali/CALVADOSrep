import json, base64, os

D = "/tmp/claude-64170/-home-sbali-BindCraft/adecb958-709c-481a-8e03-db4ccc10a3b9/scratchpad"
fig = json.load(open(os.path.join(D, "figdata.json")))

def b64(p):
    return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()

IMG = {
    "block": b64(os.path.join(D, "report_md1_block_af3_crop.png")),
    "sim":   b64(os.path.join(D, "report_md1_block_sim_crop.png")),
    "clamp": b64(os.path.join(D, "report_clamp_md1md2_state3_crop.png")),
}
HITIMG={f"h{i}": b64(os.path.join(D, f"crop_hit_{i}.png")) for i in range(5)}
HITIMG["clamp"]=b64(os.path.join(D,"crop_clamp.png"))
HITIMG["variants"]=b64(os.path.join(D,"crop_variants.png"))
SEQG=json.load(open(os.path.join(D,"hits_seqs.json")))
import pandas as _pd
_ACC=_pd.read_csv("/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/designs/sweep_af3_guess/final_design_stats.csv")
_CACC=_pd.read_csv("/home/sbali/CALVADOS/examples/PARP14_MDP/bindcraft_md/designs/guess_clamp_md1md2_state3/final_design_stats.csv")

MD1HS = [33,37,39,41,42,43,44,46,47,134,135,136,138,172,173]

TARGETS = {
 "block": dict(key="md1_block_af3", name="md1_block_af3", concept="BLOCK", src="AF3 model",
   tagline="Occupy MD1's ADP-ribose pocket", hs={"MD1":MD1HS,"MD2":[]},
   kept="1–404 (MD1 + MD2 scaffold)", lengths="70–130 aa", dates="Aug 14 – Aug 28", days=14,
   early=155, clash=46, lowconf=146, relaxed=270,
   fail={"i_pAE":5370,"i_pTM":5243,"pLDDT":1283}, scored=10,
   rosetta={"n_InterfaceUnsatHbonds":10,"ShapeComplementarity":1,"n_InterfaceHbonds":1}, accepted=0),
 "sim": dict(key="md1_block_sim", name="md1_block_sim", concept="BLOCK", src="CALVADOS MD conformer",
   tagline="Same pocket, same hotspots — simulated conformer", hs={"MD1":MD1HS,"MD2":[]},
   kept="1–404 (MD1 + MD2 scaffold)", lengths="70–130 aa", dates="Sep 9 – Sep 16", days=7,
   early=170, clash=105, lowconf=108, relaxed=91,
   fail={"i_pAE":1820,"i_pTM":1816,"pLDDT":325}, scored=0, rosetta={}, accepted=0),
 "clamp": dict(key="clamp_md1md2_state3", name="clamp_md1md2_state3", concept="CLAMP", src="TICA state 3",
   tagline="Bridge MD1 and MD2 across the cleft",
   hs={"MD1":[71,102,103,106,107,109,110],"MD2":[351,354,357,395,396,399,400]},
   kept="1–404 (MD1 + MD2, state-3 pose)", lengths="90–150 aa", dates="Aug 28 – Sep 9", days=12,
   early=179, clash=110, lowconf=96, relaxed=220,
   fail={"i_pAE":4400,"i_pTM":4319,"pLDDT":1642}, scored=0, rosetta={}, accepted=0),
}
ORDER = ["block","sim","clamp"]
for t in TARGETS.values():
    t["total"] = t["early"]+t["clash"]+t["lowconf"]+t["relaxed"]
    t["seqs"]  = t["relaxed"]*20

C = {"md1":"var(--md1)","md2":"var(--md2)","sim":"var(--sim)","loss":"var(--loss)"}
SER = {"block":"var(--md1)","sim":"var(--sim)","clamp":"var(--md2)"}

def funnel_svg(t):
    stages=[("Aborted mid-hallucination",t["early"],"loss"),("Clashing",t["clash"],"loss"),
            ("Low final pLDDT",t["lowconf"],"loss"),("Relaxed — folded &amp; in contact",t["relaxed"],"md1")]
    W,rowh,gap,labelw=620,25,9,232
    H=len(stages)*(rowh+gap); tot=t["total"]
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="trajectory fate" class="chart">']
    y=0
    for lab,v,col in stages:
        bw=(v/tot)*(W-labelw-58) if tot else 0
        o.append(f'<text x="0" y="{y+rowh*0.72}" class="svg-lab">{lab}</text>')
        o.append(f'<rect x="{labelw}" y="{y}" width="{W-labelw-58}" height="{rowh}" class="track"/>')
        o.append(f'<rect x="{labelw}" y="{y}" width="{bw:.1f}" height="{rowh}" fill="{C[col]}"><title>{lab}: {v} of {tot} ({100*v/tot:.0f}%)</title></rect>')
        o.append(f'<text x="{W}" y="{y+rowh*0.72}" class="svg-num" text-anchor="end">{v}</text>')
        y+=rowh+gap
    o.append("</svg>"); return "".join(o)

def hist_svg(tk, groups):
    d=fig[TARGETS[tk]["key"]]["buckets"]["Relaxed"]
    W,H,pl,pb,pt=620,200,34,30,10
    bins=[round(0.5*i,1) for i in range(29)]
    ser=[(g,d[g]["hist"]) for g in groups if g in d]
    maxc=max(max(h) for _,h in ser) if ser else 1
    pw,ph=W-pl-8,H-pb-pt; bw=pw/28
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="closest approach" class="chart">']
    for f in (0,.25,.5,.75,1):
        yy=pt+ph*(1-f)
        o.append(f'<line x1="{pl}" y1="{yy:.1f}" x2="{W-8}" y2="{yy:.1f}" class="grid"/>')
        o.append(f'<text x="{pl-6}" y="{yy+3:.1f}" class="svg-tick" text-anchor="end">{int(maxc*f)}</text>')
    x4=pl+(4/14)*pw
    o.append(f'<line x1="{x4:.1f}" y1="{pt}" x2="{x4:.1f}" y2="{pt+ph}" class="ref"/>')
    o.append(f'<text x="{x4+4:.1f}" y="{pt+10}" class="svg-tick">4 Å</text>')
    for gi,(g,h) in enumerate(ser):
        col=C["md1"] if g=="MD1" else C["md2"]
        if tk=="sim": col=C["sim"]
        for i,c in enumerate(h):
            if not c: continue
            bh=(c/maxc)*ph; x=pl+i*bw; y=pt+ph-bh
            o.append(f'<rect x="{x+0.6:.1f}" y="{y:.1f}" width="{bw-1.2:.1f}" height="{bh:.1f}" fill="{col}" opacity="{0.9 if gi==0 else 0.6}"><title>{g}: {c} at {bins[i]}–{bins[i]+0.5} Å</title></rect>')
    for xv in (0,2,4,6,8,10,12,14):
        o.append(f'<text x="{pl+(xv/14)*pw:.1f}" y="{H-12}" class="svg-tick" text-anchor="middle">{xv}</text>')
    o.append(f'<text x="{pl+pw/2:.1f}" y="{H-1}" class="svg-tick" text-anchor="middle">closest binder→hotspot distance (Å)</text>')
    o.append("</svg>"); return "".join(o)

def scatter_svg():
    pts=fig["clamp_md1md2_state3"]["points"]
    W,H,pl,pb,pt,pr=620,420,46,42,14,14
    pw,ph=W-pl-pr,H-pb-pt; MX=14
    sx=lambda v: pl+min(v,MX)/MX*pw
    sy=lambda v: pt+ph-min(v,MX)/MX*ph
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="MD1 vs MD2 approach" class="chart">']
    o.append(f'<rect x="{sx(0)}" y="{sy(4):.1f}" width="{sx(4)-sx(0):.1f}" height="{sy(0)-sy(4):.1f}" class="quad"/>')
    o.append(f'<text x="{sx(4)+6:.1f}" y="{sy(0)-8:.1f}" class="svg-tick">both &lt; 4 Å — clamp geometry satisfied</text>')
    for v in (0,2,4,6,8,10,12,14):
        o.append(f'<line x1="{sx(v):.1f}" y1="{pt}" x2="{sx(v):.1f}" y2="{pt+ph}" class="grid"/>')
        o.append(f'<line x1="{pl}" y1="{sy(v):.1f}" x2="{W-pr}" y2="{sy(v):.1f}" class="grid"/>')
        o.append(f'<text x="{sx(v):.1f}" y="{H-24}" class="svg-tick" text-anchor="middle">{v}</text>')
        o.append(f'<text x="{pl-7}" y="{sy(v)+3:.1f}" class="svg-tick" text-anchor="end">{v}</text>')
    bc={"R":C["md1"],"C":C["loss"],"L":C["md2"]}
    for x,y,b in pts:
        o.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.4" fill="{bc[b]}" opacity="0.6" class="dot"><title>MD1 {x} Å · MD2 {y} Å</title></circle>')
    o.append(f'<text x="{pl+pw/2:.1f}" y="{H-6}" class="svg-tick" text-anchor="middle">closest approach to MD1 hotspots (Å)</text>')
    o.append(f'<text transform="translate(11,{pt+ph/2:.1f}) rotate(-90)" class="svg-tick" text-anchor="middle">closest approach to MD2 hotspots (Å)</text>')
    o.append("</svg>"); return "".join(o)

def fate_ab_svg():
    W,rowh,gap,labelw=620,30,16,150
    rows=[("md1_block_af3 (AF3)","block"),("md1_block_sim (MD)","sim")]
    H=len(rows)*(rowh+gap)
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="fate af3 vs sim" class="chart">']
    y=0; bwtot=W-labelw-8
    for lab,k in rows:
        t=TARGETS[k]; tot=t["total"]; x=labelw
        segs=[("aborted",t["early"],"var(--loss)",.45),("clashing",t["clash"],"var(--flag)",1),
              ("low pLDDT",t["lowconf"],"var(--loss)",.75),("relaxed",t["relaxed"],SER[k],1)]
        o.append(f'<text x="0" y="{y+rowh*0.68}" class="svg-lab">{lab}</text>')
        for nm,v,col,op in segs:
            sw=v/tot*bwtot
            o.append(f'<rect x="{x:.1f}" y="{y}" width="{max(sw-2,0):.1f}" height="{rowh}" fill="{col}" opacity="{op}"><title>{nm}: {v} of {tot} ({100*v/tot:.0f}%)</title></rect>')
            if sw>46:
                o.append(f'<text x="{x+sw/2:.1f}" y="{y+rowh*0.66}" class="svg-inbar" text-anchor="middle">{100*v/tot:.0f}%</text>')
            x+=sw
        y+=rowh+gap
    o.append("</svg>"); return "".join(o)

def failnorm_svg():
    rows=[("Interface PAE too high","i_pAE ≤ 0.35"),("Interface pTM too low","i_pTM ≥ 0.50"),("Binder pLDDT too low","pLDDT ≥ 0.80")]
    keys=["i_pAE","i_pTM","pLDDT"]
    W,rowh,gap,labelw=620,19,30,200
    H=len(rows)*(rowh*3+gap)
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="MPNN sequence rejection rate" class="chart">']
    y=0
    for (lab,thr),key in zip(rows,keys):
        o.append(f'<text x="0" y="{y+rowh*0.72}" class="svg-lab">{lab}</text>')
        o.append(f'<text x="0" y="{y+rowh*1.72}" class="svg-sub">{thr}</text>')
        for i,k in enumerate(ORDER):
            t=TARGETS[k]; pct=100*t["fail"][key]/t["seqs"]
            bw=pct/100*(W-labelw-74); yy=y+i*rowh
            o.append(f'<rect x="{labelw}" y="{yy+2}" width="{W-labelw-74}" height="{rowh-5}" class="track"/>')
            o.append(f'<rect x="{labelw}" y="{yy+2}" width="{bw:.1f}" height="{rowh-5}" fill="{SER[k]}"><title>{t["name"]}: {t["fail"][key]:,} of {t["seqs"]:,} sequences ({pct:.1f}%)</title></rect>')
            o.append(f'<text x="{W}" y="{yy+rowh*0.72}" class="svg-num" text-anchor="end">{pct:.1f}%</text>')
        y+=rowh*3+gap
    o.append("</svg>"); return "".join(o)


def replay_svg():
    r=fig["_replay"]; W,H,pl,pb,pt=620,230,34,34,12
    pw,ph=W-pl-8,H-pb-pt; n=len(r["bins"]); bw=pw/n
    maxc=max(max(r["off"]),max(r["on"]))
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="replay i_pAE before and after" class="chart">']
    for f in (0,.5,1):
        yy=pt+ph*(1-f)
        o.append(f'<line x1="{pl}" y1="{yy:.1f}" x2="{W-8}" y2="{yy:.1f}" class="grid"/>')
        o.append(f'<text x="{pl-6}" y="{yy+3:.1f}" class="svg-tick" text-anchor="end">{int(maxc*f)}</text>')
    xt=pl+0.35*pw
    o.append(f'<line x1="{xt:.1f}" y1="{pt}" x2="{xt:.1f}" y2="{pt+ph}" class="ref"/>')
    o.append(f'<text x="{xt+4:.1f}" y="{pt+10}" class="svg-tick">threshold 0.35</text>')
    for key,col,op in (("on","var(--win)",0.95),("off","var(--md2)",0.75)):
        for i,c in enumerate(r[key]):
            if not c: continue
            bh=(c/maxc)*ph; x=pl+i*bw
            lbl="with guess" if key=="on" else "without"
            o.append(f'<rect x="{x+0.5:.1f}" y="{pt+ph-bh:.1f}" width="{bw-1:.1f}" height="{bh:.1f}" fill="{col}" opacity="{op}"><title>{lbl}: {c} predictions at i_pAE {r["bins"][i]:.2f}</title></rect>')
    for xv in (0,.2,.4,.6,.8,1.0):
        o.append(f'<text x="{pl+xv*pw:.1f}" y="{H-14}" class="svg-tick" text-anchor="middle">{xv:.1f}</text>')
    o.append(f'<text x="{pl+pw/2:.1f}" y="{H-2}" class="svg-tick" text-anchor="middle">interface PAE (i_pAE) — lower is better</text>')
    o.append("</svg>"); return "".join(o)

def armA_svg():
    A=fig["_armA"]
    rows=[("MPNN sequences rejected by i_pAE", 99.4, 100*A["fail"]["i_pAE"]/A["seqs"], "lower better"),
          ("Sequences reaching Rosetta scoring", 100*10/5400, 100*A["scored"]/A["seqs"], "higher better")]
    W,rowh,gap,labelw=620,20,30,236
    H=len(rows)*(rowh*2+gap)
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="baseline vs arm A" class="chart">']
    y=0
    for lab,base,arm,hint in rows:
        o.append(f'<text x="0" y="{y+rowh*0.72}" class="svg-lab">{lab}</text>')
        o.append(f'<text x="0" y="{y+rowh*1.72}" class="svg-sub">{hint}</text>')
        for i,(v,col,nm) in enumerate([(base,"var(--md1)","baseline"),(arm,"var(--win)","with guess")]):
            bwid=max(v,0.4)/100*(W-labelw-70); yy=y+i*rowh
            o.append(f'<rect x="{labelw}" y="{yy+2}" width="{W-labelw-70}" height="{rowh-5}" class="track"/>')
            o.append(f'<rect x="{labelw}" y="{yy+2}" width="{bwid:.1f}" height="{rowh-5}" fill="{col}"><title>{nm}: {v:.1f}%</title></rect>')
            o.append(f'<text x="{W}" y="{yy+rowh*0.72}" class="svg-num" text-anchor="end">{v:.1f}%</text>')
        y+=rowh*2+gap
    o.append("</svg>"); return "".join(o)



def gates_svg():
    g=fig["_queue"]["gates3"]
    order=["n_InterfaceUnsatHbonds","ShapeComplementarity","Surface_Hydrophobicity",
           "Binder_RMSD","n_InterfaceHbonds","dG"]
    W,rowh,gap,labelw=620,22,12,238
    H=len(order)*(rowh+gap)
    o=[f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Rosetta gate failure rates" class="chart">']
    y=0
    for k in order:
        v=g[k]; pct=v["pct"]
        bw=max(pct,0.5)/100*(W-labelw-64)
        col="var(--flag)" if pct>50 else ("var(--md2)" if pct>10 else "var(--win)")
        o.append(f'<text x="0" y="{y+rowh*0.72}" class="svg-lab">{k}</text>')
        o.append(f'<rect x="{labelw}" y="{y+2}" width="{W-labelw-64}" height="{rowh-4}" class="track"/>')
        o.append(f'<rect x="{labelw}" y="{y+2}" width="{bw:.1f}" height="{rowh-4}" fill="{col}"><title>{k}: {v["fail"]} of {v["n"]} fail (threshold {v["thr"]}, median {v["med"]})</title></rect>')
        o.append(f'<text x="{W}" y="{y+rowh*0.72}" class="svg-num" text-anchor="end">{pct:.0f}%</text>')
        y+=rowh+gap
    o.append("</svg>"); return "".join(o)

def hits_section():
    import re as _re
    order=["l80_s632138","l130_s735229","l128_s514541","l122_s135691","l104_s444809"]
    cards=[]
    for i,traj in enumerate(order):
        sub=_ACC[_ACC["Design"].str.contains(traj)]
        best=sub.sort_values("Average_i_pAE").iloc[0]
        cards.append(f'''<figure class="hit">
      <img src="{HITIMG[f"h{i}"]}" alt="binder from trajectory {traj} bound to MD1"/>
      <figcaption><span class="hname" style="color:{["var(--md1)","var(--md2)","var(--sim)","#eda100","#e87ba4"][i]}">&#9632;</span>
      <code>{traj}</code> · {int(best["Length"])} aa · {len(sub)} accepted<br>
      i_pAE {best["Average_i_pAE"]:.2f} · SC {best["Average_ShapeComplementarity"]:.2f} · unsat {best["Average_n_InterfaceUnsatHbonds"]:.1f} · dG {best["Average_dG"]:.0f}</figcaption>
    </figure>''')
    # sequence alignment blocks
    rows=[]
    for g in SEQG:
        seqs=[x["seq"] for x in g["seqs"]]
        if len(seqs)>1:
            diff=[a!=b for a,b in zip(seqs[0],seqs[1])]
        else:
            diff=[False]*len(seqs[0])
        lines=[]
        for x in g["seqs"]:
            spans="".join(f'<span class="{"mm" if d else ""}">{c}</span>' for c,d in zip(x["seq"],diff))
            lines.append(f'<div class="seqrow"><span class="seqid">{x["name"]}</span><span class="seq">{spans}</span></div>')
        ident=f'{g["identity"]}% identical' if g["identity"] else "single accepted design"
        rows.append(f'''<div class="seqblock">
      <div class="seqhead"><code>{g["traj"]}</code><span>{g["label"]} · {g["len"]} aa · {ident}</span></div>
      {"".join(lines)}</div>''')
    cb=_CACC.sort_values("Average_i_pAE").iloc[0]
    return f'''
<h2>What the hits actually look like</h2>
<p class="sub">Every accepted MD1 binder, one panel per source backbone, all in the same orientation. Target in
grey, the 15 hotspot side chains in violet, the designed binder coloured.</p>
<div class="hitgrid">{"".join(cards)}</div>
<p class="cap">Five independent trajectories, five different folds and lengths (80–130 aa) — and all five land on
the <b>same face</b>, packed against the hotspot patch. That convergence is the design working as intended:
nothing constrained these backbones to agree, only the hotspot definition did.</p>

<h3 style="margin-top:2.4rem">The clamp bridges both domains</h3>
<div class="panel" style="margin-top:.8rem">
  <img class="wide" src="{HITIMG["clamp"]}" alt="clamp binder bridging MD1 and MD2"/>
  <p class="cap">Best clamp design (<code>l97_s990418</code>, 97 aa, i_pAE {cb["Average_i_pAE"]:.2f},
  SC {cb["Average_ShapeComplementarity"]:.2f}, dG {cb["Average_dG"]:.0f}). MD1 in cool grey, MD2 in warm grey,
  binder in orange; contact residues shown as sticks — <span class="k1">7 on MD1</span> and
  <span style="color:var(--sim);font-weight:600">18 on MD2</span>. It lies across the cleft touching both,
  which is the whole point of a clamp and the thing a passing filter score alone would not prove. That run has since finished with <b>5 accepted designs</b> from 50 trajectories.</p>
</div>

<h3 style="margin-top:2.4rem">Same backbone, different sequences</h3>
<div class="panel" style="margin-top:.8rem">
  <div class="twocol">
    <img src="{HITIMG["variants"]}" alt="two accepted MPNN variants on one backbone"/>
    <div>
      <p style="margin-top:0">Where a backbone yields more than one accepted design, MPNN found genuinely
      different sequences for it — <b>76–86% identical</b>, so 14–24% of positions differ — and both still fold
      and dock. Structurally they are near-superimposable; the variation is chemical, not conformational.</p>
      <p style="margin-bottom:0">That matters for ordering: each backbone gives several independent sequence
      attempts at the same binding mode, rather than one fragile candidate.</p>
    </div>
  </div>
  <div class="seqwrap">{"".join(rows)}</div>
  <p class="cap">Accepted sequences grouped by source backbone. Highlighted positions differ between the two
  accepted variants of that backbone.</p>
</div>
'''

def hs_chips(t):
    h=[]
    for g in ("MD1","MD2"):
        if not t["hs"][g]: continue
        cls="c1" if g=="MD1" else "c2"
        h.append(f'<div class="hsrow"><span class="hstag {cls}">{g}</span><span class="hsnums">{" ".join(str(x) for x in t["hs"][g])}</span><span class="hsn">{len(t["hs"][g])}</span></div>')
    return "".join(h)

cards=[]
for k in ORDER:
    t=TARGETS[k]
    groups=["MD1"]+(["MD2"] if t["hs"]["MD2"] else [])
    d=fig[t["key"]]["buckets"]["Relaxed"]
    stat=" · ".join(f'{g} median <b>{d[g]["med"]} Å</b>' for g in groups if g in d)
    badge="b1" if k=="block" else ("b3" if k=="sim" else "b2")
    cards.append(f"""
<article class="card">
  <header class="card-h">
    <span class="badge {badge}">{t['concept']}</span>
    <h3>{t['name']}</h3>
    <p class="tag">{t['tagline']}</p>
  </header>
  <figure class="struct"><img src="{IMG[k]}" alt="{t['name']} target with hotspots"/>
  <figcaption>Hotspot side chains as sticks — <span class="k1">MD1</span>{' and <span class="k2">MD2</span>' if t['hs']['MD2'] else ''}. Receptor source: {t['src']}.</figcaption></figure>
  <dl class="spec">
    <div><dt>Receptor</dt><dd>{t['src']}</dd></div>
    <div><dt>Residues kept</dt><dd>{t['kept']}</dd></div>
    <div><dt>Binder length</dt><dd>{t['lengths']}</dd></div>
    <div><dt>Ran</dt><dd>{t['dates']} (~{t['days']} d)</dd></div>
  </dl>
  <div class="hs">{hs_chips(t)}</div>
  <h4>Where the {t['total']:,} attempts went</h4>
  {funnel_svg(t)}
  <h4>How close the folded binders got</h4>
  {hist_svg(k,groups)}
  <p class="cap">Relaxed only (n={t['relaxed']}). {stat}.</p>
</article>""")

def trow(label,f,cls=""):
    return "<tr><td>"+label+"</td>"+"".join(f'<td class="n {cls if f(TARGETS[k]) in ("0","—") and cls else ""}">{f(TARGETS[k])}</td>' for k in ORDER)+"</tr>"

html=f"""<title>Three Misses And A Fix</title>
<style>
:root {{
  --paper:#f5f7fa; --panel:#ffffff; --ink:#111418; --ink2:#4a5463; --muted:#8a93a3;
  --line:rgba(17,20,24,.12); --grid:rgba(17,20,24,.09);
  --md1:#2a78d6; --md2:#eb6834; --sim:#1baf7a; --loss:#96a0b0; --quad:rgba(42,120,214,.09);
  --flag:#b8442a; --flagbg:#fdeee8; --win:#4a3aa7; --winbg:#eceafa;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --paper:#0f1216; --panel:#171b21; --ink:#eef2f7; --ink2:#aab4c2; --muted:#7f8a99;
  --line:rgba(255,255,255,.14); --grid:rgba(255,255,255,.10);
  --md1:#5b9bec; --md2:#f5854f; --sim:#2cc98d; --loss:#5d6773; --quad:rgba(91,155,236,.13);
  --flag:#ff9a76; --flagbg:#2e1a13; --win:#9085e9; --winbg:#221d3d; }} }}
:root[data-theme="dark"] {{
  --paper:#0f1216; --panel:#171b21; --ink:#eef2f7; --ink2:#aab4c2; --muted:#7f8a99;
  --line:rgba(255,255,255,.14); --grid:rgba(255,255,255,.10);
  --md1:#5b9bec; --md2:#f5854f; --sim:#2cc98d; --loss:#5d6773; --quad:rgba(91,155,236,.13);
  --flag:#ff9a76; --flagbg:#2e1a13; --win:#9085e9; --winbg:#221d3d; }}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--paper);color:var(--ink);
  font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;font-size:16px;line-height:1.55;-webkit-font-smoothing:antialiased}}
.eyebrow,.svg-lab,.svg-num,.svg-tick,.svg-sub,.svg-inbar,.hstag,.hsnums,.badge,dt,.hsn,th,code
  {{font-family:ui-monospace,"SF Mono","Cascadia Code","JetBrains Mono",Menlo,Consolas,monospace}}
.wrap{{max-width:1180px;margin:0 auto;padding:clamp(28px,5vw,64px) clamp(18px,4vw,40px) 90px}}
.eyebrow{{font-size:.72rem;letter-spacing:.17em;text-transform:uppercase;color:var(--md1);display:flex;align-items:center;gap:.6em;margin:0 0 1.1rem}}
.eyebrow::before{{content:"";width:1.7em;height:1px;background:var(--md1)}}
h1{{font-size:clamp(1.9rem,3.9vw,2.9rem);line-height:1.1;letter-spacing:-.018em;margin:0 0 .7rem;text-wrap:balance;font-weight:650}}
.dek{{font-size:1.09rem;color:var(--ink2);max-width:68ch;margin:0 0 2.2rem}}
h2{{font-size:1.32rem;letter-spacing:-.01em;margin:3.4rem 0 .5rem;font-weight:640}}
h2:first-of-type{{margin-top:2.6rem}}
.sub{{color:var(--ink2);margin:0 0 1.4rem;max-width:74ch;font-size:.97rem}}
h3{{font-size:1.04rem;margin:.1rem 0 .15rem;font-weight:640}}
h4{{font-size:.74rem;letter-spacing:.11em;text-transform:uppercase;color:var(--muted);margin:1.8rem 0 .65rem;font-weight:600;font-family:ui-monospace,Menlo,monospace}}
.grid3{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:5px;padding:20px 20px 24px;min-width:0}}
.card-h{{display:flex;flex-wrap:wrap;align-items:baseline;gap:.5rem .65rem;margin-bottom:.85rem}}
.badge{{font-size:.64rem;letter-spacing:.13em;padding:.28em .6em;border-radius:3px;font-weight:600}}
.b1{{background:color-mix(in srgb,var(--md1) 15%,transparent);color:var(--md1)}}
.b2{{background:color-mix(in srgb,var(--md2) 17%,transparent);color:var(--md2)}}
.b3{{background:color-mix(in srgb,var(--sim) 17%,transparent);color:var(--sim)}}
.tag{{flex-basis:100%;margin:0;color:var(--ink2);font-size:.9rem}}
.struct{{margin:0 0 1rem;background:#fff;border:1px solid var(--line);border-radius:4px;overflow:hidden}}
.struct img{{display:block;width:100%;height:auto}}
figcaption{{font-size:.76rem;color:var(--ink2);padding:.5rem .65rem;background:var(--panel);border-top:1px solid var(--line)}}
.k1{{color:var(--md1);font-weight:600}} .k2{{color:var(--md2);font-weight:600}}
.spec{{display:grid;grid-template-columns:1fr 1fr;gap:1px;background:var(--line);border:1px solid var(--line);margin:0 0 .9rem}}
.spec>div{{background:var(--panel);padding:.55rem .7rem}}
dt{{font-size:.62rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);margin-bottom:.2rem}}
dd{{margin:0;font-size:.86rem}}
.hs{{display:flex;flex-direction:column;gap:6px}}
.hsrow{{display:flex;align-items:baseline;gap:.55rem;padding:.4rem .55rem;border:1px solid var(--line);border-radius:3px}}
.hstag{{font-size:.63rem;letter-spacing:.09em;padding:.15em .42em;border-radius:2px;flex:none}}
.c1{{background:color-mix(in srgb,var(--md1) 15%,transparent);color:var(--md1)}}
.c2{{background:color-mix(in srgb,var(--md2) 17%,transparent);color:var(--md2)}}
.hsnums{{font-size:.74rem;color:var(--ink2);line-height:1.45;flex:1;min-width:0;word-spacing:.15em}}
.hsn{{font-size:.74rem;color:var(--muted);flex:none;font-variant-numeric:tabular-nums}}
.chart{{width:100%;height:auto;display:block;overflow:visible}}
.track{{fill:var(--grid)}} .grid{{stroke:var(--grid);stroke-width:1}}
.ref{{stroke:var(--md2);stroke-width:1.5;stroke-dasharray:3 3;opacity:.75}}
.quad{{fill:var(--quad)}} .dot{{stroke:var(--panel);stroke-width:.7}}
.svg-lab{{font-size:12.5px;fill:var(--ink)}} .svg-sub{{font-size:10.5px;fill:var(--muted)}}
.svg-num{{font-size:12.5px;fill:var(--ink);font-variant-numeric:tabular-nums}}
.svg-tick{{font-size:10.5px;fill:var(--muted);font-variant-numeric:tabular-nums}}
.svg-inbar{{font-size:10.5px;fill:#fff;font-variant-numeric:tabular-nums;font-weight:600}}
.cap{{font-size:.82rem;color:var(--ink2);margin:.7rem 0 0}}
.legend{{display:flex;flex-wrap:wrap;gap:1.1rem;margin:.2rem 0 1rem;font-size:.82rem;color:var(--ink2)}}
.legend span{{display:inline-flex;align-items:center;gap:.42em}}
.sw{{width:11px;height:11px;border-radius:2px;flex:none}}
.panel{{background:var(--panel);border:1px solid var(--line);border-radius:5px;padding:22px}}
.flag{{background:var(--flagbg);border:1px solid color-mix(in srgb,var(--flag) 35%,transparent);border-left:3px solid var(--flag);border-radius:4px;padding:1rem 1.15rem;margin:1.3rem 0 0}}
.flag b{{color:var(--flag)}}
.win{{background:var(--winbg);border:1px solid color-mix(in srgb,var(--win) 35%,transparent);
  border-left:3px solid var(--win);border-radius:4px;padding:1rem 1.15rem;margin:1.3rem 0 0}}
.win b{{color:var(--win)}}
.win .lab{{font-family:ui-monospace,Menlo,monospace;font-size:.65rem;letter-spacing:.11em;
  text-transform:uppercase;color:var(--win);display:block;margin-bottom:.4rem}}
td.good{{color:var(--win);font-weight:600}}
.hitgrid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin-top:.4rem}}
.hit{{margin:0;background:var(--panel);border:1px solid var(--line);border-radius:5px;overflow:hidden}}
.hit img{{display:block;width:100%;height:auto;background:#fff}}
.hit figcaption{{font-size:.73rem;line-height:1.5;color:var(--ink2);padding:.5rem .6rem;border-top:1px solid var(--line)}}
.hname{{margin-right:.35em}}
img.wide{{display:block;width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:4px}}
.twocol{{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1.1fr);gap:20px;align-items:start}}
.twocol img{{width:100%;height:auto;background:#fff;border:1px solid var(--line);border-radius:4px}}
.seqwrap{{margin-top:1.3rem;display:flex;flex-direction:column;gap:.9rem}}
.seqblock{{border:1px solid var(--line);border-radius:4px;padding:.6rem .7rem;overflow-x:auto}}
.seqhead{{display:flex;gap:.7rem;align-items:baseline;font-size:.74rem;color:var(--muted);margin-bottom:.45rem}}
.seqhead code{{color:var(--ink);background:none;padding:0}}
.seqrow{{display:flex;gap:.6rem;align-items:baseline;white-space:nowrap}}
.seqid{{font-family:ui-monospace,Menlo,monospace;font-size:.68rem;color:var(--muted);width:4.2em;flex:none}}
.seq{{font-family:ui-monospace,Menlo,monospace;font-size:.7rem;letter-spacing:.02em;color:var(--ink2)}}
.seq .mm{{background:color-mix(in srgb,var(--win) 26%,transparent);color:var(--ink);border-radius:2px}}
@media (max-width:640px){{.twocol{{grid-template-columns:1fr}}}}
.flag .lab{{font-family:ui-monospace,Menlo,monospace;font-size:.65rem;letter-spacing:.11em;text-transform:uppercase;color:var(--flag);display:block;margin-bottom:.4rem}}
code{{font-size:.86em;background:var(--grid);padding:.1em .35em;border-radius:3px}}
.tblwrap{{overflow-x:auto;margin-top:1rem}}
table{{width:100%;border-collapse:collapse;font-size:.9rem;min-width:640px}}
th{{text-align:left;font-size:.64rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);font-weight:600;padding:0 .8rem .55rem 0;border-bottom:1px solid var(--line)}}
td{{padding:.62rem .8rem .62rem 0;border-bottom:1px solid var(--line);vertical-align:baseline}}
td.n{{font-family:ui-monospace,Menlo,monospace;font-variant-numeric:tabular-nums}}
.zero{{color:var(--flag);font-weight:600}}
footer{{margin-top:3.4rem;padding-top:1.2rem;border-top:1px solid var(--line);font-size:.8rem;color:var(--muted)}}
@media (max-width:560px){{.spec{{grid-template-columns:1fr}}}}
</style>

<div class="wrap">
<p class="eyebrow">PARP14 · BindCraft diagnostic</p>
<h1>Three misses and a fix</h1>
<p class="dek">Three PARP14 binder campaigns spent ~33 GPU-days and produced <b>zero</b> accepted designs.
Tracing every attempt from hallucination to final filter found one gate responsible for almost all of it —
and changing a single setting turned the same target into <b>9 accepted designs in 17 hours</b>.
The failure analysis is kept below in full, because it is what located the fix.</p>

<h2>The three runs</h2>
<p class="sub">Identical pipeline, identical filters. Two of them (<code>md1_block_af3</code> /
<code>md1_block_sim</code>) are the <em>same pocket and the same 15 hotspots</em> — they differ only in which
structure of MD1 was handed to BindCraft. The third asks for a two-domain clamp instead.</p>
<div class="grid3">{cards[0]}{cards[1]}{cards[2]}</div>

<h2>The receptor conformer matters enormously</h2>
<p class="sub">Same pocket, same hotspots, same binder lengths — only the input structure changes. The AF3
model is a clean idealised prediction; the “sim” target is a CALVADOS coarse-grained frame back-mapped to
all-atom and rigid-body declashed. That difference alone reshapes the run.</p>
<div class="panel">
  <div class="legend">
    <span><i class="sw" style="background:var(--loss);opacity:.45"></i>aborted</span>
    <span><i class="sw" style="background:var(--flag)"></i>clashing</span>
    <span><i class="sw" style="background:var(--loss);opacity:.75"></i>low pLDDT</span>
    <span><i class="sw" style="background:var(--md1)"></i>relaxed (af3)</span>
    <span><i class="sw" style="background:var(--sim)"></i>relaxed (sim)</span>
  </div>
  {fate_ab_svg()}
  <p class="cap">Trajectory fate as a share of all attempts. On the simulated conformer the clash rate
  jumps from <b>10% to 35%</b> and the yield of usable folded trajectories falls from <b>58% to 30%</b>.</p>
  <div class="flag">
    <span class="lab">This was predicted</span>
    The project's own target-prep notes flagged that CG-backmapped multi-domain poses “clash at atomic
    detail” and warned to <b>flag states giving uniformly bad poses</b>. That is exactly what the numbers
    show — the declashing step removes gross interpenetration but leaves a surface AF2 finds markedly
    harder to dock against.
  </div>
</div>

<h2>The clamp geometry is working</h2>
<p class="sub">The obvious hypothesis for the clamp's total failure was that spanning two domains at a fixed
pose is geometrically too hard — that binders never manage to touch both sides. The coordinates say otherwise.</p>
<div class="panel">
  <div class="legend">
    <span><i class="sw" style="background:var(--md1)"></i>Relaxed (folded OK)</span>
    <span><i class="sw" style="background:var(--loss)"></i>Clashing</span>
    <span><i class="sw" style="background:var(--md2)"></i>Low confidence</span>
  </div>
  {scatter_svg()}
  <p class="cap">All 426 clamp trajectories, by closest approach to each domain's hotspot patch. Of the 220
  that folded cleanly, <b>206 (94%) contact both MD1 and MD2 within 4 Å</b> — the shaded quadrant. The
  binders land exactly where the design asks.</p>
</div>

<h2>One filter kills almost everything</h2>
<p class="sub">After a trajectory folds, ProteinMPNN rewrites its sequence (20 per trajectory) and AlphaFold2
re-predicts the complex. Below: the share of those sequences rejected by each AF2 gate, normalised so the
three runs are comparable.</p>
<div class="panel">
  <div class="legend">
    <span><i class="sw" style="background:var(--md1)"></i>md1_block_af3</span>
    <span><i class="sw" style="background:var(--sim)"></i>md1_block_sim</span>
    <span><i class="sw" style="background:var(--md2)"></i>clamp_md1md2_state3</span>
  </div>
  {failnorm_svg()}
  <p class="cap">Rejections as a percentage of MPNN sequences generated
  (5,400 / 1,820 / 4,400 respectively).</p>
  <div class="flag">
    <span class="lab">The headline number</span>
    <b>Interface PAE rejects 100.0% of MPNN sequences in both md1_block_sim and clamp_md1md2_state3</b>
    — 1,820 of 1,820 and 4,400 of 4,400, not one survivor. In md1_block_af3 it rejects 99.4%, and those
    ~30 survivors are the entire reason that run has any scored designs at all. Whatever else differs
    between these targets, <code>i_pAE ≤ 0.35</code> is the common wall.
  </div>
</div>

<h2>Where each run stops</h2>
<p class="sub">All three end at zero, but they are stopped at different stages — which is what determines
what to change next.</p>
<div class="tblwrap">
<table>
<thead><tr><th>Stage</th><th>md1_block_af3</th><th>md1_block_sim</th><th>clamp_md1md2_state3</th></tr></thead>
<tbody>
<tr><td>Hallucination attempts</td><td class="n">617</td><td class="n">474</td><td class="n">605</td></tr>
<tr><td>Folded &amp; in contact (relaxed)</td><td class="n">270</td><td class="n">91</td><td class="n">220</td></tr>
<tr><td>Hotspot patch engaged</td><td class="n">270 (100%)</td><td class="n">91 (100%)</td><td class="n">206 both-domain (94%)</td></tr>
<tr><td>MPNN sequences generated</td><td class="n">5,400</td><td class="n">1,820</td><td class="n">4,400</td></tr>
<tr><td>Survived AF2 re-prediction</td><td class="n">~30</td><td class="n zero">0</td><td class="n zero">0</td></tr>
<tr><td>Reached Rosetta scoring</td><td class="n">10</td><td class="n zero">0</td><td class="n zero">0</td></tr>
<tr><td>Passed Rosetta filters</td><td class="n zero">0</td><td class="n">—</td><td class="n">—</td></tr>
<tr><td>Accepted</td><td class="n zero">0</td><td class="n zero">0</td><td class="n zero">0</td></tr>
</tbody>
</table>
</div>

<div class="panel" style="margin-top:20px">
  <h4 style="margin-top:0">What this means</h4>
  <p style="margin:0 0 .8rem"><b>md1_block_af3</b> is the only run that reaches the final biophysical check.
  It loses all 10 designs on one metric — buried unsatisfied H-bonds (5–9 against a limit of 4) — while AF2
  confidence passes comfortably. Geometrically real interface, chemically unsatisfying.</p>
  <p style="margin:0 0 .8rem"><b>md1_block_sim</b> never gets there. The back-mapped conformer costs it
  trajectories up front (35% clash rate) and then the surviving 91 produce 1,820 sequences that AF2 rejects
  without exception. Same pocket as af3 — worse structure to design against.</p>
  <p style="margin:0"><b>clamp_md1md2_state3</b> places binders correctly across both domains, then loses
  every sequence at the same i_pAE wall. Its bottleneck is sequence–structure agreement at a two-domain
  interface, not geometry — so loosening hotspots would not help.</p>
</div>

<h2>The fix: one setting</h2>
<p class="sub">If the failures had been near-misses, loosening the threshold would have been defensible. They were
not. So the question became whether AF2 could recover the binding mode at all given a little help. BindCraft has a
flag for exactly this — <code>predict_initial_guess</code> — which hands the validation model the binder's own
trajectory coordinates instead of asking it to fold and dock from sequence alone.</p>
<div class="panel">
  <div class="legend">
    <span><i class="sw" style="background:var(--md2)"></i>without guess (reproduces the failed runs)</span>
    <span><i class="sw" style="background:var(--win)"></i>with predict_initial_guess</span>
  </div>
  {replay_svg()}
  <p class="cap">A <b>paired</b> replay over already-generated trajectories: identical backbone, identical MPNN
  sequence, identical model — only the validation protocol differs (n=396 pairs). Median i_pAE moves
  <b>0.877 → 0.508</b>; 96.5% of pairs improve; <b>118 rescued, 0 lost</b>. Reusing banked trajectories meant this
  took hours rather than the days a fresh campaign would have needed.</p>
  <div class="flag">
    <span class="lab">The honest caveat</span>
    This does test a weaker criterion — "is this sequence consistent with this pose?" rather than "does this
    sequence find this pose?". But it is <b>not</b> a rubber stamp: 70% of predictions still fail, the median still
    sits above the threshold, and the spread stays wide. It discriminates; it just stops demanding de-novo docking,
    which at 0.2% was not a workable pipeline.
  </div>
</div>

<h2>It worked</h2>
<p class="sub">Same target, same 15 hotspots, same filters, same Rosetta gates. One flag changed.</p>
<div class="panel">
  <div class="legend">
    <span><i class="sw" style="background:var(--md1)"></i>md1_block_af3 baseline</span>
    <span><i class="sw" style="background:var(--win)"></i>with predict_initial_guess</span>
  </div>
  {armA_svg()}
  <div class="tblwrap">
  <table>
  <thead><tr><th></th><th>baseline</th><th>with guess</th></tr></thead>
  <tbody>
  <tr><td>Runtime</td><td class="n">14 days</td><td class="n">17 hours</td></tr>
  <tr><td>Trajectory attempts</td><td class="n">617</td><td class="n">21</td></tr>
  <tr><td>Designs reaching Rosetta</td><td class="n">10</td><td class="n">51</td></tr>
  <tr><td><b>Accepted designs</b></td><td class="n zero">0</td><td class="n good">9</td></tr>
  <tr><td>Distinct source trajectories</td><td class="n">1 of 270</td><td class="n good">5 of 12</td></tr>
  </tbody>
  </table>
  </div>
  <div class="win">
    <span class="lab">Why these are real</span>
    The Rosetta filters were never touched — only the AF2 gate changed. The accepted designs clear the exact
    biophysical checks that killed all 10 baseline designs: <b>1.0–4.0 buried unsatisfied H-bonds</b> against a
    limit of 4 (baseline sat at 5–9 and every one died there), with roughly twice the interface H-bonds and better
    dG. And <b>Binder_RMSD 0.83–1.63 Å</b> — the binder predicted <em>alone</em> matches its pose in the complex,
    a check computed with no initial guess at all, so the flag cannot inflate it.
    Nine designs came from <b>five different trajectories</b>, so this is a broadly higher hit rate rather than
    one lucky backbone.
  </div>
  <p class="cap" style="margin-top:1rem">The wall also moved rather than vanishing: i_pAE now rejects 58% of
  sequences instead of 99.4%, and the leading rejection reason is Rosetta shape complementarity (32) followed by
  unsatisfied H-bonds (20). That is a normal design funnel.</p>
</div>

{hits_section()}
<h2>The queue so far</h2>
<p class="sub">All nine remaining targets are running with the flag, one at a time on a single GPU, each
stopping at 8 accepted designs or 50 relaxed trajectories. Two have results.</p>
<div class="tblwrap">
<table>
<thead><tr><th>Target</th><th>Relaxed</th><th>Scored</th><th>Accepted</th><th>Status</th></tr></thead>
<tbody>
<tr><td><code>clamp_md1md2_state3</code></td><td class="n">50</td><td class="n">197</td><td class="n good">5</td><td>complete — hit the 50-trajectory cap</td></tr>
<tr><td><code>clamp_md1md2_state1</code></td><td class="n">9</td><td class="n zero">0</td><td class="n zero">0</td><td>running</td></tr>
<tr><td colspan="5" style="color:var(--muted)">7 more queued: clamp_md1md2_state5, clamp_md2md3_state3/5, md3_block ×2, md2_block ×2</td></tr>
</tbody>
</table>
</div>

<h3 style="margin-top:2.2rem">State 3: the gate moved to interface chemistry</h3>
<div class="panel" style="margin-top:.8rem">
  {gates_svg()}
  <p class="cap">Share of the 197 scored designs failing each Rosetta gate. Energetics are not the problem —
  <b>dG and interface H-bond count fail 0%</b> (medians −58.8 and 10.5). The limiters are
  <b>buried unsatisfied H-bonds (68% fail, median 5.0 against a limit of 4)</b> and shape complementarity
  (68%, median 0.58 against 0.60). Both miss by a hair.</p>
  <div class="flag">
    <span class="lab">The same wall, one stage later</span>
    Unsatisfied polar burial is exactly what killed the original <code>md1_block_af3</code> campaign. Beating the
    i_pAE gate moved designs far enough down the funnel for it to reappear. A clamp has more polar surface to
    satisfy than a pocket plug, so a systematically higher unsat count may be intrinsic to the concept rather
    than a fixable defect — 5 designs still cleared it.
  </div>
</div>

<h3 style="margin-top:2.2rem">State 1: the flag cannot rescue a bad starting interface</h3>
<div class="panel" style="margin-top:.8rem">
  <div class="tblwrap">
  <table>
  <thead><tr><th>Median at trajectory stage</th><th>state 3</th><th>state 1</th></tr></thead>
  <tbody>
  <tr><td>i_pAE (threshold 0.35)</td><td class="n good">0.345</td><td class="n zero">0.540</td></tr>
  <tr><td>i_pTM</td><td class="n">0.705</td><td class="n">0.790</td></tr>
  <tr><td>pLDDT</td><td class="n">0.810</td><td class="n">0.800</td></tr>
  <tr><td>Scored designs</td><td class="n">197</td><td class="n zero">0</td></tr>
  </tbody>
  </table>
  </div>
  <p class="cap">State 1's trajectories fold just as well (pLDDT 0.80) but start with a markedly worse interface
  — i_pAE 0.540 against state 3's 0.345. <code>predict_initial_guess</code> helps AF2 re-find a pose the sequence
  already supports; it cannot manufacture an interface that was never good.</p>
  <div class="flag">
    <span class="lab">Consistent with the TICA clustering</span>
    State 3 was the <b>dominant basin (~80% of the population)</b>; state 1 a minor one (~6%). If clamp
    designability tracks state population, that is worth knowing before spending GPU time on the remaining
    states — but 9 trajectories is far too few to call, and the 50-trajectory cap bounds the cost of finding out.
  </div>
</div>

<h2>What is still unknown</h2>
<p class="sub">Whether the fix holds across the MD2 and MD3 pockets, and whether clamp quality really does track
TICA state population. Seven targets remain. The two MD2 blockers need unified memory — their 586-residue
construct exceeds a 46 GB card — so they run last and slowest.</p>

<footer>
The three baseline runs are stopped; the predict_initial_guess queue is in progress. Built from <code>Trajectory/</code> PDB coordinates,
<code>trajectory_stats.csv</code>, <code>mpnn_design_stats.csv</code> and <code>failure_csv.csv</code> under
<code>bindcraft_md/designs/</code>. Distances are minimum heavy-atom separation between the binder chain and
the union of that domain's hotspot residues. MPNN sequence counts are relaxed trajectories × <code>num_seqs=20</code>.
</footer>
</div>
"""
out=os.path.join(D,"parp14_geometry_report.html")
open(out,"w").write(html)
print("wrote",out,os.path.getsize(out)//1024,"KB")
