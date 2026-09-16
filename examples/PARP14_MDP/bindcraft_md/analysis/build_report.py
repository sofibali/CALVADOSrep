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

html=f"""<title>Three Ways To Miss A Pocket</title>
<style>
:root {{
  --paper:#f5f7fa; --panel:#ffffff; --ink:#111418; --ink2:#4a5463; --muted:#8a93a3;
  --line:rgba(17,20,24,.12); --grid:rgba(17,20,24,.09);
  --md1:#2a78d6; --md2:#eb6834; --sim:#1baf7a; --loss:#96a0b0; --quad:rgba(42,120,214,.09);
  --flag:#b8442a; --flagbg:#fdeee8;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --paper:#0f1216; --panel:#171b21; --ink:#eef2f7; --ink2:#aab4c2; --muted:#7f8a99;
  --line:rgba(255,255,255,.14); --grid:rgba(255,255,255,.10);
  --md1:#5b9bec; --md2:#f5854f; --sim:#2cc98d; --loss:#5d6773; --quad:rgba(91,155,236,.13);
  --flag:#ff9a76; --flagbg:#2e1a13; }} }}
:root[data-theme="dark"] {{
  --paper:#0f1216; --panel:#171b21; --ink:#eef2f7; --ink2:#aab4c2; --muted:#7f8a99;
  --line:rgba(255,255,255,.14); --grid:rgba(255,255,255,.10);
  --md1:#5b9bec; --md2:#f5854f; --sim:#2cc98d; --loss:#5d6773; --quad:rgba(91,155,236,.13);
  --flag:#ff9a76; --flagbg:#2e1a13; }}
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
<h1>Three ways to miss a pocket</h1>
<p class="dek">Three binder campaigns, ~33 GPU-days, <b>zero accepted designs</b>. They do not fail for the
same reason, and none fails where I first assumed. This traces every attempt from hallucination to final
filter using the trajectory coordinates themselves — and finds one filter that kills almost everything.</p>

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

<footer>
All three runs are now stopped. Built from <code>Trajectory/</code> PDB coordinates,
<code>trajectory_stats.csv</code>, <code>mpnn_design_stats.csv</code> and <code>failure_csv.csv</code> under
<code>bindcraft_md/designs/</code>. Distances are minimum heavy-atom separation between the binder chain and
the union of that domain's hotspot residues. MPNN sequence counts are relaxed trajectories × <code>num_seqs=20</code>.
</footer>
</div>
"""
out=os.path.join(D,"parp14_geometry_report.html")
open(out,"w").write(html)
print("wrote",out,os.path.getsize(out)//1024,"KB")
