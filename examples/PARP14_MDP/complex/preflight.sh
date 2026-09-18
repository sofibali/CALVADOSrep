#!/usr/bin/env bash
# Run this ONCE on the GPU host before launching anything.
# Verifies: GPU visible, CALVADOS importable, a real 2000-step run completes.
set -u
PY="${PYTHON_EXE:-python}"
echo "== 1. GPU =="
nvidia-smi --query-gpu=index,name,memory.total --format=csv 2>/dev/null || { echo "  FAIL: no nvidia-smi"; exit 1; }
echo "== 2. OpenMM CUDA platform =="
"$PY" - <<'PYEOF' || exit 1
import openmm
try:
    p = openmm.Platform.getPlatformByName('CUDA')
    s = openmm.System(); s.addParticle(1.0)
    openmm.Context(s, openmm.VerletIntegrator(0.001), p)
    print("  OK: CUDA context created")
except Exception as e:
    print("  FAIL:", str(e).splitlines()[0]); raise SystemExit(1)
PYEOF
echo "== 3. CALVADOS =="
"$PY" -c "import calvados,os;print('  OK:',os.path.dirname(calvados.__file__))" || exit 1
echo "== 4. timestep sanity =="
"$PY" - <<'PYEOF'
import calvados.sim, inspect, re
src = inspect.getsource(calvados.sim)
m = re.search(r'LangevinMiddleIntegrator\([^)]*\)', src)
line = m.group(0) if m else '?'
if 'self.dt' in line:
    print("  patched build (dt read from config.yaml) - dt: 0.01 will be honoured")
elif '0.01' in line:
    print("  stock build (dt hardcoded 0.01 ps) - matches config dt: 0.01, OK")
else:
    print("  WARNING: unexpected integrator line ->", line)
PYEOF
echo "== 5. short real run (2000 steps, CPU) =="
d=$(mktemp -d); cp -rL binding/p9_dtx3l/rep-1/. "$d"/ 2>/dev/null || { echo "  FAIL: run from the dir containing binding/"; exit 1; }
( cd "$d" && "$PY" - <<'PYEOF'
import yaml
c=yaml.safe_load(open('config.yaml')); c.update(steps=2000,wfreq=500,platform='CPU',threads=4)
yaml.dump(c,open('config.yaml','w'),default_flow_style=False,sort_keys=False)
PYEOF
cd "$d" && "$PY" run.py > smoke.log 2>&1 && ls *.dcd >/dev/null 2>&1 \
  && echo "  OK: trajectory written" || { echo "  FAIL - see $d/smoke.log"; exit 1; } )
echo
echo "Preflight passed. Launch with:  bash binding/run_all.sh 3"
