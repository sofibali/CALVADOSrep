# FOXP Model Simulation - Change Log

## Date: 2026-02-02

### Bug Fixes to CALVADOS Core (`/home/sbali/CALVADOS/calvados/`)

#### 1. `components.py` - Line 43
**Issue:** `Component.calc_properties()` missing `comp_setup` parameter
**Fix:** Added `comp_setup: str = 'linear'` parameter to base `Component.calc_properties()` method
```python
# Before
def calc_properties(self, pH: float = 7.0, verbose: bool = False):

# After
def calc_properties(self, pH: float = 7.0, verbose: bool = False, comp_setup: str = 'linear'):
```

#### 2. `components.py` - Line ~108 (before Protein class)
**Issue:** Base `Component` class missing `init_restraint_force()` method, causing crash when generic components have restraints enabled
**Fix:** Added stub method to base `Component` class:
```python
def init_restraint_force(self, eps_lj=None, cutoff_lj=None, eps_yu=None, k_yu=None):
    """ Initialize restraint force (stub for base Component). """
    pass
```

---

### New Simulation Configurations Created

#### 1. Proteins Only Setup
**Location:** `/home/sbali/CALVADOS/examples/foxP_model/proteins_only/`

**Files:**
- `config.yaml` - Simulation parameters (5M steps, 50nm box, grid topology)
- `components.yaml` - FOXP4, FOX, chain_C with k_harmonic=2000.0
- `run.py` - Execution script

**Purpose:** Quick simulation of 3 proteins with highly restrained domains

**Run command:**
```bash
cd /home/sbali/CALVADOS/examples/foxP_model/proteins_only
conda activate CALVADOS
python run.py
```

#### 2. Proteins + DNA (as RNA) Setup
**Location:** `/home/sbali/CALVADOS/examples/foxP_model/proteins_with_DNA/`

**Files:**
- `config.yaml` - Simulation parameters
- `components.yaml` - 3 proteins + 2 DNA chains (as coarse-grained RNA)
- `sequences.fasta` - Protein sequences + RNA sequences (21 'r' characters each for chain_E/F)
- `run.py` - Execution script

**Purpose:** Simulate protein-DNA complex with DNA modeled as coarse-grained RNA

**Run command:**
```bash
cd /home/sbali/CALVADOS/examples/foxP_model/proteins_with_DNA
conda activate CALVADOS
python run.py
```

---

### Key Configuration Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `k_harmonic` | 2000.0 | Very strong restraints for rigid domains |
| `cutoff_restr` | 0.9 nm | Distance cutoff for applying restraints |
| `topol` | grid | Grid placement for multiple molecules |
| `steps` | 5,000,000 | Total simulation steps |
| `wfreq` | 5,000 | Write frequency (produces ~1000 frames) |
| `box` | 50x50x50 nm | Simulation box size |

---

### Notes

1. **CALVADOS does not support `molecule_type: ion`** - Magnesium ions (chain_D) cannot be simulated directly
2. **DNA is modeled as RNA** - CALVADOS uses coarse-grained 2-bead RNA model (backbone 'p' + base 'r')
3. **Residues file for mixed systems:** Use `residues_C2RNA.csv` which contains both amino acids and RNA beads
4. **Visualization:** Use `--skip 0` with analyze_trajectory.py if you have < 100 frames

---

### Domains Configuration
Defined in `/home/sbali/CALVADOS/examples/foxP_model/input/domains.yaml`:
- FOXP4: [118-199], [299-375], [451-551]
- FOX: [118-299], [319-375], [450-539]
- chain_C: [232-321]
- chain_E: [1-21]
- chain_F: [1-21]

Residues within these domains get harmonic restraints; disordered regions outside move freely.
