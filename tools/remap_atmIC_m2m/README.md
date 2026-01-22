# E3SM Initial Condition Remapping Workflow

Remap E3SM atmospheric initial conditions from one grid to another with proper surface pressure adjustment for topography differences.

Uses may need to edit the following input variables in the driver `remap_atmIC_m2m.sh` (m2m stands for model-to-model).

```bash
CASENAME="v3.LR.historical_0091"
TIMESTAMP="2015-01-01-00000"
SOURCE_GRID="ne30np4"
TARGET_GRID="northamericax4v1np4"
GRID_LABEL="northamericax4v1pg2"
```

The final output file is `${CASENAME}_mapped_${TARGET_GRID}-topoadj.eam.i.${TIMESTAMP}.nc`.


## Quick Start

### 1. Get the Scripts

From E3SM repository:

```bash
# Clone or navigate to E3SM repository
cd /path/to/E3SM/tools/remap_atmIC_m2m/
cp -r . /path/to/your/working/directory/
cd /path/to/your/working/directory/
```

Or download directly from:
<https://github.com/E3SM-Project/E3SM/tree/master/tools/remap_atmIC_m2m/>

Example files available at:
<https://portal.nersc.gov/cfs/e3sm/tang30/remap_atmIC_m2m/>

### 2. Prepare Input Files

You need:

- Source initial condition file: `${CASENAME}.eam.i.${TIMESTAMP}.nc`
- Remapping weights on the "np4" grids (TempestRemap is recommended.): `map_source_to_target.nc`
- Source grid topography file
- Target grid topography file  
- Vertical coordinate template: `E3SM_vert_coor_L80.nc`

### 3. Edit Configuration

Edit `remap_atmIC_m2m.sh`:

```bash
CASENAME="v3.LR.historical_0091"
TIMESTAMP="2015-01-01-00000"
SOURCE_GRID="ne30np4"
TARGET_GRID="northamericax4v1np4"
GRID_LABEL="northamericax4v1pg2"
```

### 4. Run

```bash
./remap_atmIC_m2m.sh
```

Output: `${CASENAME}_mapped_${TARGET_GRID}-topoadj.eam.i.${TIMESTAMP}.nc`

## Files

- `remap_atmIC_m2m.sh` - Main workflow script
- `adjust_ps_m2m.py` - Python script for PS adjustment
- `README.md` - This file

## What It Does

1. **Horizontal remapping** - Maps IC from source to target grid
2. **PS adjustment** - Adjusts surface pressure for topography differences
3. **Vertical remapping** - Interpolates to new vertical coordinate
4. **Cleanup** - Removes intermediate files (set `KEEP_INTERMEDIATE="yes"` to keep)

## Configuration Options

In `remap_atmIC_m2m.sh`:

| Variable            | Description          | Example                    |
|---------------------|----------------------|----------------------------|
| `CASENAME`          | Case name            | `v3.LR.historical_0091`    |
| `TIMESTAMP`         | Time stamp           | `2015-01-01-00000`         |
| `SOURCE_GRID`       | Source grid          | `ne30np4`                  |
| `TARGET_GRID`       | Target grid          | `northamericax4v1np4`      |
| `GRID_LABEL`        | For file naming      | `northamericax4v1pg2`      |
| `VERTICAL_LEVELS`   | Vertical levels      | `L80`                      |
| `PHIS_VAR`          | PHIS variable name   | `PHIS_d`                   |
| `KEEP_INTERMEDIATE` | Keep temp files      | `no` (default)             |

## Requirements

- E3SM unified environment
- Python 3 with numpy, netCDF4
- NCO tools (ncks, ncremap)

Load E3SM unified environment:

```bash
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```
