#!/bin/bash
################################################################################
# E3SM Initial Condition Remapping Workflow
#
# This script remaps E3SM atmospheric initial conditions from one grid 
# resolution to another while properly adjusting surface pressure (PS) to 
# account for topography differences.
#
# INSTRUCTIONS FOR ADAPTING TO DIFFERENT GRID COMBINATIONS:
# 1. Set CASENAME to your simulation case name
# 2. Set TIMESTAMP to your desired initialization time (format: YYYY-MM-DD-NNNNN)
# 3. Set SOURCE_GRID to your source grid name (e.g., ne30np4, ne120np4)
# 4. Set TARGET_GRID to your target grid name (e.g., northamericax4v1np4, conusx4v1np4)
# 5. Set GRID_LABEL to match your target topography file naming convention
# 6. Update file paths in "FILE PATHS" section below:
#    - ESMF_MAP_FILE: Horizontal remapping weights from source to target
#    - SOURCE_TOPO_FILE: Topography file for source grid
#    - TARGET_TOPO_FILE: Topography file for target grid
#    - VERT_TEMPLATE: Vertical coordinate template file
# 7. Verify PHIS_VAR matches the variable name in your topography files
# 8. Ensure vertical levels (L80) match your model configuration
# 9. Set KEEP_INTERMEDIATE to "yes" if you want to keep intermediate files
#
################################################################################

# Load E3SM environment
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

################################################################################
# USER CONFIGURATION - MODIFY THESE VARIABLES
################################################################################

# Case identification
CASENAME="v3.LR.historical_0321"
TIMESTAMP="2015-01-01-00000"

# Grid specifications
SOURCE_GRID="ne30np4"
TARGET_GRID="northamericax4v1np4"
GRID_LABEL="northamericax4v1pg2"  # For topography file naming

# Vertical resolution
VERTICAL_LEVELS="L80"  # Number of vertical levels

# PHIS variable name in topography files
PHIS_VAR="PHIS_d"

# Keep intermediate files? (yes/no, default: no)
KEEP_INTERMEDIATE="no"

################################################################################
# FILE PATHS - UPDATE THESE FOR YOUR CONFIGURATION
################################################################################

# Input files
SOURCE_IC_FILE="${CASENAME}.eam.i.${TIMESTAMP}.nc"
ESMF_MAP_FILE="map_${SOURCE_GRID}_to_${TARGET_GRID}_esmfbilin.20250813.nc"
SOURCE_TOPO_FILE="/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_${SOURCE_GRID}pg2_x6t-SGH.c20210614.nc"
TARGET_TOPO_FILE="USGS-gtopo30_${TARGET_GRID}pg2_oroshp_x6t.c20250813.nc"
VERT_TEMPLATE="E3SM_vert_coor_${VERTICAL_LEVELS}.nc"

# Output file naming
MAPPED_IC_FILE="${CASENAME}_mapped_${TARGET_GRID}.eam.i.${TIMESTAMP}.nc"
VERTFILE="E3SM_vert_coor_${VERTICAL_LEVELS}_${GRID_LABEL}.${TIMESTAMP}.nc"
MAPPED_PHIS_FILE="${PHIS_VAR}.${SOURCE_GRID}_mapped_${TARGET_GRID}.nc"
FINAL_OUTPUT="${CASENAME}_mapped_${TARGET_GRID}-topoadj.eam.i.${TIMESTAMP}.nc"

################################################################################
# WORKFLOW EXECUTION - DO NOT MODIFY BELOW THIS LINE
################################################################################

echo "================================================================================"
echo "E3SM Initial Condition Remapping Workflow"
echo "================================================================================"
echo "Case name: ${CASENAME}"
echo "Timestamp: ${TIMESTAMP}"
echo "Source grid: ${SOURCE_GRID}"
echo "Target grid: ${TARGET_GRID}"
echo "Vertical levels: ${VERTICAL_LEVELS}"
echo "Keep intermediate files: ${KEEP_INTERMEDIATE}"
echo "================================================================================"

# Step 1: Horizontal remapping from source to target grid
echo ""
echo "Step 1: Horizontal remapping..."
echo "  Input: ${SOURCE_IC_FILE}"
echo "  Output: ${MAPPED_IC_FILE}"
echo "  Weights: ${ESMF_MAP_FILE}"

ncks -6 -O --map "${ESMF_MAP_FILE}" \
  "${SOURCE_IC_FILE}" \
  "${MAPPED_IC_FILE}"

if [ $? -ne 0 ]; then
  echo "ERROR: Horizontal remapping failed"
  exit 1
fi
echo "  ✓ Horizontal remapping completed"

# Step 2: Create vertical coordinate file with adjusted PS for topography
echo ""
echo "Step 2: Creating vertical file with topo-adjusted PS..."

# 2.1: Copy generic vertical coordinate file
echo "  Step 2.1: Copying vertical coordinate template..."
cp "${VERT_TEMPLATE}" "${VERTFILE}"

if [ $? -ne 0 ]; then
  echo "ERROR: Failed to copy vertical coordinate file"
  exit 1
fi
echo "  ✓ Vertical coordinate file copied"

# 2.2: Remap PHIS from source to target grid
echo ""
echo "  Step 2.2: Remapping PHIS..."
echo "    Source: ${SOURCE_TOPO_FILE}"
echo "    Output: ${MAPPED_PHIS_FILE}"

ncks -O -v "${PHIS_VAR}" --map "${ESMF_MAP_FILE}" \
  "${SOURCE_TOPO_FILE}" "${MAPPED_PHIS_FILE}"

if [ $? -ne 0 ]; then
  echo "ERROR: PHIS remapping failed"
  exit 1
fi
echo "  ✓ PHIS remapping completed"

# 2.3: Adjust PS for topography differences
echo ""
echo "  Step 2.3: Adjusting PS for topography..."
echo "    IC file: ${MAPPED_IC_FILE}"
echo "    Target topo: ${TARGET_TOPO_FILE}"
echo "    Mapped PHIS: ${MAPPED_PHIS_FILE}"

python adjust_ps_m2m.py \
  --vertfile "${VERTFILE}" \
  --rawicfile "${MAPPED_IC_FILE}" \
  --topofile "${TARGET_TOPO_FILE}" \
  --mapped-phis "${MAPPED_PHIS_FILE}" \
  --phis-var "${PHIS_VAR}"

if [ $? -ne 0 ]; then
  echo "ERROR: Python script failed"
  exit 1
fi
echo "  ✓ PS adjustment completed"

# Step 3: Vertical remapping using the adjusted PS
echo ""
echo "Step 3: Vertical remapping..."
echo "  Using vertical file: ${VERTFILE}"
echo "  Output: ${FINAL_OUTPUT}"

ncremap --vrt_fl="${VERTFILE}" "${MAPPED_IC_FILE}" "${FINAL_OUTPUT}"

if [ $? -ne 0 ]; then
  echo "ERROR: ncremap failed"
  exit 1
fi
echo "  ✓ Vertical remapping completed"

# Step 4: Ensure adjusted PS is used in the final file
echo ""
echo "Step 4: Copying adjusted PS to final file..."

ncks -A -v PS "${VERTFILE}" "${FINAL_OUTPUT}"

if [ $? -ne 0 ]; then
  echo "ERROR: Failed to copy PS to final file"
  exit 1
fi
echo "  ✓ PS copied to final file"

# Step 5: Clean up intermediate files if requested
echo ""
if [ "${KEEP_INTERMEDIATE}" = "no" ] || [ "${KEEP_INTERMEDIATE}" = "No" ] || [ "${KEEP_INTERMEDIATE}" = "NO" ]; then
  echo "Step 5: Cleaning up intermediate files..."
  
  if [ -f "${MAPPED_IC_FILE}" ]; then
    rm -f "${MAPPED_IC_FILE}"
    echo "  ✓ Removed ${MAPPED_IC_FILE}"
  fi
  
  if [ -f "${VERTFILE}" ]; then
    rm -f "${VERTFILE}"
    echo "  ✓ Removed ${VERTFILE}"
  fi
  
  if [ -f "${MAPPED_PHIS_FILE}" ]; then
    rm -f "${MAPPED_PHIS_FILE}"
    echo "  ✓ Removed ${MAPPED_PHIS_FILE}"
  fi
  
  echo "  ✓ Cleanup completed"
else
  echo "Step 5: Keeping intermediate files (as requested)"
  echo "  Intermediate files preserved:"
  echo "    - ${MAPPED_IC_FILE} (horizontally remapped)"
  echo "    - ${VERTFILE} (vertical coordinate with adjusted PS)"
  echo "    - ${MAPPED_PHIS_FILE} (remapped PHIS)"
fi

echo ""
echo "================================================================================"
echo "SUCCESS: All steps completed"
echo "================================================================================"
echo "Final output file: ${FINAL_OUTPUT}"

if [ "${KEEP_INTERMEDIATE}" != "no" ] && [ "${KEEP_INTERMEDIATE}" != "No" ] && [ "${KEEP_INTERMEDIATE}" != "NO" ]; then
  echo ""
  echo "Intermediate files:"
  echo "  - ${MAPPED_IC_FILE} (horizontally remapped)"
  echo "  - ${VERTFILE} (vertical coordinate with adjusted PS)"
  echo "  - ${MAPPED_PHIS_FILE} (remapped PHIS)"
fi

echo "================================================================================"
