# MOSART Bifurcation and IBT System - Complete Reference Guide

This comprehensive reference guide contains all implementation details, examples, and technical specifications for the MOSART bifurcation and Inter-Basin Transfer (IBT) system.

## Overview

MOSART (Model for Scale Adaptive River Transport) is the river routing model component of E3SM. The bifurcation system enables river flow splitting at delta locations and managed water diversions, working exclusively with the KW (kinematic wave) routing scheme.

## Bifurcation and IBT Flag Combinations

The bifurcation system operates in different modes depending on the combination of `bifurcflag` and `ibtflag` settings:

| **`bifurcflag`** | **`ibtflag`** | **Mode** | **Description** | **Required NetCDF Variables** | **Behavior** |
|------------------|---------------|----------|-----------------|-------------------------------|--------------|
| `.false.` | `.false.` | **Standard MOSART** | Normal routing without bifurcation | None | Original MOSART routing behavior |
| `.false.` | `.true.` | ❌ **ERROR** | Invalid configuration | N/A | Model aborts: IBT requires bifurcation |
| `.true.` | `.false.` | **Pure Bifurcation** | Natural flow splitting only | `dnID(cells,downstream)` | Equal splits or custom ratios |
| `.true.` | `.true.` | **Mixed Mode** | Natural splitting + managed diversions | `dnID` + some `ibt_demand` | IBT volumes + equal splits |

### Mode Details

#### **Standard MOSART** (`bifurcflag=false, ibtflag=false`)
- **Purpose:** Original MOSART functionality
- **Routing:** Single downstream connection per cell
- **NetCDF:** Standard `dnID(cells)` or `dnID(lon,lat)` format
- **Physics:** KW or DW routing methods supported

#### **Pure Bifurcation Mode** (`bifurcflag=true, ibtflag=false`)  
- **Purpose:** Natural river bifurcation (deltas, distributaries)
- **Routing:** Multiple downstream connections with fixed ratios
- **NetCDF:** `dnID(cells,downstream)` required, `bifurc_ratio(cells,downstream)` optional
- **Physics:** KW routing only (DW incompatible)
- **Fallback:** Equal splits if no ratio data provided

**Example Namelist:**
```fortran
&mosart_inparm
  bifurcflag = .true.
  ibtflag = .false.
  RoutingMethod = 1  ! KW routing required
/
```

#### **Mixed Mode** (`bifurcflag=true, ibtflag=true`)
- **Purpose:** Natural bifurcation + managed water diversions  
- **Routing:** Heterogeneous - some cells use volumes, others use ratios
- **NetCDF:** `dnID` required, `ibt_demand` required for some cells, `bifurc_ratio` optional
- **Physics:** KW routing only, dynamic ratio calculations
- **Advanced:** Environmental flow protection, conflict detection

**Example Namelist:**
```fortran
&mosart_inparm
  bifurcflag = .true.
  ibtflag = .true.
  RoutingMethod = 1  ! KW routing required
/
```

### Configuration Requirements by Mode

#### **Pure Bifurcation Mode Requirements:**
```text
✅ Minimum: dnID(cells,downstream) connectivity data
✅ Optional: bifurc_ratio(cells,downstream) for custom ratios  
✅ Fallback: Equal splits (50-50, 33-33-33, 25-25-25-25)
✅ Compatible: Any bifurcation point can use equal splits
```

#### **Mixed Mode Requirements:**
```text
✅ Required: dnID(cells,downstream) connectivity data
✅ Required: ibt_demand(cells,downstream) for at least some bifurcation points
✅ Optional: bifurc_ratio(cells,downstream) for other bifurcation points
❌ Forbidden: Same cell cannot have both ibt_demand > 0 AND bifurc_ratio > 0
✅ Mixed Usage: Some cells use IBT volumes, others use equal splits
```

## MCT Communication Options

MOSART uses MCT (Model Coupling Toolkit) for parallel communication and supports two matrix building modes:

- **`smat_option = 'opt'`** - Distributed matrix building (each processor builds local entries)
- **`smat_option = 'Xonly'`** - Centralized matrix building (master processor builds global matrix)
- **`smat_option = 'Yonly'`** - Alternative centralized mode (same as 'Xonly')

**Bifurcation Compatibility:**
- ✅ **'opt' mode:** Full bifurcation support with dynamic matrix weights
- ✅ **'Xonly'/'Yonly' modes:** Full bifurcation support with 9-field data structure
- ✅ **Both modes:** Maintain bit-for-bit compatibility when `bifurcflag = .false.`

## Implementation Details

### Data Structures

#### Input Data Structure - NetCDF Parameter File Support (Both Mesh Types)

| **Mesh Type** | **Bifurcation Format** | **Legacy Format** | **Example Dimensions** |
|---------------|------------------------|-------------------|------------------------|
| **Structured** (`isgrid2d = true`) | `dnID(lon,lat,downstream)` | `dnID(lon,lat)` | `dnID(144,96,4)` |
| **Unstructured** (`isgrid2d = false`) | `dnID(cells,downstream)` | `dnID(cells)` | `dnID(50000,4)` |

**Smart NetCDF Reading:**
- Tries multi-dimensional format first, falls back to legacy format automatically
- Handles both structured `(lon,lat,downstream)` and unstructured `(cells,downstream)` formats
- `downstream` dimension can be any size ≤ `max_downstream` (4)

**Missing values:** Use -999 for unused downstream connections
**Automatic detection:** Bifurcation points identified by non-missing values in `downstream > 1`
**Equal split ratios:** Automatically calculated (50-50, 33-33-33, 25-25-25-25)

**Usage examples:**

```text
// Structured mesh bifurcation:
dnID(lon=123, lat=45, downstream=1) = 5001  // Primary
dnID(lon=123, lat=45, downstream=2) = 5002  // Secondary
dnID(lon=123, lat=45, downstream=3) = -999  // Unused
dnID(lon=123, lat=45, downstream=4) = -999  // Unused

// Unstructured mesh bifurcation:
dnID(cell=1500, downstream=1) = 2001  // Primary
dnID(cell=1500, downstream=2) = 2002  // Secondary
dnID(cell=1500, downstream=3) = -999  // Unused
dnID(cell=1500, downstream=4) = -999  // Unused
```

#### Memory Allocation (`RunoffMod.F90`)

```fortran
rtmCTL%dsig_all(begr:endr,max_downstream),     ! All downstream IDs
rtmCTL%rdsig_all(begr:endr,max_downstream),    ! Real version of downstream IDs  
rtmCTL%iDown_all(begr:endr,max_downstream),    ! Local downstream indices
rtmCTL%num_downstream(begr:endr),              ! Number of downstream connections per cell
rtmCTL%bifurc_ratio(begr:endr,max_downstream), ! Split ratios for each connection
rtmCTL%is_bifurc(begr:endr)                    ! Bifurcation point flag
```

**Initialization Values:**
- `dsig_all`, `rdsig_all`: Initialize to -999 (missing value)
- `iDown_all`: Initialize to 0 
- `num_downstream`: Initialize to 1 (backward compatibility)
- `bifurc_ratio`: Initialize to 0.0
- `is_bifurc`: Initialize to `.false.`

### Bifurcation Flow Splitting Logic

**Implementation Details:**
- **Location:** `MOSART_physics_mod.F90` lines 917-920, 1148-1151, 1261-1264
- **Subroutine Added:** `ApplyBifurcationSplitting()` (lines 939-974)

**Trigger:** Activates when `bifurcflag = .true.` AND `rtmCTL%is_bifurc(iunit) = .true.`
**Water Balance:** Routes split flows to ALL downstream connections (including primary)
**Flow Calculation:** `split_flow = total_erout * bifurc_ratio(iunit,j)` for each downstream `j`
**Direct Routing:** Applies split flows directly to `TRunoff%erin(downstream_unit,nt)`
**Prevent Double-Routing:** Sets `total_erout = 0` to disable normal routing mechanism

**Supported Routing Methods:**
- ✅ `Routing_KW` - Kinematic Wave (primary target and only supported method)
- ❌ `Routing_DW` - Diffusive Wave (incompatible - aborts with error)
- ❌ `Routing_DW_ocn_rof_two_way` - Ocean-coupled Diffusive Wave (incompatible - aborts with error)

### MCT Matrix Communication Implementation

**MCT Mode Support:**

1. **'opt' Mode (Distributed Matrix Building):**
   - Each processor builds matrix entries for its local cells
   - Uses existing 2-field aVect structure: `f1:f2` (source, downstream)
   - Matrix weights set dynamically based on bifurcation ratios

2. **'Xonly' Mode (Centralized Matrix Building):**
   - Master processor builds global matrix from gathered data
   - **Challenge:** Variable downstream connections (1-4 per cell) caused MCT gather size mismatch
   - **Solution:** Fixed-size 9-field data structure for consistent MPI communication
   - **Field Layout:**
     - `f1` = source gindex
     - `f2-f5` = downstream IDs (up to 4 connections)
     - `f6-f9` = corresponding weights
   - **B4B Compatibility:** Conditional logic preserves original 2-field structure when `bifurcflag = .false.`

## Inter-Basin Transfer (IBT) Extension

### Overview
The IBT extension provides flexible support for both natural bifurcation and managed water diversions, including inter-basin transfers with environmental flow protection.

### IBT Operation Modes

#### Mode 1: Pure Ratio-Based (`ibtflag = false`)
**Use Case:** Natural bifurcation only, fixed percentage diversions
- **NetCDF Input:** `bifurc_ratio(lon,lat,downstream)` - values between 0.0-1.0
- **Validation:** Ratios must sum to 1.0 (within 1e-6 tolerance), auto-normalized if needed
- **Examples:** 
  - Delta bifurcation: [0.5, 0.5, 0.0, 0.0]
  - Irrigation canal: [0.9, 0.1, 0.0, 0.0] (90% stays, 10% diverted)

#### Mode 2: Mixed Mode (`ibtflag = true`)
**Use Case:** Heterogeneous basins with both natural deltas AND managed diversions
- **NetCDF Input:** BOTH variables supported in same file:
  - `bifurc_ratio(lon,lat,downstream)` - fixed ratios for delta points
  - `ibt_demand(lon,lat,downstream)` - absolute demands for IBT points (m³/s)
- **Per-Cell Logic:** 
  - Cells with `ibt_demand > 0`: Dynamic ratios calculated from demands + environmental protection
  - Cells with `ibt_demand = 0`: Use fixed ratios from `bifurc_ratio` (or equal splits)
- **Unified Processing:** All cells use the same MCT bifurcation infrastructure

### Implementation Details

#### Data Structures (`RunoffMod.F90`)
```fortran
real(r8), pointer :: ibt_demand(:,:)     ! IBT demands in m³/s (cell, downstream)
real(r8), pointer :: bifurc_ratio(:,:)   ! Split ratios (cell, downstream)
logical, pointer :: is_bifurc(:)         ! Bifurcation point flag
```

#### NetCDF File Support
**Both Grid Types Supported:**
- **Structured:** `variable(lon, lat, downstream)`
- **Unstructured:** `variable(cells, downstream)`  

**NetCDF Variables by Mode:**
- **Pure Ratio Mode (`ibtflag=false`):** 
  - Required: `dnID` - downstream connectivity (3D format)
  - Required: `bifurc_ratio` - split ratios (0.0-1.0)
- **Mixed Mode (`ibtflag=true`):**
  - Required: `dnID` - downstream connectivity (3D format)  
  - Optional: `bifurc_ratio` - fixed ratios for delta points (0.0-1.0)
  - Optional: `ibt_demand` - absolute demands for IBT points (m³/s)
  - **Note:** Can have both variables in same file for heterogeneous basins

#### Environmental Flow Protection
- **Minimum Flow:** 10% of current `erout` (hard-coded)
- **Priority Order:** Environmental flow > Legal diversions > Remaining flow
- **Warning System:** Logs when environmental minimum is triggered
- **Water Balance:** Perfect conservation maintained

#### Unified Dynamic Ratio Calculation
**Location:** `MOSART_physics_mod.F90` - `UpdateIBTRatios()` subroutine (lines 1224-1310)

**Key Insight:** Convert IBT demands to dynamic ratios, then use existing MCT bifurcation infrastructure

```fortran
! Mixed Mode Logic:
do iunit = begr, endr
   if (is_bifurcation_cell) then
      ! Check if this cell has IBT demands
      has_ibt_demands = any(ibt_demand(cell,:) > 0)
      
      if (has_ibt_demands) then
         ! Dynamic ratio calculation from demands
         total_flow = abs(erout(cell))
         env_minimum = 0.1 * total_flow  ! 10% environmental protection
         available_flow = max(0.0, total_flow - env_minimum)
         
         ! Priority-based allocation
         do k = 2, num_downstream
            demand = ibt_demand(cell,k)
            actual_diversion(k) = min(demand, available_flow)
            available_flow = available_flow - actual_diversion(k)
         end do
         
         ! Convert to ratios for MCT system
         bifurc_ratio(cell,k) = actual_diversion(k) / total_flow
         bifurc_ratio(cell,1) = 1.0 - sum(bifurc_ratio(cell,2:max_downstream))
      else
         ! Keep existing fixed ratios unchanged (from bifurc_ratio or equal splits)
         ! No modifications needed - ratios set during initialization
      endif
   endif
end do

! MCT matrices use updated ratios for ALL cells (unified approach)
```

#### Dynamic MCT Matrix Weight Updates
**Location:** `MOSART_physics_mod.F90` - `UpdateMCTMatrixWeights()` subroutine (lines 1321-1360)

**Technical Achievement:** Successfully identified and implemented correct MCT API for runtime matrix weight updates:

```fortran
! Correct MCT Data Structure Access:
! SMatP_upstrm (type: SparseMatrixPlus)
! └── %Matrix (type: SparseMatrix) 
!     └── %data (type: AttrVect)
!         └── %rAttr(iwgt, matrix_index) = new_weight

subroutine UpdateMCTMatrixWeights()
  ! Get weight attribute index using correct API
  iwgt = mct_sMat_indexRA(SMatP_upstrm%Matrix, 'weight')
  
  ! Selective updates - only bifurcation points
  do iunit = begr, endr
     if (is_bifurc(nr)) then
        do k = 1, num_downstream(nr)
           matrix_index = matrix_idx(nr,k)  ! Stored during initialization
           if (matrix_index > 0) then
              ! Update specific matrix weight with current dynamic ratio
              SMatP_upstrm%Matrix%data%rAttr(iwgt, matrix_index) = bifurc_ratio(nr,k)
              num_updates = num_updates + 1
           endif
        end do
     endif
  end do
end subroutine
```

### Use Cases

#### 1. **Pure Natural Delta Bifurcation**
```text
ibtflag = false
bifurc_ratio(delta_cell, 1) = 0.5  ! 50% to main channel
bifurc_ratio(delta_cell, 2) = 0.5  ! 50% to secondary channel
```

#### 2. **Pure IBT Water Rights Diversion**
```text
ibtflag = true
ibt_demand(diversion_point, 1) = 0.0   ! Primary gets remainder
ibt_demand(diversion_point, 2) = 25.0  ! 25 m³/s legal right
```

#### 3. **Mixed Mode: Delta + Municipal Withdrawal**
```text
ibtflag = true
# NetCDF file contains BOTH variables:

# Natural delta (fixed ratios):
bifurc_ratio(delta_cell, 1) = 0.6    ! 60% to main channel
bifurc_ratio(delta_cell, 2) = 0.4    ! 40% to secondary channel  
ibt_demand(delta_cell, :) = 0.0      ! No IBT demands

# Municipal withdrawal (dynamic ratios):
bifurc_ratio(city_cell, :) = 0.0     ! No fixed ratios needed
ibt_demand(city_cell, 1) = 0.0       ! Primary gets remainder
ibt_demand(city_cell, 2) = 20.0      ! 20 m³/s municipal demand
```

#### 4. **Complex Mixed Basin**
```text
ibtflag = true
# Heterogeneous basin with multiple process types:

# Mississippi Delta (natural 50-50):
bifurc_ratio(delta_lon, delta_lat, 1) = 0.5
bifurc_ratio(delta_lon, delta_lat, 2) = 0.5
ibt_demand(delta_lon, delta_lat, :) = 0.0

# New Orleans Municipal (20 m³/s):
ibt_demand(city_lon, city_lat, 1) = 0.0       ! Primary gets remainder
ibt_demand(city_lon, city_lat, 2) = 20.0      ! Municipal withdrawal

# Industrial Cooling (30 m³/s):
ibt_demand(plant_lon, plant_lat, 1) = 0.0     ! Primary gets remainder  
ibt_demand(plant_lon, plant_lat, 2) = 30.0    ! Cooling water

# Inter-basin Transfer (100 m³/s):
ibt_demand(transfer_lon, transfer_lat, 1) = 0.0    ! Primary gets remainder
ibt_demand(transfer_lon, transfer_lat, 2) = 100.0  ! To distant basin
```

#### 5. **Real-World Application: California Central Valley**
```text
ibtflag = true
# Sacramento-San Joaquin Delta + Diversions:

# Natural delta splits (Sacramento River):
bifurc_ratio(sac_delta, 1) = 0.7     ! 70% to San Francisco Bay
bifurc_ratio(sac_delta, 2) = 0.3     ! 30% to southern distributaries

# State Water Project (variable demand):
ibt_demand(swp_intake, 1) = 0.0      ! Primary gets remainder
ibt_demand(swp_intake, 2) = 150.0    ! 150 m³/s to Southern California

# Central Valley Project:
ibt_demand(cvp_intake, 1) = 0.0      ! Primary gets remainder  
ibt_demand(cvp_intake, 2) = 200.0    ! 200 m³/s agricultural demand

# Local municipal (Stockton):
ibt_demand(stockton, 1) = 0.0        ! Primary gets remainder
ibt_demand(stockton, 2) = 15.0       ! 15 m³/s municipal supply
```

## NetCDF Parameter File Structure

**Unified Variable Support:**

| **Variable** | **Units** | **Range** | **Usage** | **Required** |
|--------------|-----------|-----------|-----------|--------------|
| `dnID` | dimensionless | integer IDs | Downstream connectivity | Always |
| `bifurc_ratio` | dimensionless | 0.0 - 1.0 | Fixed split ratios for delta points | Ratio mode only |
| `ibt_demand` | m³/s | ≥ 0.0 | Absolute demands for IBT points | Mixed mode optional |

**Example NetCDF File Content:**
```text
// Mixed Mode (ibtflag=.true.) - BOTH variables in same file:
dimensions:
  lon = 144 ; lat = 96 ; downstream = 4 ;
variables:
  int dnID(lon,lat,downstream) ;             // Required: connectivity
  float bifurc_ratio(lon,lat,downstream) ;   // Optional: fixed ratios for deltas  
  float ibt_demand(lon,lat,downstream) ;     // Optional: demands for IBT points
  
// Pure Ratio Mode (ibtflag=.false.):
dimensions:
  lon = 144 ; lat = 96 ; downstream = 4 ;
variables:
  int dnID(lon,lat,downstream) ;             // Required: connectivity
  float bifurc_ratio(lon,lat,downstream) ;   // Required: ratios 0.0-1.0
```

**Mixed Mode Example Data:**
```text
// Sacramento Delta (natural 70-30 split):
dnID(sac_delta, :) = [bay_cell, south_cell, -999, -999]
bifurc_ratio(sac_delta, :) = [0.7, 0.3, 0.0, 0.0]    // Fixed delta ratios
ibt_demand(sac_delta, :) = [0.0, 0.0, 0.0, 0.0]      // No IBT demands

// Municipal intake (20 m³/s withdrawal):  
dnID(city_intake, :) = [river_cell, treatment_plant, -999, -999]
bifurc_ratio(city_intake, :) = [0.0, 0.0, 0.0, 0.0]  // No fixed ratios
ibt_demand(city_intake, :) = [0.0, 20.0, 0.0, 0.0]   // 20 m³/s municipal demand
```

## Validation and Testing

### Configuration Requirements and Validation

**CRITICAL: IBT Mode Prerequisites**
IBT mode has strict configuration requirements that are automatically validated:

1. **`ibtflag` requires `bifurcflag = .true.`**
   ```
   Error: ibtflag=true requires bifurcflag=true
   IBT (Inter-Basin Transfer) mode is an extension of the bifurcation system
   ```

2. **IBT only supports KW routing (`RoutingMethod = 1`)**
   ```
   Error: Inter-basin transfer (IBT) is currently only supported with KW routing method
   Current routing method: 2 (1=KW, 2=DW)
   ```

3. **NetCDF Variables by Mode**
   - **Pure Ratio Mode (`ibtflag=false`):** Must contain `bifurc_ratio(lon,lat,downstream)`
   - **Mixed Mode (`ibtflag=true`):** Can contain both `ibt_demand` AND `bifurc_ratio` variables
   - **All Modes:** Must contain `dnID(lon,lat,downstream)` for connectivity

### Validation Rules by Mode

| **Mode** | **Requirements** | **Behavior** | **Error Conditions** |
|----------|------------------|--------------|---------------------|
| **Pure Bifurcation** (`ibtflag=false`) | `dnID` connectivity only | ✅ Equal splits if no ratio data<br>✅ Partial ratio support | ❌ Missing `dnID` variable<br>❌ Ratios sum > 1.0 |
| **Mixed Mode** (`ibtflag=true`) | `dnID` + some `ibt_demand` data | ✅ Mixed IBT + ratio points<br>✅ Partial ratio support | ❌ No IBT data found<br>❌ Same point has both data types<br>❌ Ratios sum > 1.0 |

### Enhanced Validation Features

#### **Partial Ratio Support** (New Feature)
- ✅ **Input Flexibility:** Users can specify only known ratios, MOSART calculates remainder
- ✅ **Automatic Adjustment:** Primary downstream gets remaining ratio when sum < 1.0
- ✅ **Error Prevention:** Aborts simulation when ratios sum > 1.0 (impossible allocation)
- ✅ **Clear Messaging:** Detailed warnings explain what adjustments were made
- ✅ **Backward Compatible:** Existing complete ratio files work unchanged

#### **Robust Error Detection**
- ✅ **Connectivity Validation:** Checks all downstream IDs for validity
- ✅ **Circular Reference Detection:** Prevents self-referencing connections
- ✅ **Water Balance Enforcement:** Ensures all ratios sum to 1.0 within tolerance
- ✅ **IBT Conflict Detection:** Prevents mixed data types for same cell

**For detailed scenarios and examples, see [Partial Ratio Handling and Smart Validation](#partial-ratio-handling-and-smart-validation) section.**

### Error Examples and Solutions

**1. IBT Mode Without Volume Data:**
```text
ERROR: IBT mode enabled (ibtflag=true) but no IBT volume data found
       When ibtflag=true, NetCDF file must contain ibt_demand variable
       with volume data (m³/s) for at least some bifurcation points
       Either:
         1. Add ibt_demand(lon,lat,downstream) variable to NetCDF file, OR
         2. Set ibtflag=false for pure bifurcation mode (allows equal splits)
```

**2. Conflicting Data Types:**
```text
ERROR: Bifurcation point conflict at grid cell 12345
       Cell has both IBT volume data AND ratio data
       Each bifurcation point must use either:
         - IBT volumes: ibt_demand(cell,downstream) > 0, OR
         - Fixed ratios: bifurc_ratio(cell,downstream) > 0
       But not both for the same cell
```

**3. Successful Validation Summary:**
```text
ValidateBifurcationIBTData: Bifurcation/IBT validation summary:
  Total bifurcation points found: 25
  IBT points with volume data: 8
  Regular bifurcation points (equal splits): 17
All bifurcation/IBT data passed validation
```

## Partial Ratio Handling and Smart Validation

### Overview
MOSART now provides flexible handling of partial ratio inputs, allowing users to provide incomplete ratio data while maintaining water conservation. This feature eliminates the need for users to manually calculate exact ratios for all downstream connections.

### Enhanced Validation Logic

**Implementation Location:** `RunoffMod.F90` - `ValidateAndAdjustRatios()` subroutine (lines 865-968)

The new validation system intelligently handles various partial ratio scenarios:

### Scenario 1: N Downstream Connections with M Ratios (N > M)

**Condition:** More downstream connections than defined ratios  
**Example:** 3 downstream connections but only 2 ratios provided

#### Sub-case 1A: Ratios Sum ≤ 1.0 (Valid Partial Input)
**Behavior:** ✅ **WARNING** - Assigns remaining ratio to primary downstream  
**Water Balance:** Perfect conservation maintained

**Example:**
```text
Input NetCDF Data:
dnID(cell_123, :) = [5001, 5002, 5003, -999]        # 3 downstream connections
bifurc_ratio(cell_123, :) = [0.0, 0.3, 0.4, 0.0]   # Only connections 2&3 defined (sum=0.7)

MOSART Processing:
ValidateAndAdjustRatios WARNING: 3 downstream connections with 2 ratios found for point 12345
- assigning remaining ratio 0.300000 to primary downstream

Final Ratios:
bifurc_ratio(cell_123, :) = [0.3, 0.3, 0.4, 0.0]   # Primary gets remaining 30%
```

**Message Format:**
```text
ValidateAndAdjustRatios WARNING: N downstream connections with M ratios found for point XXXXX
- assigning remaining ratio R.RRRRRR to primary downstream
```

#### Sub-case 1B: Ratios Sum > 1.0 (Invalid Input)
**Behavior:** ❌ **ERROR** - Simulation aborts  
**Reason:** Cannot assign negative ratio to primary downstream

**Example:**
```text
Input NetCDF Data:
dnID(cell_456, :) = [5001, 5002, 5003, -999]        # 3 downstream connections  
bifurc_ratio(cell_456, :) = [0.0, 0.6, 0.7, 0.0]   # Connections 2&3 sum to 1.3

MOSART Processing:
ValidateAndAdjustRatios ERROR: 3 downstream connections with 2 ratios found for point 45678
but ratios sum to 1.300000 > 1.0

Result: Simulation aborts
```

**Message Format:**
```text
ValidateAndAdjustRatios ERROR: N downstream connections with M ratios found for point XXXXX
but ratios sum to S.SSSSSS > 1.0
```

### Scenario 2: Primary Downstream Pre-defined

**Condition:** Primary downstream ratio already specified in input  
**Behavior:** ✅ **INFO** - Reports configuration, no adjustments made

**Example:**
```text
Input NetCDF Data:
dnID(cell_789, :) = [5001, 5002, 5003, -999]        # 3 downstream connections
bifurc_ratio(cell_789, :) = [0.5, 0.3, 0.2, 0.0]   # All ratios defined, sum=1.0

MOSART Processing:
ValidateAndAdjustRatios INFO: 3 downstream connections with 3 ratios found for point 78901
- ratios sum to 1.000000

Result: No changes needed
```

### Scenario 3: No Ratios Defined

**Condition:** All ratios are zero or missing  
**Behavior:** ✅ **WARNING** - Sets primary downstream to 1.0, others to 0.0

**Example:**
```text
Input NetCDF Data:
dnID(cell_999, :) = [5001, 5002, -999, -999]        # 2 downstream connections
bifurc_ratio(cell_999, :) = [0.0, 0.0, 0.0, 0.0]   # No ratios defined

MOSART Processing:
ValidateAndAdjustRatios WARNING: No ratios defined for bifurcation point 99901
- setting primary downstream to 1.0

Final Ratios:
bifurc_ratio(cell_999, :) = [1.0, 0.0, 0.0, 0.0]   # All flow to primary
```

### Scenario 4: Perfect Ratio Input

**Condition:** All downstream connections have ratios that sum to 1.0  
**Behavior:** ✅ **PASS** - Standard validation, no messages

**Example:**
```text
Input NetCDF Data:
dnID(cell_111, :) = [5001, 5002, -999, -999]        # 2 downstream connections
bifurc_ratio(cell_111, :) = [0.6, 0.4, 0.0, 0.0]   # Perfect sum=1.0

MOSART Processing:
(No special messages - passes standard validation)

Result: Ratios used as-is
```

### Scenario 5: More Ratios Than Connections

**Condition:** Defined ratios exceed number of downstream connections  
**Behavior:** ❌ **ERROR** - Logic error in data structure

**Example:**
```text
ValidateAndAdjustRatios ERROR: 2 downstream connections but 3 ratios found for point 11122
```

**Note:** This scenario should not occur with proper NetCDF data structure but is checked for robustness.

### Practical Usage Examples

#### Example 1: Simple Delta with Unknown Primary Split
**Use Case:** Know secondary channel takes 40%, don't want to calculate primary

```text
# NetCDF Input:
dnID(delta, :) = [bay_outlet, south_channel, -999, -999]
bifurc_ratio(delta, :) = [0.0, 0.4, 0.0, 0.0]      # Only specify secondary

# MOSART Result:
bifurc_ratio(delta, :) = [0.6, 0.4, 0.0, 0.0]      # Primary gets remaining 60%
# Warning message logged
```

#### Example 2: Triple Bifurcation with Two Known Splits
**Use Case:** Complex delta with engineering studies for two channels

```text
# NetCDF Input:
dnID(complex_delta, :) = [main_bay, east_channel, west_channel, -999]
bifurc_ratio(complex_delta, :) = [0.0, 0.25, 0.35, 0.0]  # Know east=25%, west=35%

# MOSART Result:
bifurc_ratio(complex_delta, :) = [0.4, 0.25, 0.35, 0.0]  # Main gets remaining 40%
# Warning message logged
```

#### Example 3: Irrigation Canal with Fixed Diversion
**Use Case:** Agricultural diversion with known withdrawal percentage

```text
# NetCDF Input:
dnID(irrigation, :) = [river_continue, farm_canal, -999, -999]  
bifurc_ratio(irrigation, :) = [0.0, 0.15, 0.0, 0.0]     # Farm takes 15%

# MOSART Result:
bifurc_ratio(irrigation, :) = [0.85, 0.15, 0.0, 0.0]    # River keeps 85%
# Warning message logged
```

### Implementation Benefits

#### 1. **User-Friendly Input Preparation**
- **Before:** Users had to manually calculate `primary_ratio = 1.0 - sum(other_ratios)`
- **After:** Users specify only known/measured ratios, MOSART calculates remainder
- **Reduces Errors:** Eliminates manual calculation mistakes in input files

#### 2. **Flexible Data Integration**
- **Engineering Studies:** Can directly use measured secondary channel ratios
- **Observational Data:** Use field-measured split percentages without conversion
- **Iterative Modeling:** Easy to test different secondary channel scenarios

#### 3. **Robust Error Detection**
- **Catches Over-allocation:** Prevents ratios summing >1.0 (physically impossible)
- **Clear Messages:** Informative warnings explain exactly what adjustments were made
- **Conservation Guarantee:** All scenarios maintain perfect water balance

#### 4. **Backward Compatibility**
- **Existing Files:** No changes needed for current ratio input files
- **Standard Validation:** Files with complete ratios work exactly as before
- **Migration Path:** Users can gradually adopt partial ratio approach

### Configuration Recommendations

#### When to Use Partial Ratios:
- ✅ **Delta Studies:** Secondary channel ratios from field measurements
- ✅ **Irrigation Systems:** Known diversion percentages from water rights
- ✅ **Observational Modeling:** Incorporating measured flow splits
- ✅ **Sensitivity Analysis:** Testing different secondary channel scenarios

#### When to Use Complete Ratios:
- ✅ **Theoretical Studies:** All splits calculated from hydraulic equations
- ✅ **Balanced Deltas:** Equal splits or predetermined proportions
- ✅ **Model Verification:** Exact control over all flow distributions

### Technical Implementation Details

#### Validation Threshold: `1.e-12_r8`
- **Purpose:** Distinguish between "zero" and "very small" ratios
- **Robustness:** Handles floating-point precision issues
- **Conservative:** Treats near-zero values as intentionally unspecified

#### Tolerance for Sum Checking: `1.e-6_r8`
- **Purpose:** Allow for minor floating-point rounding in complete ratio sets
- **Standard:** Consistent with existing MOSART validation tolerances
- **Practical:** Accommodates typical NetCDF precision limitations

#### Memory Efficiency:
- **No Additional Storage:** Uses existing `bifurc_ratio` arrays
- **In-Place Modification:** Adjusts ratios during validation phase
- **Zero Overhead:** No runtime performance impact after initialization

## Bifurcation Output Variables

### BIFUR_AMOUNT - Water Transfer Tracking
**Purpose:** Track water diversions and transfers through the bifurcation/IBT system

**Variable Definition:**
- **Name:** `BIFUR_AMOUNT_LIQ` (liquid water tracer)
- **Units:** m³/s
- **Type:** History output variable (default='inactive')
- **Data Structure:** `rtmCTL%bifur_amount(:,:)` - 2D array for multiple tracers

**Value Interpretation:**
- **Negative values:** Water diverted away from natural flow at bifurcation/IBT points
- **Positive values:** Water received from diversions at destination cells
- **Zero values:** 
  - Non-bifurcation cells (normal flow routing)
  - Primary downstream cells (natural flow continuation, not tracked as "diversion")

**Water Balance:** Global sum equals zero, ensuring perfect conservation

**Examples:**

**50-50 Natural Bifurcation:**
```text
Bifurcation point:     BIFUR_AMOUNT = -50 m³/s  (50% diverted to secondary channel)
Primary downstream:    BIFUR_AMOUNT = 0 m³/s    (natural flow, not tracked)
Secondary downstream:  BIFUR_AMOUNT = +50 m³/s  (receives diverted water)
```

**IBT Diversion (20 m³/s):**
```text
Diversion point:       BIFUR_AMOUNT = -20 m³/s  (water taken for IBT)
Main river:           BIFUR_AMOUNT = 0 m³/s    (remainder flows naturally)
Destination basin:    BIFUR_AMOUNT = +20 m³/s  (receives transferred water)
```

## Grid Dimension Handling

**IMPORTANT:** The variable names `rtmlon` and `rtmlat` are misleading and do NOT always represent longitude/latitude dimensions.

### How `ncd_inqfdims` Works:
- **Purpose:** Auto-detects grid structure and returns dimensions through subroutine parameters
- **For Structured Grids (`isgrid2d = .true.`):**
  - `rtmlon` = longitude dimension (e.g., 144)
  - `rtmlat` = latitude dimension (e.g., 96)
  - NetCDF dimensions: `lon`, `lat` or `ni`, `nj`

- **For Unstructured Grids (`isgrid2d = .false.`):**
  - `rtmlon` = total number of grid cells (e.g., 50,000)
  - `rtmlat` = 1 (forced to 1, acts as dummy dimension)
  - NetCDF dimensions: `gridcell` (1D array)

### NetCDF Reading Implementation:
Both grid types use same loop structure:
```fortran
do j=1,rtmlat  ! For unstructured: j=1 only
do i=1,rtmlon  ! For unstructured: i=1 to number_of_gridcells
   n = (j-1)*rtmlon + i  ! Maps to 1D global index
   dnID_global(n,k) = itempr(i,j)
end do
end do
```

## WRM Unstructured Mesh Support

### Overview
The Water Resource Management (WRM) module supports both structured and unstructured mesh reservoir parameter files, providing the same grid flexibility as the main MOSART parameter file reading.

### Grid Type Detection
- **Independent Detection:** WRM files are analyzed independently of main MOSART parameter file grid type
- **Automatic Detection:** Uses `ncd_inqfdims()` to determine if reservoir file is structured or unstructured
- **Variable:** `iswm2d` determines WRM input file grid type (similar to `isgrid2d` for main MOSART)

### Supported File Formats

| **Grid Type** | **NetCDF Structure** | **Example** | **Dimensions** |
|---------------|---------------------|-------------|----------------|
| **Structured** | `DamInd_2d(lon, lat)` | Traditional lat-lon grid | `ndims = 2` |
| **Unstructured** | `DamInd_2d(gridcell)` | MPAS-style unstructured | `ndims = 1` |

### Usage Examples

**Unstructured Reservoir File:**
```text
dimensions:
    gridcell = 18559 ;
    Dams = 38 ;
    Month = 12 ;

variables:
    int DamInd_2d(gridcell) ;  // 1D unstructured
        DamInd_2d:long_name = "Dam index in the 2d global domain" ;
```

**Structured Reservoir File:**
```text
dimensions:
    lon = 144 ;
    lat = 96 ;
    Dams = 38 ;

variables:
    int DamInd_2d(lon, lat) ;  // 2D structured
        DamInd_2d:long_name = "Dam index in the 2d global domain" ;
```

## Mixed Mode Benefits - Key Innovation

**Unified Architecture Achievement:**
The mixed mode implementation represents a significant advancement that elegantly unifies natural and managed water processes:

1. **✅ True Heterogeneous Basin Support**
   - Same river system can have both natural delta bifurcation AND managed diversions
   - Each bifurcation point operates according to its appropriate physics
   - No artificial separation between "natural" and "human" water management

2. **✅ Unified Technical Infrastructure**  
   - Single MCT matrix system handles all bifurcation types
   - IBT demands → dynamic ratios → existing proven bifurcation logic
   - Leverages robust, tested MCT sparse matrix infrastructure
   - No duplicate code paths or separate physics implementations

3. **✅ Operational Flexibility**
   - Natural deltas: Fixed proportional splits regardless of flow conditions
   - Municipal/industrial: Absolute demands with environmental protection
   - Water rights: Legal allocations with automatic curtailment during droughts
   - Inter-basin transfers: Large-scale water movement between distant basins

4. **✅ Real-World Applicability**
   - California Central Valley: Delta splits + State Water Project + agricultural diversions
   - Mississippi River: Natural distributaries + municipal intakes + industrial cooling
   - Colorado River: Glen Canyon releases + municipal diversions + agricultural allocations
   - Any major river system with both natural and managed flow splitting

5. **✅ Robust Environmental Protection**
   - 10% minimum flow automatically enforced for all IBT points
   - Fixed ratio splits maintain natural flow patterns at deltas
   - Dynamic demand curtailment during low-flow periods
   - Perfect water balance conservation across all splitting types

**Technical Elegance:**
The key insight was recognizing that both fixed ratios and dynamic demands ultimately need to be converted to ratios for the MCT system. By implementing this conversion in the physics module, we achieved true mixed mode support without duplicating the complex MCT infrastructure.

## Scientific Applications

- **Seasonal Delta Dynamics:** Flow splitting ratios change based on discharge variability
- **Freshwater Delivery:** Track how seasonal flows affect ocean freshwater inputs
- **Sediment Transport:** Same framework applicable to sediment bifurcation research
- **Climate Impacts:** Investigate how changing hydrology affects delta function

## Future Extensions

- **Seasonal Demands:** Monthly varying diversion schedules
- **Reservoir Integration:** Coordinate with WRM module for storage releases
- **Return Flows:** Model return flow from diversions back to river system
- **Water Quality:** Track pollutant concentrations through diversions
- **Economic Optimization:** Allocate water based on economic value
- **Climate Adaptation:** Dynamic demands based on drought indices