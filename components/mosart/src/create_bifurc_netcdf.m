%function create_bifurc_netcdf()
% Create MOSART bifurcation NetCDF file with downstream dimension
% Supports both structured (2D->3D) and unstructured (1D->2D) mesh formats
% Adds bifurcation_ratio and ibt_demand variables

% Input and output files
pre = '\\compy01';
%input_file = [pre '/compyfs/zhou014/datasets/E3SM_inputs/MOSART_Global_half_20210422.nc'];
%output_file = [pre '/compyfs/zhou014/datasets/E3SM_inputs/MOSART_Global_half_20210422_bifurc2.nc'];

input_file = [pre '/compyfs/zhou014/ICoM/dataset/MOSART_Mid-Atlantic_MPAS_c220107.nc'];
output_file = [pre '/compyfs/zhou014/ICoM/dataset/MOSART_Mid-Atlantic_MPAS_c220107_4.nc'];

fprintf('Reading original NetCDF file: %s\n', input_file);

% Detect mesh type by checking cell ID dimensions
ID_info = ncinfo(input_file, 'ID');
ndims = length(ID_info.Size);
is_structured = (ndims == 2);

% Get the lat long for each cell
lat = ncread(input_file,'latixy');
lon = ncread(input_file,'longxy');
ID = ncread(input_file,'ID');

if is_structured
    % Structured mesh: dnID(lon(row), lat(col))
    dnID_2d = ncread(input_file, 'dnID');
    [nrow, ncol] = size(dnID_2d);
    ncells = nrow * ncol;
    fprintf('Structured mesh detected: dnID dimensions %d x %d (total %d cells)\n', nrow, ncol, ncells);
else
    % Unstructured mesh: dnID(cells)
    dnID_1d = ncread(input_file, 'dnID');
    ncells = length(dnID_1d);
    nrow = ncells;  % For unstructured, treat as 1D array
    ncol = 1;       % Dummy dimension
    fprintf('Unstructured mesh detected: dnID dimensions %d cells\n', ncells);
end

% Create arrays with downstream dimension
max_downstream = 2;

if is_structured
    % Structured: Create 3D arrays
    dnID_2d_expand_to_3d = -9999 * ones(nrow, ncol, max_downstream);
    bifurc_ratio_3d = zeros(nrow, ncol, max_downstream);
    ibt_demand_3d = zeros(nrow, ncol, max_downstream);

    % Copy original 2D data to first layer (primary downstream)
    dnID_2d_expand_to_3d(:, :, 1) = dnID_2d;

else
    % Unstructured: Create 2D arrays
    dnID_1d_expand_to_2d = -9999 * ones(ncells, max_downstream);
    bifurc_ratio_2d = zeros(ncells, max_downstream);
    ibt_demand_2d = zeros(ncells, max_downstream);

    % Copy original 1D data to first column (primary downstream)
    dnID_1d_expand_to_2d(:, 1) = dnID_1d;
end

% Copy original file to new file
fprintf('Copying original file to: %s\n', output_file);
copyfile(input_file, output_file);

% Open the copied file for modification
ncid = netcdf.open(output_file, 'WRITE');

try
    % Enter define mode to modify structure
    netcdf.reDef(ncid);

    % Create downstream dimension
    downstream_dimid = netcdf.defDim(ncid, 'downstream', max_downstream);

    % Get original dnID variable info
    dnID_varid = netcdf.inqVarID(ncid, 'dnID');
    [~, dnID_type, ~, ~] = netcdf.inqVar(ncid, dnID_varid);

    if is_structured
        % Structured mesh: Find existing dimension IDs
        lon_dimid = netcdf.inqDimID(ncid, 'lon');
        lat_dimid = netcdf.inqDimID(ncid, 'lat');

        % Rename original dnID variable to keep it
        netcdf.renameVar(ncid, dnID_varid, 'dnID_original');

        % Create new 3D variables: dnID(lon, lat, downstream)
        dnID_new_varid = netcdf.defVar(ncid, 'dnID', dnID_type, ...
            [lon_dimid, lat_dimid, downstream_dimid]);
        bifurc_ratio_varid = netcdf.defVar(ncid, 'bifurc_ratio', 'double', ...
            [lon_dimid, lat_dimid, downstream_dimid]);
        ibt_demand_varid = netcdf.defVar(ncid, 'ibt_demand', 'double', ...
            [lon_dimid, lat_dimid, downstream_dimid]);

        % Add variable attributes
        netcdf.putAtt(ncid, bifurc_ratio_varid, 'long_name', 'Bifurcation split ratios');
        netcdf.putAtt(ncid, bifurc_ratio_varid, 'units', 'dimensionless');
        netcdf.putAtt(ncid, bifurc_ratio_varid, 'valid_range', [0.0, 1.0]);

        netcdf.putAtt(ncid, ibt_demand_varid, 'long_name', 'Inter-basin transfer demand');
        netcdf.putAtt(ncid, ibt_demand_varid, 'units', 'm3/s');
        netcdf.putAtt(ncid, ibt_demand_varid, 'valid_min', 0.0);

    else
        % Unstructured mesh: Find existing dimension ID
        cells_dimid = netcdf.inqDimID(ncid, 'gridcell');  % Standard MPAS name

        % Rename original dnID variable to keep it
        netcdf.renameVar(ncid, dnID_varid, 'dnID_original');

        % Create new 2D variables: dnID(cells, downstream)
        dnID_new_varid = netcdf.defVar(ncid, 'dnID', dnID_type, ...
            [cells_dimid, downstream_dimid]);
        bifurc_ratio_varid = netcdf.defVar(ncid, 'bifurc_ratio', 'double', ...
            [cells_dimid, downstream_dimid]);
        ibt_demand_varid = netcdf.defVar(ncid, 'ibt_demand', 'double', ...
            [cells_dimid, downstream_dimid]);

        % Add variable attributes
        netcdf.putAtt(ncid, bifurc_ratio_varid, 'long_name', 'Bifurcation split ratios');
        netcdf.putAtt(ncid, bifurc_ratio_varid, 'units', 'dimensionless');
        netcdf.putAtt(ncid, bifurc_ratio_varid, 'valid_range', [0.0, 1.0]);

        netcdf.putAtt(ncid, ibt_demand_varid, 'long_name', 'Inter-basin transfer demand');
        netcdf.putAtt(ncid, ibt_demand_varid, 'units', 'm3/s');
        netcdf.putAtt(ncid, ibt_demand_varid, 'valid_min', 0.0);
    end

    % End define mode
    netcdf.endDef(ncid);

    % Write the data
    if is_structured
        fprintf('Writing 3D structured mesh data...\n');
        netcdf.putVar(ncid, dnID_new_varid, dnID_2d_expand_to_3d);
        netcdf.putVar(ncid, bifurc_ratio_varid, bifurc_ratio_3d);
        netcdf.putVar(ncid, ibt_demand_varid, ibt_demand_3d);
        fprintf('Structured format: dnID(lon,lat,downstream), bifurc_ratio(lon,lat,downstream), ibt_demand(lon,lat,downstream)\n');
    else
        fprintf('Writing 2D unstructured mesh data...\n');
        netcdf.putVar(ncid, dnID_new_varid, dnID_1d_expand_to_2d);
        netcdf.putVar(ncid, bifurc_ratio_varid, bifurc_ratio_2d);
        netcdf.putVar(ncid, ibt_demand_varid, ibt_demand_2d);
        fprintf('Unstructured format: dnID(cells,downstream), bifurc_ratio(cells,downstream), ibt_demand(cells,downstream)\n');
    end

    fprintf('Successfully created bifurcation NetCDF file!\n');
    fprintf('Note: Original dnID preserved as "dnID_original"\n');

catch ME
    fprintf('Error: %s\n', ME.message);
    netcdf.close(ncid);
    rethrow(ME);
end

% Close the file
netcdf.close(ncid);

%%
% Add values to the new file
% read in the bifurction/IBT file

csv_file = '\\connie-1\qpeople/zhou014/E3SMv3/code/MOSART_bifurcation/components/mosart/src\DRB_BIF.csv';
% csv_file = '/qfs/people/zhou014/E3SMv3/code/MOSART_bifurcation/components/mosart/src/DRB_IBT.csv';
bif = readtable(csv_file);
fprintf('Processing CSV file: %s\n', csv_file);
fprintf('Found %d bifurcation/IBT entries\n', height(bif));

% Initialize arrays for cell IDs to be filled in CSV
cell_ids = zeros(height(bif), 1);
dncell1_ids = zeros(height(bif), 1); 
dncell2_ids = zeros(height(bif), 1);

% Process each row in the CSV
for i = 1:height(bif)
    fprintf('\nProcessing row %d: %s\n', i, bif.note{i});
    
    % Extract data from CSV row
    split_lat = bif.cell_lat(i);
    split_lon = bif.cell_long(i);
    is_ibt = bif.ibt(i);
    
    % Convert longitude from 0-360 to -180-180 format if needed
   % if split_lon > 180
   %     split_lon = split_lon - 360;
   % end
    
    % Find closest cell for splitting point
    [split_cell_idx, split_cell_id] = find_closest_cell(lat, lon, ID, split_lat, split_lon, is_structured);
    if split_cell_idx == -1
        error('Could not find cell for splitting point at lat=%.4f, lon=%.4f', split_lat, split_lon);
    end
    
    fprintf('  Split point: lat=%.4f, lon=%.4f -> cell_id=%d\n', split_lat, split_lon, split_cell_id);
    cell_ids(i) = split_cell_id;
    
    % Handle primary downstream cell
    if bif.dncell1_lat(i) ~= -9999 && bif.dncell1_long(i) ~= -9999
        % Custom primary downstream specified
        dn1_lat = bif.dncell1_lat(i);
        dn1_lon = bif.dncell1_long(i);

     %   if dn1_lon > 180
     %       dn1_lon = dn1_lon - 360;
     %   end
        
        [dn1_cell_idx, dn1_cell_id] = find_closest_cell(lat, lon, ID, dn1_lat, dn1_lon, is_structured);
        if dn1_cell_idx == -1
            error('Could not find primary downstream cell at lat=%.4f, lon=%.4f', dn1_lat, dn1_lon);
        end
        
        fprintf('  Primary downstream: lat=%.4f, lon=%.4f -> cell_id=%d\n', dn1_lat, dn1_lon, dn1_cell_id);
        dncell1_ids(i) = dn1_cell_id;
        
        % Update primary downstream in arrays
        if is_structured
            [split_row, split_col] = get_structured_indices(split_cell_idx, nrow, ncol);
            dnID_2d_expand_to_3d(split_row, split_col, 1) = dn1_cell_id;
        else
            dnID_1d_expand_to_2d(split_cell_idx, 1) = dn1_cell_id;
        end
    else
        % Use existing primary downstream
        if is_structured
            [split_row, split_col] = get_structured_indices(split_cell_idx, nrow, ncol);
            existing_dnID = dnID_2d_expand_to_3d(split_row, split_col, 1);
        else
            existing_dnID = dnID_1d_expand_to_2d(split_cell_idx, 1);
        end
        
        dncell1_ids(i) = existing_dnID;
        fprintf('  Primary downstream: using existing cell_id=%d\n', existing_dnID);
    end
    
    % Handle secondary downstream cell
    dn2_lat = bif.dncell2_lat(i);
    dn2_lon = bif.dncell2_long(i);
  %  if dn2_lon > 180
  %      dn2_lon = dn2_lon - 360;
  %  end
    
    [dn2_cell_idx, dn2_cell_id] = find_closest_cell(lat, lon, ID, dn2_lat, dn2_lon, is_structured);
    if dn2_cell_idx == -1
        error('Could not find secondary downstream cell at lat=%.4f, lon=%.4f', dn2_lat, dn2_lon);
    end
    
    fprintf('  Secondary downstream: lat=%.4f, lon=%.4f -> cell_id=%d\n', dn2_lat, dn2_lon, dn2_cell_id);
    dncell2_ids(i) = dn2_cell_id;
    
    % Set secondary downstream connection and ratios/demands
    if is_structured
        [split_row, split_col] = get_structured_indices(split_cell_idx, nrow, ncol);
        dnID_2d_expand_to_3d(split_row, split_col, 2) = dn2_cell_id;
        
        if is_ibt == 1
            % IBT mode: use demand values
            ibt_demand_3d(split_row, split_col, 1) = 0.0;  % Primary gets remainder
            demand2_val = bif.demand2(i);
            if demand2_val == -9999
                demand2_val = 0.0;  % Treat missing as 0
            end
            ibt_demand_3d(split_row, split_col, 2) = demand2_val;
            fprintf('  IBT demand: %.2f m³/s\n', demand2_val);
        else
            % Regular bifurcation: use ratio values
            ratio1 = bif.ratio1(i);
            ratio2 = bif.ratio2(i);
            
            if ratio1 == -9999
                ratio1 = 0.0;
            end
            if ratio2 == -9999
                ratio2 = 0.0;
            end
            
            % Validate ratios
            if ratio1 + ratio2 > 1.0001  % Allow small floating point tolerance
                error('Ratios sum to %.4f > 1.0 for row %d', ratio1 + ratio2, i);
            end
            
            bifurc_ratio_3d(split_row, split_col, 1) = ratio1;
            bifurc_ratio_3d(split_row, split_col, 2) = ratio2;
            fprintf('  Bifurcation ratios: %.3f, %.3f\n', ratio1, ratio2);
        end
    else
        dnID_1d_expand_to_2d(split_cell_idx, 2) = dn2_cell_id;
        
        if is_ibt == 1
            % IBT mode: use demand values
            ibt_demand_2d(split_cell_idx, 1) = 0.0;  % Primary gets remainder
            demand2_val = bif.demand2(i);
            if demand2_val == -9999
                demand2_val = 0.0;  % Treat missing as 0
            end
            ibt_demand_2d(split_cell_idx, 2) = demand2_val;
            fprintf('  IBT demand: %.2f m³/s\n', demand2_val);
        else
            % Regular bifurcation: use ratio values
            ratio1 = bif.ratio1(i);
            ratio2 = bif.ratio2(i);
            
            if ratio1 == -9999
                ratio1 = 0.0;
            end
            if ratio2 == -9999
                ratio2 = 0.0;
            end
            
            % Validate ratios
            if ratio1 + ratio2 > 1.0001  % Allow small floating point tolerance
                error('Ratios sum to %.4f > 1.0 for row %d', ratio1 + ratio2, i);
            end
            
            bifurc_ratio_2d(split_cell_idx, 1) = ratio1;
            bifurc_ratio_2d(split_cell_idx, 2) = ratio2;
            fprintf('  Bifurcation ratios: %.3f, %.3f\n', ratio1, ratio2);
        end
    end
end

% Write the updated arrays to NetCDF file
fprintf('\nWriting updated arrays to NetCDF file...\n');
if is_structured
    ncwrite(output_file,'bifurc_ratio',bifurc_ratio_3d);
    ncwrite(output_file,'ibt_demand',ibt_demand_3d);
    ncwrite(output_file,'dnID',dnID_2d_expand_to_3d);
else
    ncwrite(output_file,'bifurc_ratio',bifurc_ratio_2d);
    ncwrite(output_file,'ibt_demand',ibt_demand_2d);
    ncwrite(output_file,'dnID',dnID_1d_expand_to_2d);
end

% Fill in the cell IDs in the CSV table
bif.cell_id = cell_ids;
bif.dncell1_id = dncell1_ids;
bif.dncell2_id = dncell2_ids;

% Write updated CSV file
output_csv = strrep(csv_file, '.csv', '_updated.csv');
writetable(bif, output_csv);
fprintf('Updated CSV file written to: %s\n', output_csv);

%%
% Verify the result
fprintf('\nVerifying result...\n');
verify_bifurc_file(output_file, height(bif), is_structured);

%end  % End of create_bifurc_netcdf function

%%

function [cell_idx, cell_id] = find_closest_cell(lat_array, lon_array, ID_array, target_lat, target_lon, is_structured)
% Find the closest cell to the target latitude and longitude
% Returns cell index and cell ID

min_distance = inf;
cell_idx = -1;
cell_id = -1;

if is_structured
    % For structured mesh, search through 2D arrays
    [nrow, ncol] = size(lat_array);
    for i = 1:nrow
        for j = 1:ncol
            % Calculate distance using Haversine formula approximation
            lat_diff = lat_array(i, j) - target_lat;
            lon_diff = lon_array(i, j) - target_lon;
            distance = sqrt(lat_diff^2 + lon_diff^2);
            
            if distance < min_distance
                min_distance = distance;
                cell_idx = (i-1)*ncol + j;  % Convert to 1D index
                cell_id = ID_array(i, j);
            end
        end
    end
else
    % For unstructured mesh, search through 1D arrays
    ncells = length(lat_array);
    for i = 1:ncells
        % Calculate distance
        lat_diff = lat_array(i) - target_lat;
        lon_diff = lon_array(i) - target_lon;
        distance = sqrt(lat_diff^2 + lon_diff^2);
        
        if distance < min_distance
            min_distance = distance;
            cell_idx = i;
            cell_id = ID_array(i);
        end
    end
end

% Check if we found a reasonable match (within ~1 degree)
if min_distance > 1.0
    fprintf('Warning: Closest cell is %.4f degrees away from target (lat=%.4f, lon=%.4f)\n', ...
        min_distance, target_lat, target_lon);
end

fprintf('    Found closest cell: distance=%.4f degrees, cell_id=%d\n', min_distance, cell_id);
end

function [row, col] = get_structured_indices(cell_idx, nrow, ncol)
% Convert 1D cell index to 2D row/col indices for structured mesh
col = mod(cell_idx - 1, ncol) + 1;
row = floor((cell_idx - 1) / ncol) + 1;
end

function verify_bifurc_file(filename, num_bifurc_points, is_structured)
% Verify the created file has correct structure and data

fprintf('Verification of %s:\n', filename);

% Read dimensions
info = ncinfo(filename);
for i = 1:length(info.Dimensions)
    if strcmp(info.Dimensions(i).Name, 'downstream')
        fprintf('  downstream dimension: %d\n', info.Dimensions(i).Length);
    end
end

% Read dnID and check structure
dnID_data = ncread(filename, 'dnID');
bifurc_ratio = ncread(filename, 'bifurc_ratio');
ibt_demand = ncread(filename, 'ibt_demand');

if is_structured
    % Structured mesh verification
    [nlon, nlat, ndownstream] = size(dnID_data);
    fprintf('  Structured mesh - dnID dimensions: %d x %d x %d\n', nlon, nlat, ndownstream);

    % Count non-missing values in each layer
    for k = 1:ndownstream
        layer = dnID_data(:, :, k);
        valid_count = sum(layer(:) > 0 & layer(:) ~= -9999);
        fprintf('  Layer %d: %d valid connections\n', k, valid_count);
    end

    % Count bifurcation points (cells with secondary downstream connections)
    layer2 = dnID_data(:, :, 2);
    bifurc_count = sum(layer2(:) > 0 & layer2(:) ~= -9999);
    fprintf('  Bifurcation points found: %d\n', bifurc_count);

else
    % Unstructured mesh verification
    [ncells, ndownstream] = size(dnID_data);
    fprintf('  Unstructured mesh - dnID dimensions: %d x %d\n', ncells, ndownstream);

    % Count non-missing values in each column
    for k = 1:ndownstream
        column = dnID_data(:, k);
        valid_count = sum(column > 0 & column ~= -9999);
        fprintf('  Column %d: %d valid connections\n', k, valid_count);
    end

    % Count bifurcation points (cells with secondary downstream connections)
    column2 = dnID_data(:, 2);
    bifurc_count = sum(column2 > 0 & column2 ~= -9999);
    fprintf('  Bifurcation points found: %d\n', bifurc_count);
end

% Check if we have the expected number of bifurcation points
if bifurc_count == num_bifurc_points
    fprintf('  ✓ Expected %d bifurcation points, found %d\n', num_bifurc_points, bifurc_count);
else
    fprintf('  ✗ Expected %d bifurcation points, found %d\n', num_bifurc_points, bifurc_count);
end

% Count non-zero ratio and demand values
nonzero_ratios = sum(bifurc_ratio(:) > 0);
nonzero_demands = sum(ibt_demand(:) > 0);
fprintf('  Non-zero bifurc_ratio values: %d\n', nonzero_ratios);
fprintf('  Non-zero ibt_demand values: %d\n', nonzero_demands);

% Summary of variables created
fprintf('  Variables created:\n');
fprintf('    - dnID: Downstream connectivity with bifurcation support\n');
fprintf('    - bifurc_ratio: Split ratios (0.0-1.0, zeros for most cells)\n');
fprintf('    - ibt_demand: IBT demands (m³/s, zeros for most cells)\n');
fprintf('    - dnID_original: Backup of original dnID\n');

fprintf('Verification complete.\n');

end