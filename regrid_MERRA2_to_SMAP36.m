function [X_on_SMAP_grid, SMAP_centers_lat,SMAP_centers_lon,...
    min_lat_SMAP,max_lat_SMAP] = regrid_MERRA2_to_SMAP36(X)
% The function regrid_MERRA2_to_SMAP36 regrids data from MERRA2 to the 36km
% EASEgrid2 grid used by the SMAP mission. MERRA2 data, when read from
% their standard NetCDF files is a 576(lon) x 361(lat) matrix, and
% regrid_MERRA2_to_SMAP expects the input X to be this size (or a 3D array
% with size(X,1) = 576, size(X,2) = 361.
%
%


if ~exist('M2_to_SMAP.mat','file')
    error('The file M2_to_SMAP.mat was not found. Is it in your path? Aborting!');
end

if size(X,1) == 361 && size(X,2) == 576
    warning('The input X to the function regrid_MERRA2_to_SMAP36 seems to be the transpose of what would be expected. Transposing and continuing...');
    X = permute(X,[2,1,3]);
end
    
% We want to be able to apply this to either 2D X matrices (which should be
% 361 x 576) or 3D arrays (which should be 361 x 576 x d):
if (size(X,1) ~= 576) || (size(X,2) ~= 361)
    error('The function regrid_MERRA2_to_SMAP36 expects the input X to be a 576x361 matrix or a 576x361xd array on the MERRA2 grid. Aborting!');
end

load('M2_to_SMAP.mat'); % Loads the M2_to_SMAP matrix

V = reshape(permute(X,[2,1,3]),[size(X,1)*size(X,2),size(X,3)]);
X_on_SMAP_grid = flipud(reshape(M2_to_SMAP * V, [406,964,size(X,3)]));

if nargout > 1
   if ~exist('SMAP36_centers_latlon.mat','file')
       error('The file SMAP36_centers_latlon.mat was not found. Is it in your path? Aborting!');
   end
   load('SMAP36_centers_latlon.mat');
   SMAP_centers_lat = SMAP36_centers_lat;
   SMAP_centers_lon = SMAP36_centers_lon;
end

% In case anyone wants the ranges of the max/min grid cells:
min_lat_SMAP = -85.0445664076212;
max_lat_SMAP = -min_lat_SMAP;

end % function