function calculate_MODIS_to_SMAP36km_regridding_matrix()

if isempty(gcp)
    parpool(20);
end

% The basic plan here is to create a sparse matrix that you can multiply
% against a vector of MODIS data to get SMAP-gridded data (i.e.,
% MODIS_on_SMAP_grid = MODIS_to_SMAP_mat * MODIS_data(:) )

% MODIS is on a 0.05 deg lat x 0.05 deg lon grid, with corners at -180 lon,
% 180 lon, -90 lat, and 90 lat.
MOD_centers_lat = (-89.975:0.05:89.975)';
MOD_centers_lon = -179.975:0.05:179.975;
MOD_centers_lat = repmat(MOD_centers_lat,1,length(MOD_centers_lon));
MOD_centers_lon = repmat(MOD_centers_lon,size(MOD_centers_lat,1),1);

% SMAP has the lats going from highest to lowest, so we're going to flip
% the M2 data to match:
MOD_centers_lat = flipud(MOD_centers_lat);
MOD_centers_lon = flipud(MOD_centers_lon);

% Now, on to SMAP!
% Use saved lat lon data if available:
if ~exist('SMAP_NW_corners36.mat','file')
    [cols_36km, rows_36km] = meshgrid((1:965)-0.5,(1:407)-0.5);
    [SMAP_corners_lat, SMAP_corners_lon] = EASE22LatLon_M36KM(rows_36km, cols_36km); % These are the NW corners
    if exist('SMAPlatlonP_no_nan.mat','file')
        load('SMAPlatlonP_no_nan.mat');
        SMAP_centers_lat = lat;
        SMAP_centers_lon = lon;
    else
        SMAP_centers_lat = (SMAP_corners_lat + circshift(SMAP_corners_lat,-1))/2;
        SMAP_centers_lon = (SMAP_corners_lon + circshift(SMAP_corners_lon,[0,-1]))/2;
        SMAP_centers_lat(end,:) = [];
        SMAP_centers_lat(:,end) = [];
        SMAP_centers_lon(end,:) = [];
        SMAP_centers_lon(:,end) = [];
    end
    save('SMAP_NW_corners36.mat','SMAP_corners_lat','SMAP_corners_lon',...
        'SMAP_centers_lat','SMAP_centers_lon');
else
    load('SMAP_NW_corners36.mat');
end

areas_vec = areaquad( reshape(SMAP_corners_lat(2:end,1:(end-1)),[numel(SMAP_centers_lat),1]),...
    reshape(SMAP_corners_lon(2:end,1:(end-1)),[numel(SMAP_centers_lat),1]),...
    reshape(SMAP_corners_lat(1:(end-1),2:end),[numel(SMAP_centers_lat),1]),...
    reshape(SMAP_corners_lon(1:(end-1),2:end),[numel(SMAP_centers_lat),1]),...
    referenceEllipsoid('wgs84','kilometers'),'degrees');
SMAP_areas = reshape(areas_vec,size(SMAP_centers_lat));


% Now, we need to loop over every SMAP grid cell, and get the area-averaged
% values that go into it

% Define a distance that will be sure to include all possibly overlapping
% cells around a SMAP centerpoint:
% max_SMAP_L = max(distance(SMAP_centers_lat(1:(end-1),1),SMAP_centers_lon(1:(end-1),1),...
%     SMAP_centers_lat(2:end,2),SMAP_centers_lon(2:end,2),referenceEllipsoid('wgs84','kilometers'),'degrees'));
% 
% max_M2_L = max(distance(M2_centers_lat(1:(end-1),1),M2_centers_lon(1:(end-1),1),...
%     M2_centers_lat(2:end,2),M2_centers_lon(2:end,2),referenceEllipsoid('wgs84','kilometers'),'degrees'));
% 
% r = 1*(max_SMAP_L + max_M2_L); % Anything closer than this should be checked to see if it overlaps with the SMAP box

max_MOD_lat_spacing = max(max(abs(diff(MOD_centers_lat,1,1))));
max_MOD_lon_spacing = max(max(abs(diff(MOD_centers_lon,1,2))));

max_SMAP_lat_spacing = max(max(abs(diff(SMAP_centers_lat,1,1))));
max_SMAP_lon_spacing = max(max(abs(diff(SMAP_centers_lon,1,2))));

r_lat = max_MOD_lat_spacing+max_SMAP_lat_spacing;
r_lon = max_MOD_lon_spacing+max_SMAP_lon_spacing;

% We're going to make a version of the lat and lon matrices that have some
% additional columns added on both sides: these columns will be duplicates
% so that instead of running from -180 to 180, it will be more like -200 to
% 200 lon. This will allow us to more easily find neigboring grid cells.
N_lon_cells_added_MOD = sum(MOD_centers_lon(1,:) > (180-2*r_lon)); % We'll add this number on to both sides of the M2 lon matrices

MOD_centers_lon_extended = [MOD_centers_lon(:, (end-N_lon_cells_added_MOD+1):(end) )-360, ...
    MOD_centers_lon, ...
    MOD_centers_lon(:, 1:(N_lon_cells_added_MOD) )+360];
MOD_centers_lat_extended = [MOD_centers_lat(:, (end-N_lon_cells_added_MOD+1):(end) ), ...
    MOD_centers_lat, ...
    MOD_centers_lat(:, 1:(N_lon_cells_added_MOD) )];

% SMAP_corners_lon360 = SMAP_corners_lon;
% SMAP_corners_lon360(SMAP_corners_lon360 < 0) = SMAP_corners_lon360(SMAP_corners_lon360 < 0) + 180;
% SMAP_centers_lon360 = SMAP_centers_lon;
% SMAP_centers_lon360(SMAP_centers_lon360 < 0) = SMAP_centers_lon360(SMAP_centers_lon360 < 0) + 180;
% 
% M2_centers_lon360 = M2_centers_lon;
% M2_centers_lon360(M2_centers_lon360 < 0) = M2_centers_lon360(M2_centers_lon360 < 0) + 180;

C = cell(numel(SMAP_centers_lat),1);

ref_elipsoid = referenceEllipsoid('wgs84','kilometers');

for i = 1:size(SMAP_centers_lat,1)
    fprintf('i = %d / %d...\n',i,size(SMAP_centers_lat,1));
    for j = 1:size(SMAP_centers_lat,2)
        fprintf('i = %d / %d, j = %d / %d...\n',i,size(SMAP_centers_lat,1),...
            j,size(SMAP_centers_lat,2));
        
        SMAP_poly_latlon = [SMAP_corners_lat(i,j),SMAP_corners_lon(i,j); ...
            SMAP_corners_lat(i,j+1),SMAP_corners_lon(i,j+1); ...
            SMAP_corners_lat(i+1,j+1),SMAP_corners_lon(i+1,j+1); ...
            SMAP_corners_lat(i+1,j),SMAP_corners_lon(i+1,j)]; % Going CW from NW corner
        
        % Find possible overlappers:
        % The distance command takes forever, so we'll give it an initial
        % speedup:
        close_MODs = abs(MOD_centers_lat_extended(:) - SMAP_centers_lat(i,j)) < 1.05*r_lat ...
            & abs(MOD_centers_lon_extended(:) - SMAP_centers_lon(i,j)) < 1.05*r_lon;
        
        % Make a matrix out of all of the close_M2s:
        MOD_tmp_lat = [MOD_centers_lat_extended(close_MODs)+.025, MOD_centers_lat_extended(close_MODs)+.025, ...
            MOD_centers_lat_extended(close_MODs)-.025, MOD_centers_lat_extended(close_MODs)-.025]'; % 4xN
        MOD_tmp_lon = [MOD_centers_lon_extended(close_MODs)-.025, MOD_centers_lon_extended(close_MODs)+.025, ...
            MOD_centers_lon_extended(close_MODs)+.025, MOD_centers_lon_extended(close_MODs)-.025]'; % 4xN
        
        
        
        MOD_indices = find(close_MODs);
        % For some really annoying reason, polybool won't let me do this
        % all at once, so another for loop it is:
        N = sum(close_MODs);
        
        [i_vec,j_vec] = ind2sub(size(MOD_centers_lat_extended),MOD_indices);
        a_vec = zeros(N,1); % The areas will go in here...
        for n = 1:N
            %[lon_intersect, lat_intersect] = polybool('intersection',SMAP_poly_latlon(:,2)+lon_offset,SMAP_poly_latlon(:,1),...
            %    M2_tmp_lon(:,n),M2_tmp_lat(:,n));
            [lon_intersect, lat_intersect] = my_polybool(SMAP_poly_latlon(:,2),SMAP_poly_latlon(:,1),...
                MOD_tmp_lon(:,n),MOD_tmp_lat(:,n));
            
            if ~isempty(lat_intersect)
                a_vec(n) = areaquad(min(lat_intersect),min(lon_intersect),...
                    max(lat_intersect),max(lon_intersect),ref_elipsoid,'degrees');
            end
        end
        
        % Should probably check here to make sure that the total overlap
        % area is close to the SMAP area:
        ratio = sum(a_vec)/SMAP_areas(i,j);
        if ratio <0.9 || ratio > 1.1
            error('The areas are not adding up!\n Aborting!');
        end
        
        % Normalize by total area:
        a_vec = a_vec/SMAP_areas(i,j);
        C{sub2ind(size(SMAP_centers_lat),i,j)} = [i_vec,j_vec,a_vec];
    end

end

delete(gcp);

% Now time to make a sparse matrix

save('tmp_areas_cells.mat','C');

% So the idea is that now we have a massive number of matrices, and we need
% to combine them into one giant sparse matrix that has numel(SMAP) rows
% and numel(M2) columns.
% As a first step, we need to turn the j_vec values (column 2) in all of
% the C{:} matrices into the proper values for the M2_centers matrices --
% right now, they are for the M2_extended matrices:

N = numel(C);
[nrowMOD, ncolMOD] = size(MOD_centers_lat);

n_sparse_elements = 0; % We're going to use this value to initialize the vectors of values going into the sparse matrix

% The M2_extended matrices have 2*N_lon_cells_added_M2 + ncolM2 columns
fprintf('Re-calculating proper column indices...\n');
for n = 1:N
    if mod(n,1000) == 0;
        fprintf('%d / %d completed...\n',n,N);
    end
    j_ext = C{n}(:,2);
    j_orig = zeros(size(j_ext));
    % Three cases: in the left add-on wing, in the middle, or in the right
    % add-on wing:
    % Case 1: Middle:
    mids = (j_ext > N_lon_cells_added_MOD) & (j_ext <= N_lon_cells_added_MOD+ncolMOD);
    j_orig(mids) = j_ext(mids) - N_lon_cells_added_MOD;
    
    % Case 2: Left wing:
    lefts = (j_ext <= N_lon_cells_added_MOD);
    j_orig(lefts) = j_ext(lefts) - N_lon_cells_added_MOD + ncolMOD;
    
    % Case 3: Right wing:
    rights = (j_ext > N_lon_cells_added_MOD+ncolMOD);
    j_orig(rights) = j_ext(rights) - N_lon_cells_added_MOD - ncolMOD;
    
    indxs = sub2ind( [nrowMOD, ncolMOD], C{n}(:,1) + 0, j_orig );
    C{n} = [indxs, C{n}(:,3) ]; % Now just 2 columns, and the indxs vectors correspond to the indices of the full MOD matrix
    n_sparse_elements = n_sparse_elements + length(indxs);
end

% So now we have a cell array of two-column matrices, where the first
% column is the index of the MOD matrix (which will correspond to a column
% number in our sparse matrix), and the second column is the area weight
% (i.e., the value that will go in the matrix). Sparse automatically adds
% together any elements that have duplicate row/column pairs, so we don't
% even need to make sure that the elements in C{n}(:,1) are unique

% We're going to loop through one more time and make our i_sparse, j_sparse
% and v_sparse vectors:
i_sparse = zeros(n_sparse_elements,1);
j_sparse = zeros(n_sparse_elements,1);
v_sparse = zeros(n_sparse_elements,1);

fprintf('Re-arranging into sparse matrix form...\n');

counter = 1;
for n = 1:N
    if mod(n,1000) == 0;
        fprintf('%d / %d completed...\n',n,N);
    end
    L = size(C{n},1);
    i_sparse(counter:(counter+L-1)) = n;
    j_sparse(counter:(counter+L-1)) = C{n}(:,1);
    v_sparse(counter:(counter+L-1)) = C{n}(:,2);
    counter = counter+L;
end

numel_SMAP = N;
numel_MOD = (nrowMOD+2)*ncolMOD;

fprintf('Creating sparse matrix...\n');
MOD_to_SMAP = sparse(i_sparse,j_sparse,v_sparse,numel_SMAP,numel_MOD);

save('MOD_to_SMAP.mat','MOD_to_SMAP');
fprintf('MOD_to_SMAP matrix saved in MOD_to_SMAP.mat\n');

end % function


function [lat, lon] = EASE22LatLon_M36KM(row, col)

%% WGS_1984
a = 6378137;            % Semimajor axis of ellipsoid
f = 1/298.257223563;    % Flattening
e = sqrt(2*f-f^2);      % Eccentricity

%% EASE-Grid 2.0
phi0    =  0.0;
lambda0 =  0.0;
phi1    = 30.0;

rows = 406;
cols = 964;
res  = 36032.220840584;

r0 = (cols + 1)/2;
s0 = (rows + 1)/2;

x = (col-r0) * res;
y = (s0-row) * res;
k0 = cosd(phi1)/sqrt(1-e^2*sind(phi1)^2);
qp = (1-e^2)*(sind(90)/(1-e^2*sind(90)^2) - 1/(2*e)*log((1-e*sind(90))/(1+e*sind(90))));
beta = asind(2*y*k0/a/qp);

% phi = beta + (e^2/3 + 31*e^4/180 + 517*e^6/5040)*sind(2*beta) ...
%            + (23*e^4/360 + 251*e^6/3780)*sind(4*beta) ...
%            + (761*e^6/45360)*sind(6*beta);

q = qp*sind(beta);

phi = asind(q/2);
% delta = phi;
% while any(delta >= 1e-15)
for ii = 1:2000
    phi_new = phi + ((1-e^2*sind(phi).^2).^2)./(2*cosd(phi)) .* ...
              (q/(1-e^2) - sind(phi)./(1-e^2*sind(phi).^2) + ...
              1/(2*e)*log((1-e*sind(phi))./(1+e*sind(phi))));
%     delta = abs(phi_new-phi);
    phi = phi_new;
end
lambda = lambda0 + x/a/k0;

lat = phi;
lon = lambda/pi*180;

return;

end



