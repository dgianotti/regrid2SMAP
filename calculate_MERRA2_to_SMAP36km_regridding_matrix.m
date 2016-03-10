function calculate_MERRA2_to_SMAP36km_regridding_matrix()


% The basic plan here is to create a sparse matrix that you can multiply
% against a vector of MERRA2 data to get SMAP-gridded data (i.e.,
% M2_on_SMAP_grid = M2_to_SMAP_mat * MERRA2_data(:) )

% MERRA2 is on a 0.5 deg lat x 0.625 deg lon grid
M2_centers_lat = (-90:0.5:90)';
M2_centers_lon = -180:0.625:179.375;
M2_centers_lat = repmat(M2_centers_lat,1,length(M2_centers_lon));
M2_centers_lon = repmat(M2_centers_lon,size(M2_centers_lat,1),1);

% SMAP has the lats going from highest to lowest, so we're going to flip
% the M2 data to match:
M2_centers_lat = flipud(M2_centers_lat);
M2_centers_lon = flipud(M2_centers_lon);

% We're going to just throw out the northern-most and southern-most rows of
% data points (all of which are at the south or north poles... confusing
% how these can be considered "centers"... set them to NaN. Thus, the
% corners are easy to find as just being 1/2 of a lat-step and 1/2 of a
% lon-step off from the centers:
M2_centers_lat( [1,end],: ) = [];
M2_centers_lon( [1,end],: ) = [];

areas_vec = areaquad(M2_centers_lat(:)-0.25, M2_centers_lon(:)-0.3125,...
    M2_centers_lat(:)+0.25, M2_centers_lon(:)+0.3125,...
    referenceEllipsoid('wgs84','kilometers'),'degrees');
M2_areas = reshape(areas_vec,size(M2_centers_lat));

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
max_SMAP_L = max(distance(SMAP_centers_lat(1:(end-1),1),SMAP_centers_lon(1:(end-1),1),...
    SMAP_centers_lat(2:end,2),SMAP_centers_lon(2:end,2),referenceEllipsoid('wgs84','kilometers'),'degrees'));

max_M2_L = max(distance(M2_centers_lat(1:(end-1),1),M2_centers_lon(1:(end-1),1),...
    M2_centers_lat(2:end,2),M2_centers_lon(2:end,2),referenceEllipsoid('wgs84','kilometers'),'degrees'));

r = 1*(max_SMAP_L + max_M2_L); % Anything closer than this should be checked to see if it overlaps with the SMAP box

max_M2_lat_spacing = max(max(abs(diff(M2_centers_lat,1,1))));
max_M2_lon_spacing = max(max(abs(diff(M2_centers_lon,1,2))));

max_SMAP_lat_spacing = max(max(abs(diff(SMAP_centers_lat,1,1))));
max_SMAP_lon_spacing = max(max(abs(diff(SMAP_centers_lon,1,2))));

r_lat = max_M2_lat_spacing+max_SMAP_lat_spacing;
r_lon = max_M2_lon_spacing+max_SMAP_lon_spacing;

C = cell(numel(SMAP_centers_lat),1);

% Make a copy of the data that has 0-360 lon instead of -180-180 lon so
% that we cal look for neighbors more easily!
SMAP_corners_lon360 = SMAP_corners_lon;
SMAP_corners_lon360(SMAP_corners_lon360 < 0) = SMAP_corners_lon360(SMAP_corners_lon360 < 0) + 180;
SMAP_centers_lon360 = SMAP_centers_lon;
SMAP_centers_lon360(SMAP_centers_lon360 < 0) = SMAP_centers_lon360(SMAP_centers_lon360 < 0) + 180;
M2_centers_lon360 = M2_centers_lon;
M2_centers_lon360(M2_centers_lon360 < 0) = M2_centers_lon360(M2_centers_lon360 < 0) + 180;

for i = 1:size(SMAP_centers_lat,1)
    for j = 1:size(SMAP_centers_lat,2)
        fprintf('i = %d / %d, j = %d / %d...\n',i,size(SMAP_centers_lat,1),...
            j,size(SMAP_centers_lat,2));
        
        if abs(SMAP_centers_lon(i,j)) < 160 % If not near the dateline...
            SMAP_poly_latlon = [SMAP_corners_lat(i,j),SMAP_corners_lon(i,j); ...
                SMAP_corners_lat(i,j+1),SMAP_corners_lon(i,j+1); ...
                SMAP_corners_lat(i+1,j+1),SMAP_corners_lon(i+1,j+1); ...
                SMAP_corners_lat(i+1,j),SMAP_corners_lon(i+1,j)]; % Going CW from NW corner
            
            % Find possible overlappers:
            % The distance command takes forever, so we'll give it an initial
            % speedup:
            close_M2s = abs(M2_centers_lat(:) - SMAP_centers_lat(i,j)) < 1.05*r_lat ...
                & abs(M2_centers_lon(:) - SMAP_centers_lon(i,j)) < 1.05*r_lon;
                        
            % Make a matrix out of all of the close_M2s:
            M2_tmp_lat = [M2_centers_lat(close_M2s)+.25, M2_centers_lat(close_M2s)+.25, ...
                M2_centers_lat(close_M2s)-.25, M2_centers_lat(close_M2s)-.25]'; % 4xN
            M2_tmp_lon = [M2_centers_lon(close_M2s)-.3125, M2_centers_lon(close_M2s)+.3125, ...
                M2_centers_lon(close_M2s)+.3125, M2_centers_lon(close_M2s)-.3125]'; % 4xN
            
            
        else % If close to the dateline, use the 0-360 lon values!
            SMAP_poly_latlon = [SMAP_corners_lat(i,j),SMAP_corners_lon360(i,j); ...
                SMAP_corners_lat(i,j+1),SMAP_corners_lon360(i,j+1); ...
                SMAP_corners_lat(i+1,j+1),SMAP_corners_lon360(i+1,j+1); ...
                SMAP_corners_lat(i+1,j),SMAP_corners_lon360(i+1,j)]; % Going CW from NW corner
            
            % Find possible overlappers:
            % The distance command takes forever, so we'll give it an initial
            % speedup:
            close_M2s = abs(M2_centers_lat(:) - SMAP_centers_lat(i,j)) < 1.05*r_lat ...
                & abs(M2_centers_lon360(:) - SMAP_centers_lon360(i,j)) < 1.05*r_lon;
            
            % Make a matrix out of all of the close_M2s:
            M2_tmp_lat = [M2_centers_lat(close_M2s)+.25, M2_centers_lat(close_M2s)+.25, ...
                M2_centers_lat(close_M2s)-.25, M2_centers_lat(close_M2s)-.25]'; % 4xN
            M2_tmp_lon = [M2_centers_lon360(close_M2s)-.3125, M2_centers_lon360(close_M2s)+.3125, ...
                M2_centers_lon360(close_M2s)+.3125, M2_centers_lon360(close_M2s)-.3125]'; % 4xN
            
        end
        
        M2_indices = find(close_M2s);
        % For some really annoying reason, polybool won't let me do this
        % all at once, so another for loop it is:
        N = sum(close_M2s);
        
        [i_vec,j_vec] = ind2sub(size(M2_centers_lat),M2_indices);
        a_vec = zeros(N,1); % The areas will go in here...
        for n = 1:N
            %[lon_intersect, lat_intersect] = polybool('intersection',SMAP_poly_latlon(:,2)+lon_offset,SMAP_poly_latlon(:,1),...
            %    M2_tmp_lon(:,n),M2_tmp_lat(:,n));
            [lon_intersect, lat_intersect] = my_polybool(SMAP_poly_latlon(:,2),SMAP_poly_latlon(:,1),...
                M2_tmp_lon(:,n),M2_tmp_lat(:,n));
            
            if ~isempty(lat_intersect)
                a_vec(n) = areaquad(min(lat_intersect),min(lon_intersect),...
                    max(lat_intersect),max(lon_intersect),referenceEllipsoid('wgs84','kilometers'),'degrees');
            end
        end
        
        % Should probably check here to make sure that the total overlap
        % area is close to the SMAP area:
        ratio = sum(a_vec)/SMAP_areas(i,j);
        if ratio <0.9 || ratio > 1.1
            error('The areas are not adding up!\n Aborting!');
        end
        
        % Normalize by total area to force area matching:
        a_vec = a_vec/ratio;
        C{sub2ind(size(SMAP_centers_lat),i,j)} = [i_vec,j_vec,a_vec];
    end

end

% Now time to make a sparse matrix

i = 1;


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



