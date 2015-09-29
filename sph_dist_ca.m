%%% sph_dist(lon1,lon2,lat1,lat2) computes distance between two points 
%%% on earths' surface in m using central angle formula

function val = sph_dist_ca(lon1,lon2,lat1,lat2)

R = 6371*1000;

lon1(lon1 < 0) = lon1(lon1 < 0) + 360;
lon2(lon2 < 0) = lon2(lon2 < 0) + 360;

dlon = (lon1 - lon2);

lon1 = lon1*pi/180;
lon2 = lon2*pi/180;
lat1 = lat1*pi/180;
lat2 = lat2*pi/180;
dlon = dlon*pi/180;

if (numel(lon1)==numel(lon2) && ((numel(lon2) == numel(lat1)) && numel(lat1) == numel(lat2)))

 ang = sin(lat1).*sin(lat2) + cos(lat1).*cos(lat2).*cos(dlon);
 ang(ang > 1) = 1;
 val = R*acos(ang);
else
    disp('Error - inputs must have equal dimensions');
end

return
