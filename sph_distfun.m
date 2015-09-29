%%% sph_dist(XI,XJ,lat1,lat2) computes distance in m using the 
%%% Haversine Formula, assuming spherical earth of radius 6378.137km.  
%%% Distances are probably good to better than 1% of the "true" distance 
%%% on the ellipsoidal earth.
%%% XI = [lonI latI]; XJ = [XI lat1; XJ lat2; ...; lonM latM];


function val = sph_distfun(XI,XJ)

R = 6371*1000;

XI(XI(:,1) < 0,1) = XI(XI(:,1) < 0,1) + 360;
XJ(XJ(:,1) < 0,1) = XJ(XJ(:,1) < 0,1) + 360;

XI = XI*pi/180;
XJ = XJ*pi/180;

dlon = bsxfun(@minus,XI(1),XJ(:,1));
lat1 = repmat(XI(2),[length(XJ(:,2)) 1]);
lat2 = XJ(:,2);

ang = sin(lat1).*sin(lat2) + cos(lat1).*cos(lat2).*cos(dlon);
ang(ang > 1) = 1;
val = R*acos(ang);

return
