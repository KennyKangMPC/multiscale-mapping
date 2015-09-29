%%% RBF contruction

function EM = constructrbf_F(lon,lat,srbfparam,loni,lati)
g = srbfparam(2);
maxN = srbfparam(1);

% define rbf form
rE = 6371*1000;

% rbf = @(e,r) exp(-(e*r).^2); % gaussian 
rbf = @(g,t) 1./(rE*(1 + g^2 - 2*g*cos(t)).^.5);  % spherical reciprocal MQ

% find centers via cluster analysis
X = [lon lat];
distfun = @(Xi,Xj)sph_distfun(Xi,Xj);
N = maxN;

T = clusterdata(X,'distance',distfun,'criterion','distance','maxclust',N,'linkage','average');

lonC = NaN(N,1); latC = NaN(N,1);
for i = 1:N
    
    Xmn = meanm(X(T == i,2),X(T == i,1));
    
    lonC(i) = Xmn(2);
    latC(i) = Xmn(1);
    
end

% calculate distance matrix between grid points and centers

lon1 = repmat(loni,[1 N]);
lon2 = repmat(lonC',[length(loni) 1]);
lat1 = repmat(lati,[1 N]);
lat2 = repmat(latC',[length(lati) 1]);

dm_grid = sph_dist_ca(lon1,lon2,lat1,lat2);

% evaluate rbf at grid points

EM = rbf(g,dm_grid/rE);

return
