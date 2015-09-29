%%% RBF contruction

function [CM,L,g] = constructrbf_L(lon,lat,maxN)

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

% calculate min distance between centers

lon1 = repmat(lonC,[1 N]);
lat1 = repmat(latC,[1 N]);

dm_centers = sph_dist_ca(lon1,lon1',lat1,lat1');
dm_centers_sort = sort(dm_centers); 
L = min(dm_centers_sort(2,:))*10;

% rE = 6371;
gfun = @(g) ((1 + g.^2 - 2*g.*cos(L/rE))./((1 + g.^2 - 2*g)*exp(2)) - 1);
g = fzero(gfun,0.7);
L = L/1000;

% calculate distance matrix between data points and centers

lon1 = repmat(lon,[1 N]);
lon2 = repmat(lonC',[length(lon) 1]);
lat1 = repmat(lat,[1 N]);
lat2 = repmat(latC',[length(lat) 1]);

dm_data = sph_dist_ca(lon1,lon2,lat1,lat2);

% evaluate rbf at data points

CM = rbf(g,dm_data/rE);

return
