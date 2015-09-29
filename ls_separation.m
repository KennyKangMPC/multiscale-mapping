function ls_sep = ls_separation(prof)

timevec = datevec(prof.time);

maxind = 1;
for month = 1:12
    ind = find(timevec(:,2) == month);
    maxind = max([length(ind); maxind]);
end

ls_sep = struct;
ls_sep.dx = NaN(maxind^2,12);
ls_sep.dt = NaN(maxind^2,12);

for month = 1:12

    disp(month);
    ind = find(timevec(:,2) == month);
    lonm = prof.lon(ind);
    latm = prof.lat(ind);
    timem = prof.time(ind);

   dxi = NaN(length(lonm)); 
   dti = NaN(length(lonm)); 
    
for i = 1:length(lonm)
    dxi(:,i) = sph_dist_ca(repmat(lonm(i),[length(lonm) 1]),lonm,repmat(latm(i),[length(lonm) 1]),latm)/1000; % in km
    dti(:,i) = abs(repmat(timem(i),[length(timem) 1])-timem);
end

ls_sep.dx(1:length(dxi(:)),month) = dxi(:);
ls_sep.dt(1:length(dti(:)),month) = dti(:);
    
end

% output is dx/dt matrices: columns correspond to 12 months; get the cov
% matrix between all pts in a given month, regardless of year, by taking
% that column and forming a square matrix of the non-NaN values.

return
