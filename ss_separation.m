function [ss_sep_A,ss_sep_C] = ss_separation(prof,grid,xt_cutoff)

% prof = structure with lat, lon, time of data
% grid = structure with lat, lon, time of grid
% xt_cutoff = [maxds maxdt]

%%%% A
ss_sep_A = struct;

ss_sep_A.dx = NaN(300*length(prof.lat),1);
ss_sep_A.dt = NaN(300*length(prof.lat),1);
ss_sep_A.di = NaN(300*length(prof.lat),2);

latn = prof.lat; i = 1;

for n = 1:length(prof.lat)
    disp(n);
    
    lon = prof.lon(n); lat = latn(n); t = prof.time(n);
    
    s = sph_dist_ca(repmat(lon,length(prof.lon),1),prof.lon,repmat(lat,length(prof.lon),1),latn)/1000; % in km
    
    dt = abs(t-prof.time);
    dx = abs(s); clear s
    
    ind = ~isnan(dx) & dx <= xt_cutoff(1) & dt <= 2*xt_cutoff(2);
    dx = dx(ind); dt = dt(ind);
    
    ss_sep_A.dx(i:i+length(dx)-1) = dx;
    ss_sep_A.dt(i:i+length(dx)-1) = dt;
    ss_sep_A.di(i:i+length(dx)-1,1) = n;
    ss_sep_A.di(i:i+length(dx)-1,2) = find(ind);
    
    i = i+length(dx);
    
    latn(n) = NaN;
end

ss_sep_A.dx = ss_sep_A.dx(1:i-1);
ss_sep_A.dt = ss_sep_A.dt(1:i-1);
ss_sep_A.di = ss_sep_A.di(1:i-1,:);

%%%% C
ss_sep_C = struct;

for t = 1:length(grid.time)
    
    disp(t);
    
    ind = find(abs(prof.time - grid.time(t)) <= xt_cutoff(2));
    
    name = strcat('time',num2str(t));
    ss_sep_C.(name).dx = NaN(300*length(ind),1);
    ss_sep_C.(name).dt = NaN(300*length(ind),1);
    ss_sep_C.(name).di = NaN(300*length(ind),2);
    ss_sep_C.(name).ind = ind;
    
    i = 1;
    
    data = [prof.lat prof.lon prof.time (1:length(prof.lon))'];
    data = data(ind,:);
    
    for n = 1:length(data)
       
        lonn = data(n,2); latn = data(n,1); tn = data(n,3); nn = data(n,4);
        
        s = sph_dist_ca(repmat(lonn,length(grid.lon),1),grid.lon,repmat(latn,length(grid.lon),1),grid.lat)/1000; % in km
        
        dx = abs(s); clear s
        
        ind = ~isnan(dx) & dx <= xt_cutoff(1);
        
        ss_sep_C.(name).dx(i:i+sum(ind)-1) = dx(ind);
        ss_sep_C.(name).dt(i:i+sum(ind)-1) = abs(tn - grid.time(t));
        ss_sep_C.(name).di(i:i+sum(ind)-1,1) = nn;
        ss_sep_C.(name).di(i:i+sum(ind)-1,2) = find(ind);
        
        i = i+sum(ind);
        
    end
    
    ss_sep_C.(name).dx = ss_sep_C.(name).dx(1:i-1);
    ss_sep_C.(name).dt = ss_sep_C.(name).dt(1:i-1);
    ss_sep_C.(name).di = ss_sep_C.(name).di(1:i-1,:);
    
end

return


