function ss_sep = ss_separation_cov(prof,semivargrid)

dinc = semivargrid(3);
tinc = semivargrid(4);

dlag = 0:dinc:semivargrid(1);
tlag = 0:tinc:semivargrid(2);

DX = NaN(3000*length(prof.lon),1);
DT = NaN(3000*length(prof.lon),1);
DDi = NaN(3000*length(prof.lon),2);

latn = prof.lat; i = 1;

for n = 1:length(prof.lon)
    
    disp(n);
    
    lon = prof.lon(n); lat = latn(n); t = prof.time(n);
    
    s = sph_dist_ca(repmat(lon,length(prof.lon),1),prof.lon,repmat(lat,length(prof.lon),1),latn)/1000; % in km
    
    dt = abs(t-prof.time);
    dx = abs(s); clear s
    
    ind = ~isnan(dx) & dx <= max(dlag) & dt <= max(tlag) & ~(dx < 5e-4 & dt == 0);
    dx = dx(ind); dt = dt(ind);
    
    DX(i:i+length(dx)-1) = dx;
    DT(i:i+length(dx)-1) = dt;
    DDi(i:i+length(dx)-1,1) = n;
    DDi(i:i+length(dx)-1,2) = find(ind);
    
    i = i+length(dx);
        
    latn(n) = NaN;
end

DX = DX(1:i-1);
DT = DT(1:i-1);
DDi = DDi(1:i-1,:);

clear i ind lat lon n t dt dx latn

%%%% sort by DX
[DX, ind] = sort(DX); 
DT = DT(ind);
DDi = DDi(ind,:);

%%%% bin by DX
ind_dx = NaN(length(dlag),1);

for k = 1:length(dlag)
    disp(k);
    ind_dx(k) = find(DX <= dlag(k) + dinc/2,1,'last');
end

DX_dx = struct;
DT_dx = struct;
DDi_dx = struct;

for k = 1:length(dlag)
    disp(k);
    name = strcat('dx',num2str(k));
    if k == 1
    DX_dx.(name) = DX(1:ind_dx(k));
    DT_dx.(name) = DT(1:ind_dx(k));
    DDi_dx.(name) = DDi(1:ind_dx(k),:);
    else
    DX_dx.(name) = DX(ind_dx(k-1)-1:ind_dx(k));
    DT_dx.(name) = DT(ind_dx(k-1)-1:ind_dx(k));
    DDi_dx.(name) = DDi(ind_dx(k-1)-1:ind_dx(k),:);
    end
end

clear k name DX DT DDi

%%% sort by DT

for k = 1:length(dlag)
    disp(k);
    name = strcat('dx',num2str(k));
    dt = DT_dx.(name);
    [dt,idt] = sort(dt);
    
    DT_dx.(name) = dt;
    ind_dt.(name) = idt;
    
    dx = DX_dx.(name);
    dx = dx(idt);
    DX_dx.(name) = dx;
    
    dd = DDi_dx.(name);
    dd = dd(idt,:);
    DDi_dx.(name) = dd;
end

clear dd dt dx idt k name ind

%%%% bin by DT

ind_dx_dt = struct;

for k = 1:length(dlag)
    disp(k);
    name = strcat('dx',num2str(k));
    dt = DT_dx.(name);
    idt = NaN(length(tlag),1);
    
    for m = 1:length(tlag)
        idt(m) = find(dt <= tlag(m) + tinc/2,1,'last');
    end
    ind_dx_dt.(name) = idt;
end

clear m idt k name dt ind_dx

DX_dx_dt = struct;
DT_dx_dt = struct;
DDi_dx_dt = struct;

for k = 1:length(dlag)
    disp(k);
    name = strcat('dx',num2str(k));
    
    dx = DX_dx.(name);
    dt = DT_dx.(name);
    dd = DDi_dx.(name);
    ii = ind_dx_dt.(name);
    
    for m = 1:length(tlag)
        name_dt = strcat('dt',num2str(m));
        if m == 1
            DX_dx_dt.(name).(name_dt) = dx(1:ii(m));
            DT_dx_dt.(name).(name_dt) = dt(1:ii(m));
            DDi_dx_dt.(name).(name_dt) = dd(1:ii(m),:);
        else
            DX_dx_dt.(name).(name_dt) = dx(ii(m-1):ii(m));
            DT_dx_dt.(name).(name_dt) = dt(ii(m-1):ii(m));
            DDi_dx_dt.(name).(name_dt) = dd(ii(m-1):ii(m),:);
        end
    end
    
end
clear DX_dx k name dx ii m name_dt DT_dx DDi_dx

%%% remove dx==0/dt==0 pairs (for which dx == 9.5e-5)

dd = DDi_dx_dt.dx1.dt1;
ii = find(~(dd(:,1) == dd(:,2)));
dd = dd(ii,:);
DDi_dx_dt.dx1.dt1 = dd;

dd = DX_dx_dt.dx1.dt1;
dd = dd(ii);
DX_dx_dt.dx1.dt1 = dd;

dd = DT_dx_dt.dx1.dt1;
dd = dd(ii);
DT_dx_dt.dx1.dt1 = dd;

ss_sep = struct;
ss_sep.dx = DX_dx_dt;
ss_sep.dt = DT_dx_dt;
ss_sep.di = DDi_dx_dt;

return


