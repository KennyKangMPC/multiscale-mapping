function [covparam0, srbfparam] = ls_0(prof,var,ss_sep,semivarfxn,semivargrid,cp0)
options = optimoptions('lsqnonlin','TolFun',1e-16,'MaxFunEvals',10000,'MaxIter',5000,'TolX',1e-16); %,'Algorithm','levenberg-marquardt');

indvar = find(~isnan(var));

maxN = round(length(indvar)/20/12,2,'significant');

timevec = datevec(prof.time); timevec = timevec(:,2);

anom = NaN(length(prof.lon),1);
G = NaN(12,1);
L = NaN(12,1);

for month = 1:12
    disp(month)
    ind = find(timevec == month & ~isnan(var));
    
    profm = struct('lat',prof.lat(ind),'lon',prof.lon(ind),'time',prof.time(ind));
    
    [ai,l,g] = srbf_0_L(profm,var(ind),maxN);
    
    anom(ind) = ai;
    G(month) = g;
    L(month) = l;
end

clear ai g l profm ind month

dinc = semivargrid(3);
tinc = semivargrid(4);
dlag = 0:dinc:semivargrid(1);
tlag = 0:tinc:semivargrid(2);

semivar=zeros(length(dlag),length(tlag));
num=zeros(length(dlag),length(tlag));
distx=zeros(length(dlag),length(tlag));
distt=zeros(length(dlag),length(tlag));

for k = 1:length(dlag)
    for j = 1:length(tlag)
%         disp([k j]);
        namek = strcat('dx',num2str(k));
        namej = strcat('dt',num2str(j));
        dx = ss_sep.dx.(namek).(namej);
        dt = ss_sep.dt.(namek).(namej);
        
        ddi = ss_sep.di.(namek).(namej);
        dd1 = anom(ddi(:,1));
        dd2 = anom(ddi(:,2));
        
        dd = (dd1 - dd2).^2;
        ind = find(~isnan(dd));
        
        num(k,j) = length(dx(ind));
        distx(k,j) = sum(dx(ind));
        distt(k,j) = sum(dt(ind));
        semivar(k,j) = sum(dd(ind));
        
    end
end

clear DDi_dx_dt DX_dx_dt DT_dx_dt ddi dd1 dd2 dd name* j k

variogram = semivar./num./2;
dx = distx./num;
dt = distt./num;
clear distx distt

weights = sqrt((num./2./variogram.^2));
DX = [dx(:) dt(:) weights(:) variogram(:)];

[covparam0,~,~,~,~] = lsqnonlin(semivarfxn,cp0,[],[],options,DX);
clear DX

srbfparam = [maxN; mean(G); mean(L)];


return