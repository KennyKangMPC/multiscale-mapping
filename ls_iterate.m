function [covparam,file_val] = ls_iterate(prof,varz,ls_sep,ss_sep,tolerance,covfxn,semivarfxn,semivarfxn_cf,semivargrid,covparam0,FF,e,xt_cutoff)

file_val = 1;
stop_val = 0;

cp = covparam0;
options = optimoptions('lsqnonlin','TolFun',1e-16,'MaxFunEvals',10000,'MaxIter',5000,'TolX',1e-16); %,'Algorithm','levenberg-marquardt');

timevec = datevec(prof.time); timevec = timevec(:,2);

while stop_val == 0
    anom = NaN(length(varz),1);
    
    for month = 1:12
        
        ind = timevec == month;
        varm = varz(ind);
        
        indls = ~isnan(ls_sep.dx(:,month));
        ls_sepm = struct('dx',ls_sep.dx(indls,month),'dt',ls_sep.dt(indls,month));
        clear indls
        
        FFm = FF.(strcat('month',num2str(month)));
        
        ai = srbf_iterate(ls_sepm,varm,covfxn,cp,e,FFm,xt_cutoff);
        
        anom(ind) = ai;
        
    end
    clear ai month ind ls_sepm FFm varm
    
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
    
    [cp_new,~,~,~,~] = lsqnonlin(semivarfxn,cp,[],[],options,DX);
    clear DX
    
    if  max(abs(semivarfxn_cf(cp_new,dx(:),dt(:))-semivarfxn_cf(cp,dx(:),dt(:)))) < tolerance*max([(semivarfxn_cf(cp_new,dx(:),dt(:))); semivarfxn_cf(cp,dx(:),dt(:))])
        stop_val = 1;
    end
    disp(file_val)
    disp(cp_new)
    
    file_val = file_val+1;
    
    cp = cp_new;
    
end

covparam = cp;

return


