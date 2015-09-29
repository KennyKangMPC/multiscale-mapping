function [ls_var,ss_anom,ls_err] = ls_map(prof,varz,ls_sep,covfxn,cp,FF,F,e,Ngrid,xt_cutoff)

timevec = datevec(prof.time); timevec = timevec(:,2);
ss_anom = NaN(length(varz),1);
ls_var = NaN(Ngrid,12);
ls_err = NaN(Ngrid,12);

for month = 1:12
    
    ind = timevec == month;
    varm = varz(ind);
    
    indls = ~isnan(ls_sep.dx(:,month));
    ls_sepm = struct('dx',ls_sep.dx(indls,month),'dt',ls_sep.dt(indls,month));
    clear indls
    
    FFm = FF.(strcat('month',num2str(month)));
    Fm = F.(strcat('month',num2str(month)));
    
    [vi,ai,ei] = srbf_map(ls_sepm,varm,covfxn,cp,e,FFm,Fm,xt_cutoff);

    ss_anom(ind) = ai;
    ls_var(:,month) = vi;
    ls_err(:,month) = ei;
    
end

return