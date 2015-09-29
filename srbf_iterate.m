%%%% compute mean signal and residuals using (isotropic) space-time covariance

function anom = srbf_iterate(ls_sep,var,covfxn,cp,e,FF,xt_cutoff)

x = sqrt(length(ls_sep.dx));
ind = find(ls_sep.dx <= xt_cutoff(1) & ls_sep.dt <= xt_cutoff(2));
dx = ls_sep.dx(ind); dt = ls_sep.dt(ind);

AA = covfxn(dx,dt,cp);
clear dx dt

[ii,jj] = ind2sub([x x],ind); clear ind

A = sparse(ii,jj,AA,x,x); clear AA

A = A + cp(1)*e*speye(x);

indvar = ~isnan(var);
A = A(indvar,:); A = A(:,indvar);

P = A\FF;
clear A

Z = FF'*P;
B = Z\P';
clear P

bb = FF*B;
clear FF

anom_tmp = var(indvar)-bb*var(indvar);
anom = NaN(length(indvar),1);
anom(indvar) = anom_tmp;

return


