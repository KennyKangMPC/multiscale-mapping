%%%% compute mean signal and residuals for test field (cartesian
%%%% coordinates)

function [anom,L,g] = srbf_0_L(prof,var,maxN)

% construct radial basis functions at grid points (F) and data points (FF)

[FF,L,g] = constructrbf_L(prof.lon,prof.lat,maxN);

% compute coefficients anomalies at data points

Z = FF'*FF;
B = Z\FF';

bb = FF*B;
clear B FF

anom = var-bb*var;

return


