%%% example script for multi-scale mapping
% this example maps profiles, i.e. data from multiple depths at scattered
% data locations (lat,lon,time)
% for each depth, a large-scale monthly climatology is computed using all
% data in that month, regardless of year.  then a small-scale field is
% calculated for each individual month within the time range of interest


%% load data

load data.mat 
% data file should have five variables - lat, lon, time, var (profiles of variable to be mapped), and z (depth/pressure)

% make a structure with data locations
prof = struct;
prof.lat = lat; prof.lon = lon; prof.time = time;
clear time lat lon;

% make longitude in 0-360 range 
prof.lon(prof.lon > 360) = prof.lon(prof.lon > 360) - 360;
prof.lon(prof.lon < 0) = prof.lon(prof.lon < 0) + 360;

%% calculate initial fields and parameters 

%%%% specify space-time grid for covariance/semivariogram estimation [max-dx, dx-inc; max-dt, dt-inc];
% this gives the maximum separations over which the covariance is
% estimated and the bin size for both space and time
semivargrid = [3000 50; 1200 10];

%%%% calculate space-time separations used in covariance estimation
ss_sep = ss_separation_cov(prof,semivargrid);

%%%% calculate space-time separtions used in large-scale field mapping (within each month)
ls_sep = ls_separation(prof);

%%%% assumed structures of semivariogram and covariance function
% this covariance function is gaussian in space, exponential in time
% cp (covariance parameter) = [variance x-length scale t-length scale]
% x = spatial separation (isotropic), t = absolute value of temporal separation
covfxn = @(x,t,cp)(cp(1)*exp(-(x/cp(2)).^2 - (t/cp(3))));  
% this semivariogram corresponds to the above covariance function and is used in fitting to the estimated semivariogram; 
% semivarfxn = weights*((predicted variance-predicted covariance) - estimated semivariogram)
% DX = [spatial-separation  temporal-separation  weights  estimated-semivariogram] 
semivarfxn = @(cp,DX)(DX(:,3).*(cp(1)*(1 - exp(-DX(:,1).^2/cp(2)^2 - DX(:,2)/cp(3)))-DX(:,4)));
% same semivariogram but this form is used in comparing two semivariograms from different iterations to determine if convergence has been reached
semivarfxn_cf = @(cp,dx,dt)(L(1)*(1 - exp(-dx.^2/L(2)^2 - dt/L(3))));

% initial guess for covariance parameters, here in units of [(var units)^2 km days]
cp0 = [1 100 20];

% variance in sub-grid scales (as a fraction of variance of small-scale field)
% here we use 50% as a conservative value
e = .5;

% tolerance level for covariance convergence
tolerance = 1e-3;

% load grid and mask
load mask.mat
% mask file should have three variables corresponding to mapped grid:
% long-vector of longitudes, latg-vector of latitudes, mask-2D logical matrix with 0 at locations of land and 1 at ocean locations
[maplon,maplat] = meshgrid(long,latg);
maplat = maplat(mask);
maplon = maplon(mask);
Ngrid = length(maplat(:));

% cutoff for space-time covariance [maxds maxdt]
% for computational efficiency, covariance is set to 0 if the spatial separation is greater than xt_cutoff(1) OR the temporal separation is greater than xt_cutoff(2)
% adjusting this will greatly affect computational time - the larger the cutoff, the more time to do the mapping
xt_cutoff = [3000 120]; % units are same as for cp (here, [km days])

% times of small-scale estimates in matlab datenum format (here, midpoint of individual months from Jan 2011 to Dec 2014 - change to be appropriate for your data set)
ss_time = [2011*ones(12,1); 2012*ones(12,1); 2013*ones(12,1); 2014*ones(12,1)];
ss_time = [ss_time repmat((1:12)',4,1) 15*ones(length(ss_time),1)];
ss_time = datenum(ss_time);

% make structure with grid variables
grid = struct;
grid.lat = maplat(:);
grid.lon = maplon(:);
grid.time = ss_time(:);
clear ss_time

% calculate small scale separations used in forming covariance matrices for small-scale mapping
% (A=data-data covariance, C=grid-data covariance)
[ss_sep_A,ss_sep_C] = ss_separation(prof,grid,xt_cutoff);

% specify measurement uncertainty of variables, in same units as var -- change to be appropriate for your variable
err = 2;

% all of these variables are saved here in case mapping needs to be restarted at any point
save initial_setup.mat -v7.3


%% mapping

%%%% if restarting, run these lines
% clear; load initial_setup.mat

% for each depth, perform mapping
for d = 1:length(z)
    
    % select variable at that depth
    varz = squeeze(var(:,z == z(d)));
    
    % add random error to each measurement
    errz = -err + 2*err*rand(size(varz));
    varz = varz+errz;
    
    % initial large-scale calculation - assume small-scale is random and uncorrelated, calculate srbf parameters
    [covparam0, srbfparam] = ls_0(prof,varz,ss_sep,semivarfxn,semivargrid,cp0);
    
    % calculate srbf's at data locations using compute srbf parameters
    FF = struct;
    timevec = datevec(prof.time); timevec = timevec(:,2);
    for month = 1:12
        disp(month)
        ind = timevec == month & ~isnan(varz);
        
        lon = prof.lon(ind);
        lat = prof.lat(ind);
        
        ff = constructrbf_FF(lon,lat,srbfparam);
        
        FF.(strcat('month',num2str(month))) = ff;
        
    end
    clear lon lat ind month ff timevec
    
    % iterate large-scale field-covariance calculation to determine covariance parameters
    covparam = ls_iterate(prof,varz,ls_sep,ss_sep,tolerance,covfxn,semivarfxn,semivarfxn_cf,semivargrid,covparam0,FF,e,xt_cutoff);
    
    % compute srbf's at grid points
    F = struct;
    timevec = datevec(prof.time); timevec = timevec(:,2);
    for month = 1:12
        disp(month)
        ind = timevec == month & ~isnan(varz);
        
        lon = prof.lon(ind);
        lat = prof.lat(ind);
        
        f = constructrbf_F(lon,lat,srbfparam,maplon(:),maplat(:));
        
        F.(strcat('month',num2str(month))) = f;
        
    end
    clear lon lat ind month f timevec
    
    % compute final large-scale field, small-scale anomalies, error on large-scale field
    [ls_var,ss_anom,ls_err] = ls_map(prof,varz,ls_sep,covfxn,covparam,FF,F,e,Ngrid,xt_cutoff);
    
    % compute small scale anomalies 
    ss_var = ss_map(prof,ss_anom,covfxn,covparam,e,grid,ss_sep_A,ss_sep_C,xt_cutoff);
    
    % compute small scale errors -- this step can take a long time so comment out these lines for quicker run time and run separately using parallel toolbox (shown below)
    [~,allmon] = datevec(prof.time); % vector of months of all data points
    [ss_err_anom,ss_err_anls] = ss_err_map(prof,varz,covfxn,covparam,e,grid,ss_sep_A,ss_sep_C,xt_cutoff,F,FF,allmon,srbfparam);
    
    % save mapped variables and covariance function (if not calculating small-scale errors here, remove last two variables from this line)
    save(strcat('mapped_var_',num2str(d),'.mat'),'ss_var','covparam','covfxn','e','xt_cutoff','ls_var','ls_err','srbfparam','ss_err_anom','ss_err_anls');
end

%% small-scale errors (optional, if not done above)

clear; load initial_setup.mat

% start parallel computing with 4 processors (adjust number for your machine)
parpool(4);

for d = 1:length(z)
    
    % select variable at that depth
    varz = squeeze(var(:,z == z(d)));
    
    % load mapped variables
    load(strcat('mapped_var_',num2str(d),'.mat'));
    
    % calculate srbf's at data locations
    FF = struct;
    timevec = datevec(prof.time); timevec = timevec(:,2);
    for month = 1:12
        disp(month)
        ind = timevec == month & ~isnan(varz);
        
        lon = prof.lon(ind);
        lat = prof.lat(ind);
        
        ff = constructrbf_FF(lon,lat,srbfparam);
        
        FF.(strcat('month',num2str(month))) = ff;
        
    end
    clear lon lat ind month ff timevec
    
    % compute srbf's at grid points
    F = struct;
    timevec = datevec(prof.time); timevec = timevec(:,2);
    for month = 1:12
        disp(month)
        ind = timevec == month & ~isnan(varz);
        
        lon = prof.lon(ind);
        lat = prof.lat(ind);
        
        f = constructrbf_F(lon,lat,srbfparam,maplon(:),maplat(:));
        
        F.(strcat('month',num2str(month))) = f;
        
    end
    clear lon lat ind month f timevec

    % vector of months of all data points 
    [~,allmon] = datevec(prof.time); 

    % compute small scale errors 
    [ss_err_anom,ss_err_anls] = ss_err_map(prof,varz,covfxn,covparam,e,grid,ss_sep_A,ss_sep_C,xt_cutoff,F,FF,allmon,srbfparam);
       
    % save mapped variables and covariance function
    save(strcat('mapped_var_sserr_',num2str(d),'.mat'),'ss_err_anom','ss_err_anls');
end

    
