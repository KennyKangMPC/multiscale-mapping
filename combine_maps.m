%% combine mapped fields
clear;

load initial_setup.mat z latg long mask grid

% initialize matrix for all the mapped variables
mapped_var = NaN(length(latg),length(long),length(z),length(grid.time));

% number of years in the data set
num_yr = 4;

for d = 1:length(z)
    load(strcat('mapped_var_',num2str(d),'.mat'),'ls_var','ss_var');
    
    o2 = ss_var + repmat(ls_var,1,num_yr);
    tmp = NaN(length(latg),length(long));
    
    for t = 1:num_yr*12
        tmp(mask) = o2(:,t);
        mapped_var(:,:,d,t) = tmp;
    end
    
end

% compute mean field
mean_var = mean(mapped_var,4);

time = [grid.time];
lat = latg;

save mapped_var_all.mat lat long mapped_var mean_var z time

%% combine mapping uncertainties
clear;

load initial_setup.mat z latg long mask grid

% number of years in the data set
num_yr = 4;

mapped_var_sqerr = NaN(length(latg),length(long),length(z),length(grid.time));
mean_var_sqerr = NaN(length(latg),length(long),length(z));

for d = 1:length(z)
    load(strcat('mapped_var_err_',num2str(d),'.mat'),'ss_err_anom','ss_err_anls');
    load(strcat('mapped_var_',num2str(d),'.mat'),'ls_err','covparam');
    
    err_var = covparam(1) - ss_err_anom + repmat(ls_err,1,4) + ss_err_anls;
    tmp = NaN(length(latg),length(long));
    
    for t = 1:num_yr*12
        tmp(mask) = err_var(:,t);
        mapped_var_sqerr(:,:,d,t) = tmp;
    end
    
    tmp(mask)= sum(ls_err,2)./12^2 + sum((covparam(1) - ss_err_anom + ss_err_anls),2)./(length(grid.time))^2;
    mean_var_sqerr(:,:,d) = tmp;
end


time = [grid.time];
lat = latg;

save mapped_var_error_all.mat lat long mapped_var_sqerr mean_var_sqerr z time

