function ss_var = ss_map(prof,ss_anom,covfxn,covparam,e,grid,ss_sep_A,ss_sep_C,xt_cutoff)

ss_var = NaN(length(grid.lat),length(grid.time));

for i = 1:length(grid.time)
    
    disp(i);
    
    % A
    A_ind = find(abs(prof.time-grid.time(i)) <= xt_cutoff(2));
    jj = ismember(ss_sep_A.di(:,1),A_ind) & ismember(ss_sep_A.di(:,2),A_ind);
    
    dx = ss_sep_A.dx(jj);
    dt = ss_sep_A.dt(jj);
    di = ss_sep_A.di(jj,:);
    
    M = zeros(length(prof.time),1);
    M(A_ind) = (1:length(A_ind));
    di = M(di);
    
    aa = covfxn(dx,dt,covparam);
    clear dx dt 
    
    A = sparse(di(:,1),di(:,2),aa,length(A_ind),length(A_ind));
    A = triu(A,1) + sparse(di(:,2),di(:,1),aa,length(A_ind),length(A_ind));
    A = A + e*covparam(1)*speye(length(A_ind));
    clear di aa
    
    % b
    ss = ss_anom(A_ind);
    indvar = ~isnan(ss);
    
    A = A(indvar,:); A = A(:,indvar);
    
    b = A\ss(indvar);
    clear A ss
    
    % C
    name = strcat('time',num2str(i));
    
    dx = ss_sep_C.(name).dx;
    dt = ss_sep_C.(name).dt;
    di = ss_sep_C.(name).di;
    
    M = zeros(length(prof.lon),1);
    M(A_ind) = (1:length(A_ind));
    di(:,1) = M(di(:,1));
    clear M
    
    aa = covfxn(dx,dt,covparam);
    clear dx dt 
    
    C = sparse(di(:,2),di(:,1),aa,length(grid.lat),length(A_ind));
    C = C(:,indvar);
    
    ss_i = C*b;
    
    ss_var(:,i) = ss_i;
 
end


return
