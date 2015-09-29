function [ss_err_anom,ss_err_anls] = ss_err_map(prof,ss_anom,covfxn,covparam,e,grid,ss_sep_A,ss_sep_C,xt_cutoff,F,FF,allmon,srbfparam)

ss_err_anom = NaN(length(grid.lat),length(grid.time));
ss_err_anls = NaN(length(grid.lat),length(grid.time));

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
    
    % C
    name = strcat('time',num2str(i));
    %     C_ind = ismember(ss_sep_C.(name).di(:,1),A_ind);
    
    dx = ss_sep_C.(name).dx;
    dt = ss_sep_C.(name).dt;
    di = ss_sep_C.(name).di;
    clear name
    %     clear C_ind
    
    M = zeros(length(prof.lon),1);
    M(A_ind) = (1:length(A_ind));
    di(:,1) = M(di(:,1));
    clear M
    
    aa = covfxn(dx,dt,covparam);
    clear dx dt
    
    C = sparse(di(:,2),di(:,1),aa,length(grid.lat),length(A_ind));
    C = C(:,indvar);
    
    % anomaly error from objective mapping
    [L,~,P] = chol(A, 'lower') ;
   
    E = zeros (length(grid.lat),1) ;
    parfor k = 1:length(grid.lat)
        Ek = L \ (P' * C (k,:)') ;
        E(k) = Ek' * Ek ;
    end
    
    ei_anom = E;
    clear E
    
    % load FF for all the data points
    m = prof.time; m = m(A_ind);
    [~,m] = datevec(m); mon = unique(m);
    
    FFm = NaN(length(A_ind),srbfparam(1));
    
    for mm = 1:length(mon)
        disp(mm);
        
        namem = strcat('month',num2str(mon(mm)));
        iim = find(m == mon(mm));
        
        indx_ii = A_ind(iim);
        
        iiF = find(allmon == mon(mm));
        FFmm = FF.(namem);
        
        for j = 1:length(iim)
            FFm(iim(j),:) = FFmm(iiF == indx_ii(j),:);
        end
        
    end
    clear FFmm iiF indx_ii iim
    
    FFm = FFm(indvar,:);
    
    % small-scale error due to large-scale error
    P = A\FFm;
    clear A
    
    Z = FFm'*P;
    clear FFm
    
    G = C*P;
    clear C P
    GG = Z\G';
    clear Z
    
    [~,month] = datevec(grid.time(i));
    namem = strcat('month',num2str(month));
    Fm = F.(namem);
    
    ei_anls = sum(G.*GG',2) - 2*sum(Fm.*GG',2);
    
    ss_err_anom(:,i) = ei_anom;
    ss_err_anls(:,i) = ei_anls;
    
end
 
return


