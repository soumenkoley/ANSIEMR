function [L] = build_L_and_residuals_eikonal(rayAtt,dh,nX,nY,...
    step_frac, max_steps, tol_dist)
                   
%% parameters
% rayAtt is the rayAttributeStoreFull variable

%%
stnUnique = unique(rayAtt(:,1));
rC = 1;
for stnNo = 1:1:length(stnUnique)
    % find the entries with a common source stn
    stnInd = find(rayAtt(:,1)==stnUnique(stnNo,1));
    % get the source location
    xs = rayAtt(stnInd(1,1),4);
    ys = rayAtt(stnInd(1,1),5);

    % FMM in TRUE model
    T_true_src = fmm2d_cont_source(sTrue, xs, ys, 0, 0, dh);
    nRays = length(stnInd);
    
    for kk = 1:nRays
        
        xr = rayAttributeStoreFull(stnInd(kk,1),6);
        yr = rayAttributeStoreFull(stnInd(kk,1),7);
        
        ray = backtrack_ray_hybrid(T_true_src,0,0,dh,xs,ys,xr,yr, ...
                                step_frac, max_steps, tol_dist);

        L(rC,:) = ray_to_lengths(ray, 0, 0, dh, nY, nX);
        rC = rC+1;
    end
end
end