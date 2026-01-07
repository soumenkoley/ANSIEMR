% this script performs the checkerboard test using eikonal tomography
clear; close all;

%%
% load all cartesian node locations
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

dh = 200; % spacing along x and y coordinates
fGet = 1.8; % frequency to extract group or phase velocities
fExt = fGet;
fGetInd = find(fExt >= fGet,1,'first');
vBG = 800; % background velocity of the model in m/s
vAmp = 0.3; % percentage variation in the background for checkerboard
L = 1000; % length of the perturbation in m
sigmaFrac = 0.1;   % 3% noise to the theoretical travel times
step_frac = 0.5;
max_steps = 2000;
tol_dist  = 0.5*dh;
nIter = 4;  % number of nonlinear iterations for Eikonal tomography
%%
% first load the travel times along with the rayAttribute files
% loads the variables rayAttributeStoreFull, timePickStoreFull
% rayAttributeStoreFull = (8 x nRays), col 1, stnId (integer)
% col2 = stnId2 (integer), col3 = distance between stations, col4-5 = (x,y) loc stn1, col(5-6) = (x,y) stn2
% col8 = azimuth between stations
% all locations are in meters
load(['C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\LowVelAug\SNR7\grpTimeNew',num2str(fGet),'.mat']);
% where to store, the ouput maps
fPathStore= 'B:\LimburgBigSurvey1CC-Pair\VelMapsApril\';
% specify the name of the file
fStoreName = ['grpVel',num2str(fGet),'.mat'];

%% first relocate the cartesian node locations on a grid so that
% all node coordinates are positive
% since I do this, I have to  update the stn1 and stn2 locations
% in ratAttributeStoreFull as well, done after this
[newLoc,xLimit,yLimit] = setUpGeometry(dh,dh);
ifPlot = 1;
if(ifPlot)
    figure(1)
    subplot(1,2,1)
    plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
    hold on;
end

% get transpose to make it (nRays x 8)
rayAttributeStoreFull = rayAttributeStoreFull';
[s1,s2] = size(rayAttributeStoreFull);
nX = floor(xLimit/dh)+1; nY = floor(yLimit/dh)+1;
nBoxes = nX*nY;

% plot the rays in figure 1
for rayNo = 1:1:s1
    if(~isnan(phaseVelStoreFull(fGetInd,rayNo)))
        nodeA = rayAttributeStoreFull(rayNo,1); nodeB = rayAttributeStoreFull(rayNo,2);
        nodeAInd = find(newLoc(:,1)==nodeA); nodeBInd = find(newLoc(:,1)==nodeB);
        nodeALoc = newLoc(nodeAInd,2:3); nodeBLoc = newLoc(nodeBInd,2:3);
        rayStart = [nodeALoc(1,1),nodeBLoc(1,1)]';
        rayEnd = [nodeALoc(1,2),nodeBLoc(1,2)]';
        % update the station locations based on new geometry
        rayAttributeStoreFull(rayNo,4:5) = nodeALoc;
        rayAttributeStoreFull(rayNo,6:7) = nodeBLoc;
        plot(rayStart,rayEnd,'k','LineWidth',1);
    end
end
hold off;

%% prepare the checkerboard
[sTrue,x,y] = getCheckBoard(vBG,L,vAmp,nX,nY,dh);
% plot the checker board
figure(2); hold on;
imagesc(x,y,1./sTrue);
plot(newLoc(:,2),newLoc(:,3),'ko','MarkerSize',6,'MarkerFaceColor','k');
hold off;

%% construct the theoretical travel times
% so determine the number of unique entries in the source station
% it is vasically a repeated scenario of 1-2, 1-3, etc
stnUnique = unique(rayAttributeStoreFull(:,1));
rC = 1;
for stnNo = 1:1:length(stnUnique)
    % find the entries with a common source stn
    stnInd = find(rayAttributeStoreFull(:,1)==stnUnique(stnNo,1));
    % get the source location
    xs = rayAttributeStoreFull(stnInd(1,1),4);
    ys = rayAttributeStoreFull(stnInd(1,1),5);

    % FMM in TRUE model
    T_true_src = fmm2d_cont_source(sTrue, xs, ys, 0, 0, dh);
    nRays = length(stnInd);
    
    for kk = 1:nRays
        
        xr = rayAttributeStoreFull(stnInd(kk,1),6);
        yr = rayAttributeStoreFull(stnInd(kk,1),7);

        % synthetic travel time = true travel time at receiver
        Tobs(rC,1) = interp_bilinear(T_true_src, xr, yr, 0, 0, dh);
        rC = rC+1;
    end
end
%% add some noise to the travl times
noise = sigmaFrac * Tobs .* randn(size(Tobs));
Tobs_noisy = Tobs + noise;
Npaths = length(Tobs_noisy);
sigma = sigmaFrac * Tobs_noisy;   % or Tobs, they’re similar
Wd = spdiags(1./sigma, 0, Npaths, Npaths);
%% perform the iterative inversion
sBG = 1/vBG;
s_inv = sBG*ones(nY,nX); % start with an uniform model
Ls = build_Ls_laplacian(nY, nX, dh);

for it = 1:nIter

    fprintf('Nonlinear iter %d\n', it);

    % Build L and data residuals for CURRENT model
    [L, dT] = build_L_and_residuals_eikonalN( ...
        rayAttributeStoreFull, s_inv, Tobs_noisy, dh, nX, nY, ...
        step_frac, max_steps, tol_dist);
    % Regularization matrices (you can reuse your Barmin-style)
    % Ls: smoothing operator (e.g., discrete Laplacian)
    % Wd: data weights (e.g., diag(1/sigma))
    % Setup normal equations:

    Ncells = nY * nX;
    %Wd = speye(Npaths);     % or 1/sigma weighting
    
    Ad = (Wd*L).'*(Wd*L);
    As = Ls.'*Ls;
    
    scale_d = median(diag(Ad));
    scale_s = median(diag(As));
    
    %target_ratio = 0.000004;  % 20% of data term
    %lambda_s = target_ratio * (scale_d / scale_s);
    %lambda_m = 100*scale_d;
    raylen  = full(sum(L, 1)).';
    cov     = raylen / max(raylen);

    epsilon = 0.05;
    Wcov = 1 ./ (cov + epsilon);
    Wcov_mat = spdiags(Wcov, 0, Ncells, Ncells);
    
    lambda_m = 1e-3*scale_d;              % small ridge
    alpha    = 1000;                        % smoothing strength
    lambda_s = alpha * (scale_d / scale_s); % gives ~5% of data term
    lambda_cov = 0.1 * scale_d; 
    
    %lambda_m = 1e-3 * scale_d;
    %lambda_s = 0.9 * scale_d;   % scan alpha = 0..0.2
    % Ls should be N_smooth x Ncells (we can define this separately)
    % Example: Ls is 2D Laplacian on the grid.

    %A = Ad + lambda_s * As + lambda_m * speye(Ncells);
    A = Ad + lambda_s*As + lambda_m*speye(Ncells) + lambda_cov * Wcov_mat;
    b = (Wd*L).' * (Wd * dT);

    % Solve for slowness perturbation
    delta_s = A \ b;
    
    s_old   = s_inv(:);
    max_rel = max(abs(delta_s)) / mean(abs(s_old));
    if max_rel > 0.2
        delta_s = 0.2/max_rel*delta_s;
    end
    
    s_vec = s_inv(:) + delta_s;
    s_inv = reshape(s_vec, nY, nX);
    errT(it,1) = sqrt(mean(dT.^2));
    rms_misfit = sqrt(mean(dT.^2));
    fprintf('Iter %d: RMS misfit = %.4f s\n', it, rms_misfit);
end
V_rec = 1 ./ s_inv(:);
cov_thresh = 0.01;
mask = cov >= cov_thresh;
V_rec(~mask) = NaN;
V_mat = reshape(V_rec, nY, nX);

figure(4);
imagesc(x,y,V_mat);
set(gca,'YDir','normal');
%axis equal tight;
colorbar;
%caxis([1000 2000]);
title('Recovered velocity (masked, coverage-weighted)');