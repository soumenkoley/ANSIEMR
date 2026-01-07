function [L, dT] = build_L_and_residuals_eikonalN( ...
        rayAtt, s_model, Tobs, dh, nX, nY, ...
        step_frac, max_steps, tol_dist)

% rayAtt: [Npaths x 8] rayAttributeStoreFull (as in your script)
% s_model: current slowness model [nY x nX]
% Tobs: synthetic observed travel times [Npaths x 1]
% returns:
%   L   : Npaths x (nX*nY) sparse sensitivity matrix
%   dT  : Npaths x 1 data residual vector (Tobs - Tpred)

Ncells = nX * nY;
Npaths = numel(Tobs);

% Preallocate sparse triplets (rough guess: 50 nonzeros per path)
max_nnz = Npaths * 50;
I = zeros(max_nnz,1);
J = zeros(max_nnz,1);
V = zeros(max_nnz,1);
nnz_count = 0;

dT  = zeros(Npaths,1);
row = 0;

x0 = 0; y0 = 0;  % you’re using 0,0 as grid origin

stnUnique = unique(rayAtt(:,1));

for sIdx = 1:numel(stnUnique)

    stnID = stnUnique(sIdx);

    % indices in rayAtt where this station is source
    stnInd = find(rayAtt(:,1) == stnID);
    nRays  = numel(stnInd);

    % source coordinates (already updated in your script)
    xs = rayAtt(stnInd(1),4);
    ys = rayAtt(stnInd(1),5);

    % FMM in CURRENT model
    T_src = fmm2d_cont_source(s_model, xs, ys, x0, y0, dh);

    for kk = 1:nRays
        row = row + 1;

        xr = rayAtt(stnInd(kk),6);
        yr = rayAtt(stnInd(kk),7);

        % predicted travel time in current model
        Tpred = interp_bilinear(T_src, xr, yr, x0, y0, dh);

        % residual
        dT(row) = Tobs(row) - Tpred;

        % backtrack ray
        ray = backtrack_ray_hybrid(T_src, x0, y0, dh, ...
                                   xs, ys, xr, yr, ...
                                   step_frac, max_steps, tol_dist);

        % path lengths per cell
        Lrow = ray_to_lengths(ray, x0, y0, dh, nY, nX);  % 1 x Ncells

        nz = find(Lrow ~= 0);
        nn = numel(nz);

        if nnz_count + nn > max_nnz
            % grow storage
            I = [I; zeros(max_nnz,1)];
            J = [J; zeros(max_nnz,1)];
            V = [V; zeros(max_nnz,1)];
            max_nnz = max_nnz * 2;
        end

        I(nnz_count+1:nnz_count+nn) = row;
        J(nnz_count+1:nnz_count+nn) = nz;
        V(nnz_count+1:nnz_count+nn) = Lrow(nz);

        nnz_count = nnz_count + nn;
    end
end

% trim and build sparse L
I = I(1:nnz_count);
J = J(1:nnz_count);
V = V(1:nnz_count);

L = sparse(I, J, V, row, Ncells);
end