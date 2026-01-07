function Ls = build_Ls_laplacian(nY, nX, dh)

Ncells = nY * nX;

% preallocate (5 entries per cell in interior)
I = zeros(5*Ncells,1);
J = zeros(5*Ncells,1);
V = zeros(5*Ncells,1);
nnz_count = 0;

idx_from_ij = @(i,j) sub2ind([nY, nX], i, j);

for j = 1:nX
    for i = 1:nY
        idx = idx_from_ij(i,j);

        % center
        nnz_count = nnz_count + 1;
        I(nnz_count) = idx;
        J(nnz_count) = idx;
        V(nnz_count) = -4;

        % up
        if i > 1
            nnz_count = nnz_count + 1;
            I(nnz_count) = idx;
            J(nnz_count) = idx_from_ij(i-1,j);
            V(nnz_count) = 1;
        end

        % down
        if i < nY
            nnz_count = nnz_count + 1;
            I(nnz_count) = idx;
            J(nnz_count) = idx_from_ij(i+1,j);
            V(nnz_count) = 1;
        end

        % left
        if j > 1
            nnz_count = nnz_count + 1;
            I(nnz_count) = idx;
            J(nnz_count) = idx_from_ij(i,j-1);
            V(nnz_count) = 1;
        end

        % right
        if j < nX
            nnz_count = nnz_count + 1;
            I(nnz_count) = idx;
            J(nnz_count) = idx_from_ij(i,j+1);
            V(nnz_count) = 1;
        end
    end
end

I = I(1:nnz_count);
J = J(1:nnz_count);
V = V(1:nnz_count) / (dh^2);   % include 1/h^2 factor

Ls = sparse(I, J, V, Ncells, Ncells);
end