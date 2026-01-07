function [H] = getWeightingMat(rayCountVec,lambda)
%computes the smoothing matrix with given sigma
%nBoxes = length(gridCoord(:,1));

H = diag(exp(-lambda*rayCountVec));
end

