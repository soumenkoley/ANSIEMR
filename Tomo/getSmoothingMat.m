function [F] = getSmoothingMat(gridCoord,smoothFactor)
%computes the smoothing matrix with given sigma
nBoxes = length(gridCoord(:,1));
F = zeros(nBoxes,nBoxes);
for i = 1:1:nBoxes
    sumFij = 0;
    for j = 1:1:nBoxes
        if(i~=j)
            rij = sum((gridCoord(i,1:2)-gridCoord(j,1:2)).^2);
            F(i,j) = exp(-rij/2/(smoothFactor^2));
            sumFij = sumFij+F(i,j);
        end
    end
    F(i,:) = -F(i,:)/sumFij;
    F(i,i) = 1;
end
end