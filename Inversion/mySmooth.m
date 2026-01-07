function [outVec] = mySmooth(inpVec,sSpan)
% similar to matlab smooth, but takes care of NaN entries
%firstNonNaNInd = find(~isnan(inpVec), 1);
%lastNonNaNInd = find(~isnan(inpVec), 1,'last');

%nanVec = firstNonNaNInd:lastNonNaNInd;
nanVec = find(~isnan(inpVec));
% there could be NaN in between as well

nanC = 1;
if(~isempty(nanVec))
    nanLimits(1,1) = nanVec(1,1);
    if(length(nanVec)>1)
        for i = 1:1:(length(nanVec)-1)
            diffNow = nanVec(i+1)-nanVec(i);
            if(diffNow>1)
                nanLimits(nanC,2) = nanVec(i);
                nanC = nanC+1;
                nanLimits(nanC,1) = nanVec(i+1);
            else
                nanLimits(nanC,2) = nanVec(i+1);
            end
        end
    else
        nanLimits(1,1) = nanVec(1,1);
        nanLimits(1,2) = nanVec(1,1);
    end
        
end
outVec = inpVec;

if(~isempty(nanVec))
    for i = 1:1:length(nanLimits(:,1))
        inpNew = inpVec(nanLimits(i,1):nanLimits(i,2));
        if(length(inpNew)<sSpan)
            inpNew = smooth(inpNew);
        else
            inpNew = smooth(inpNew,sSpan);
        end
        
        outVec(nanLimits(i,1):nanLimits(i,2)) = inpNew;
    end
end
end

