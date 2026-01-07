function [ccOut,ccFlag] = readEndStn(mainPath,subPath,endStnVal,refStnVal)
% this function reads the mat file corresponding to the endstn val and
% tries to extract the CC wrt to the ref stn if present
% endStnVal is a character array;
% refStnVal is a string
fullEndPath = [mainPath,subPath,endStnVal,'.mat'];
ccOut = zeros(1001,1);
ccFlag = 0;
if(exist(fullEndPath))
    A  = load(fullEndPath);
    B = string(A.stnEnd); % convert to string array
    % try to find the refstn Str
    c = find(B==refStnVal,1,'first');
    if(~isempty(c))
        ccOut(:,1) = flipud(A.ccStore(:,c));
        ccFlag = 1;
    else
        ccFlag = 0;
    end
end

end

