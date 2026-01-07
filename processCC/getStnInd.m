function [z] = getStnInd(stnList,stnName)
%this function was written to get the station index from the list
    z = 0;
for i = 1:1:length(stnList(:,1))
    stnString = convertCharsToStrings(stnList(i,:));

    if(strcmp(stnString,stnName))
        z  = i;
        break;
    end
    

end

