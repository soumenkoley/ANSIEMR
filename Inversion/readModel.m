function [h,vp,vs,rho,relErr] = readModel(fName)
% this function reads a .model file output from the code
% extractVelModelGeopsy.m
% outputs a matrix of Vp, Vs, Rho, and depth, corresponding to the
% N best models
% relErr is a column vector corresponding to the relative error of each
% model

fId = fopen(fName);
tLine = fgetl(fId);
% the first line gives information about total number of models
i= 1;
%while(ischar(tLine))
while(~feof(fId))
    tLine = fgetl(fId);
    % find index of "="
    eqInd = strfind(tLine,'=');
    relErr(i,1) = str2double(tLine((eqInd+1):end));
    % number of layers in the model
    tLine = fgetl(fId);
    nLayer = str2double(tLine);
    
    for lNo = 1:1:nLayer
        
        tLine = fgetl(fId);
        sInfo = split(tLine,' ');
        
        h(lNo,i) = str2double(sInfo{1,1});
        vp(lNo,i) = str2double(sInfo{2,1});
        vs(lNo,i) = str2double(sInfo{3,1});
        rho(lNo,i) = str2double(sInfo{4,1});
        %disp('doing!');
    end
    h(:,i) = cumsum(h(:,i));
    i=i+1;
end

fclose(fId);
if(max(h(nLayer-1,:))<800)
    h(nLayer,:) = 1000; % setting the last layer depth to 400 m
else
    h(nLayer,:) = max(h(nLayer-1,:))+100;
end
    
% now plot all the Vs models


end