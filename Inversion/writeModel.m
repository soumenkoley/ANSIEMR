function [] = writeModel(h,vp,vs,rho,savePath)
%this function creates a model file to input into gpdc

fId = fopen(savePath,'w');

nLayer = length(h(:,1));
nModel = length(h(1,:));

%hThick = zeros(nLayers,1);


for i = 1:1:nModel
    hThick = [0;h(:,i)];
    hThick(end,1) = hThick(end-1,1);
    hThickNow = hThick(2:end,1)-hThick(1:(end-1),1);
    fprintf(fId,'%d\n',nLayer);
    for j = 1:1:nLayer
        fprintf(fId,'%4.4f\t',hThickNow(j,1));
        fprintf(fId,'%4.4f\t',vp(j,i));
        fprintf(fId,'%4.4f\t',vs(j,i));
        fprintf(fId,'%4.4f\n',rho(j,i));
    end
end
fclose(fId);
end

