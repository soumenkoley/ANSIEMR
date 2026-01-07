function [fNow,dispAll] = readDispModel(dispPath)
%this function reads the group or phase velocity dispersion files
% obtained using gpdc
fId = fopen(dispPath);

tLine = fgetl(fId);
oneDisp = 0;
fC = 1;
dC = 1;
while(~feof(fId))
    tLine = fgetl(fId);
    if(tLine(1,1)=='#')
        newRead = 0;
        if(oneDisp)
            dispNowIntp = interp1(dispNow(:,1),dispNow(:,2),[1.0:0.05:5]','linear','extrap');
            %dispAll(:,dC) = dispNow(:,2);
            dispAll(:,dC) = dispNowIntp;
            dC = dC+1;
            fC = 1;
            dispNow = [];
            oneDisp = 0;
            %disp(['Saving disp ', num2str(dC)]);
            
        end
    else
        if(~newRead)
            spInd = find(tLine==' ');
            fNow = tLine(1:(spInd-1));
            sNow = tLine((spInd+1):end);
            dispNow(fC,1:2) = [str2double(fNow),str2double(sNow)];
            fC = fC+1;
            oneDisp = 1;
            
            %disp('check!');
        end
    end
end
fNow = dispNow(:,1);
%disp('check!')

end

