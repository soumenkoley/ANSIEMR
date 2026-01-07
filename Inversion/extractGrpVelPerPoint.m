% this script was written to extract the group velocity per grid point
% from the group velocity maps
% also the xyLimits boundaries are used to make points NaN
% we will prepare the group velocity per grid point for two scenarios
% a) only the group velocity extracted at that point
% b) use the group velocities from the neighboring grid points as well

% a mat file with list of points which lie inside the boundaries
% output will be a mat file which will have the group velocity,
% and the grid coordinate
% also the dispersion curves will be written to a disp file
% to be used in the gptarget

clear; close all;
fExt = 1.0:0.05:5.0;
nFreq = length(fExt);
fSt = 1.0;
fStInd = find(fExt>=fSt,1,'first');
fAll = (fExt(fStInd:end))';
c0 = 4000;
k0 = 2*pi*fSt/c0;

fPath = ['B:\LimburgBigSurvey1CC-Pair\VelMapsApril\'];
currDir = pwd;
fPathSave = ['B:\LimburgBigSurvey1CC-Pair\VelGridInvertApril\'];
% load all the x,y coordinate vectors
load('AllXYGridPoints.mat');
%xCoord, yCoord all loaded 
nX = length(xCoord);
nY = length(yCoord);
% the model matrix will be (nY x nX) matrix
% travel across x, then y to create all set of points

gC = 1;
for i = 1:1:nY
    for j = 1:1:nX
        xyAllGrid(gC,1) = xCoord(j);
        xyAllGrid(gC,2) = yCoord(i);
        gC = gC+1;
    end
end

totPoints = length(xyAllGrid(:,1));
% load the boundaries of the polygon, grid points lying within will be
% accepted

%load('xyLimits.mat');

% gInd is a column vector of logicals
%gInd = inpolygon(xyAllGrid(:,1),xyAllGrid(:,2),xyPts(:,1),xyPts(:,2));

% load all the group velocity maps per frequency one by one
% populate the matrix with nCol = totGrid points, nRow = nFreq

fullVelMat = NaN*ones(nFreq,totPoints); 
for i = 1:1:nFreq
    fName = ['grpVel',num2str(fExt(i)),'.mat'];
    fNameFull = [fPath,fName];
    load(fNameFull);
    % change the matrix to column vector
    % each row becomes a column and stacked over each other
    tempMat = mFinalMatNew';
    tempMat = 1./tempMat;
    velVec = tempMat(:);
    fullVelMat(i,:) = velVec';
end

% check for NaN to get gInd
for i = 1:1:totPoints
    nanInd = isnan(sum(fullVelMat(:,i)));
    if(nanInd)
        gInd(i,1) = 0;
    else
        gInd(i,1) = 1;
    end
end
%%
% plot to check how good the extraction was
% lets say we want to extract 3 Hz,
% fInd = find(fExt>=1.6,1,'first');
% 
% minV = min(fullVelMat(fInd,:));
% maxV = max(fullVelMat(fInd,:));
% allCol = jet(100);
% figure(1);
% hold on;
% for i = 1:1:totPoints
%     if(gInd(i,1))
%         colInd = floor((fullVelMat(fInd,i)-minV)/(maxV-minV)*100)+1;
%         if(colInd>100)
%             colInd = 100;
%         end
%         plot(xyAllGrid(i,1),xyAllGrid(i,2),'ko','MarkerFaceColor',allCol(colInd,:));
%     end
% end
% hold off;
% colorbar; colormap('jet'); caxis([minV,maxV]);
% checked extraction was good
%% 
%%
% a folder will be created for every grid point
% two mat files will be saved, one with just a single group vel curve
% the other will also adjacent group curves

% a check will be done to verify that the right neighbours are selected
% check has been done
for i = 1:1:totPoints
    fPathSaveFull = [fPathSave,'grid',num2str(i),'\'];
    
    % check if not nan values then save
    if(gInd(i,1))
        
        % check if it exists
        if(~exist(fPathSaveFull))
            % create it
            mkdir(fPathSaveFull);
        end
        %figure(1);
        %plot(xyAllGrid(gInd,1),xyAllGrid(gInd,2),'b*');
        %hold on;
        
        fNameVel = [fPathSaveFull,'grpVel.mat'];
        velNow = fullVelMat(:,i);
        xyNow = xyAllGrid(i,:);
        % saving only one dispersion curve, fExt, and the grid coord
        save(fNameVel,'fExt','velNow','xyNow');
        velSaveTxt = [fExt',velNow];
        save([fPathSaveFull,'velNow.txt'],'velSaveTxt','-ascii');
        velNowSmooth = smooth(velNow);
        
        % get the phase velocity as well
        % fAll = fAll';
        % integrate the group velocity with respect to frequency
        k = 2*pi*cumtrapz(fAll(1:end,1),1./(1*velNowSmooth(fStInd:end)));
        kC = k0-k(1,1);
        k = kC+k;
        phVelOut = 2*pi*fAll(1:end,1)./k;
        
        velSaveTxtSmooth = [fExt',velNowSmooth];
        phVelSaveTxtSmooth = [fAll,phVelOut];
        save([fPathSaveFull,'velNowSmooth.txt'],'velSaveTxtSmooth','-ascii')
        save([fPathSaveFull,'phVelNowSmooth.txt'],'phVelSaveTxtSmooth','-ascii')
        
        % now saving for the case with neighboring group curves
        % first decide which grid points are nearby, and if they are
        % not Nan
        % top and bottom
        gr1 = i+nX; gr2 = i-nX;
        % left and right
        gr3 = i-1; gr4 = i+1;
        kk = 1;
        xyNbour = [];
        vN = [];
        if(gr1<=totPoints)
           % it means it has not exceeded the limits
           % check for NaN
           if(~isnan(fullVelMat(1,gr1)))
               vN(:,kk) = fullVelMat(:,gr1);
               xyNbour(kk,:) = xyAllGrid(gr1,:);
               %plot(xyNbour(kk,1),xyNbour(kk,2),'ro');
               kk = kk+1;
           end
        end
        
        if(gr2>0)
           % it means it has not exceeded the limits
           if(~isnan(fullVelMat(1,gr2)))
               vN(:,kk) = fullVelMat(:,gr2);
               xyNbour(kk,:) = xyAllGrid(gr2,:);
               %plot(xyNbour(kk,1),xyNbour(kk,2),'ro');
               kk = kk+1;
           end
        end
        
        % get the quotients
        qOrig = floor((i/(nX+1)));
        q3 = floor((gr3/(nX+1))); q4 = floor((gr4/(nX+1)));
        
        if(q3==qOrig) % means they are in the same row
            if(~isnan(fullVelMat(1,gr3)))
               vN(:,kk) = fullVelMat(:,gr3);
               xyNbour(kk,:) = xyAllGrid(gr3,:);
               %plot(xyNbour(kk,1),xyNbour(kk,2),'ro');
               kk = kk+1; 
            end
        end
        
        if(q4==qOrig) % means they are in the same row
           if(~isnan(fullVelMat(1,gr4)))
               vN(:,kk) = fullVelMat(:,gr4);
               xyNbour(kk,:) = xyAllGrid(gr4,:);
               %plot(xyNbour(kk,1),xyNbour(kk,2),'ro');
               kk = kk+1; 
            end
        end
        
        %disp('Check!');
        %clf;
        
        % save the group velocities for the neighboring scenario
        fNameVelNbour = [fPathSaveFull,'grpVelNbour.mat'];
        save(fNameVelNbour,'fExt','velNow','xyNow','vN','xyNbour');
        
    end
end
