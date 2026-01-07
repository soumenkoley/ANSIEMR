% extract green's function for all possible station pairs under some
% constraints ofcourse, because of the omni-directional nature of noise and
% hence the station paths oriented at particular azimuths are only
% considered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Node Locations Data here

clear;
close all;
tic;

currDate = 04;
endDate = 15;
%addpath('/home/soumen/Dropbox/EinsteinTelescopeSurvey/DataAnalysisPassive1/' )
addpath('C:\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
addPath('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunction\RayPathCompute\')
load('nodeLocations.mat');
load('HdHP1HzLS.mat');
colors = distinguishable_colors(20);

data = [];
fSamp = 250;
fDownSamp = 25;
DSRatio = fSamp/fDownSamp;
maxLag = 40;
startHr = 00; startMin = 00; startSec = 0; %like 5:25:00 in the morning
endHr = 05; endMin = 59; endSec = 0;

startMinWindows = 00;
windowMin = 10; % so in this case the length will be 2 mins, 5:25->5:26 AM
endMinWindows = 59;
tempNormWin = 2; % units in seconds

noiseDirSlope = 80;
noiseDirSlopeSpread = 30;
noiseDirLB = noiseDirSlope-noiseDirSlopeSpread;
noiseDirUB = noiseDirSlope+noiseDirSlopeSpread;

nodeList = [105,94,32,87];
%nodeList = [44,168,63,103,94,87,96,91];
%nodeList = [63,127,119];
%nodeList = [44,138,15,177,33,94];
% nodeList = [44;194;31;102;92;6;181;40;127;174;126;14;171;136;15;159;...
%              131;120;175;110;187;185;81;88;138;8;112;182;146;150;38;149;29;109;...
%               76;177;78;96;35;119;91;24;139;168;63;156;33;36;23;39;196;...
%         113;183;172;25;34;57;54;178;5;84;164;130;32;87;61;186;93;105;103;104;94]';
%nodeList = [24,29,38,94,127,138,139,168,185];
refNodeNumber = 105;

hourGlass = 60;
minuteGlass = 60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hourNo = 1;
lNodeList = length(nodeList);

Sen5Hz = 80;
S5Hz = Sen5Hz*(40000/(40000+1850));
Vref = 2.5; Gain = 16;

% do parts with the sensor locations
refNodeIndex = find(nodeLocations(:,4)==refNodeNumber);
refNodeLatLong = nodeLocations(refNodeIndex,1:2);

% define node distances vector
nodeLocationsCartesian = zeros(length(nodeLocations(:,1)),3);

offset = zeros(lNodeList,1);
figure(1)
hold on;
for i =1:1:lNodeList
    nodeLocationsCartesian(i,1) = nodeLocations(i,3);
    nodeInd = find(nodeLocations(:,4)==nodeList(1,i));
    nodeLocationUse = nodeLocations(nodeInd,1:2);
    [nodeLocationsCartesian(i,2),nodeLocationsCartesian(i,3)] = ...
        calculatedist(nodeLocationUse(1,1:2),refNodeLatLong,6378);
    nodeLocationsCartesian(i,2:3) = nodeLocationsCartesian(i,2:3)*1000;
    offset(i,1) = sqrt(sum(nodeLocationsCartesian(i,2:3).^2,2));
    if(~(isnan(nodeLocations(i,3))))
        plot(nodeLocationsCartesian(i,2),nodeLocationsCartesian(i,3),'*','color','b');
        str1 = ['---',num2str(nodeList(1,i))];
        text(nodeLocationsCartesian(i,2),nodeLocationsCartesian(i,3),str1);
    end
end

% now decide which pairs to consider.
for i = 1:1:(lNodeList-1)
    for j = (i+1):1:lNodeList 
        nodeAInd = find(nodeLocationsCartesian(:,1)==nodeList(1,i));
        nodeBInd = find(nodeLocationsCartesian(:,1)==nodeList(1,j));
        nodeALoc = nodeLocationsCartesian(nodeAInd,2:3);
        nodeBLoc = nodeLocationsCartesian(nodeBInd,2:3);
        raySlope = findSlopeNew(nodeALoc,nodeBloc);
        if((raySlope>noiseDirLB) && (raySlope<noiseDirUB))
            % this good now check ray length
            rayAttribute(rayCounter,1:2) = [nodeList(1,i),nodeList(1,j)];
            rayAttribute(rayCounter,3) = sqrt(sum((nodeALoc-nodeBLoc).^2));
            rayAttribute(rayCounter,4:5) = nodeALoc;
            rayAttribute(rayCounter,6:7) = nodeBLoc;
            figure(1)
            hold on;
            plot(nodeALoc,nodeBLoc,'k');
            hold off;
        end
            
    end
end
hold off;
title('Nodes Deployed but all not good');
psdCorrAvg = zeros(2*maxLag*fDownSamp+1,(lNodeList-1));
superCounter = 1;
while (currDate <= endDate)
    dayStr = num2str(currDate);
    if(currDate < 10)
        dayStr = ['0' dayStr];
    end
    dateStrActual = ['1711' dayStr];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    checkerman = 1;
    hrVals = [];
    for fileNo = 1:1:4
        sT = 0 + (fileNo-1)*6; eT = fileNo*6;
        load(['C:\Dropbox\EinsteinTelescopeSurvey\Beamformer\goodTimes\day',...
            num2str(currDate),'Hr',num2str(sT),'-',num2str(eT),'.mat']);
        hrVals = [hrVals;dayHrInfo];
    end
        
    hrVals(:,2) = hrVals(:,2)-1;
    
    startHrNo = 1;
    while(startHrNo <= length(hrVals(:,1)))
        %% ndo
        firstLetter = 17;
        firstLetter = firstLetter*10^7; %like 170000117
        startHr = hrVals(startHrNo,2);
        if(startHr == endHr)
            endMinWindows = endMin;
        end
        
        lengthTotTrace = (endMinWindows-startMin+1)*minuteGlass*fSamp;
        AllData = zeros(lengthTotTrace,lNodeList);
        outGData5HzDec =  zeros(lengthTotTrace/DSRatio,lNodeList);
        
        ii = 1;
        
        while(ii <= lNodeList)
            %nodeList = (nodeLocationsCartesian(:,1))';
            %pathName = 'D:\LimburgPassiveSurvey1\SEGY\';
            %pathName = 'C:\Dropbox\EinsteinTelescopeSurvey\DataExtract\SEGY\';
            %pathName = '/media/soumen/Soumen/LimburgPassive1/SEGY/';
            pathName = 'D:\LimburgPassive1\SEGY\';
            %superCounter = 1;
            
            dateStr = dateStrActual; %yymmdd
            f1 = firstLetter + nodeList(1,ii);
            f1Str = num2str(f1);
            f2Str = ['0',f1Str];
            fileName = [pathName,f2Str,'\',dateStr,'\',f2Str,'_',...
                dateStr,'_'];
            if(startHr<10)
                startHrStr = ['0',num2str(startHr)];
                %startMinStr = ['0',num2str(startMin)];
                startSecStr = ['0',num2str(startSec)];
            else
                startHrStr = num2str(startHr);
                %startMinStr = num2str(startMin);
                startSecStr = ['0',num2str(startSec)];
            end
            
            if(startMin<10)
                startMinStr = ['0',num2str(startMin)];
            else
                startMinStr = num2str(startMin);
            end
            
            fileNameTot = [fileName,startHrStr,startMinStr,startSecStr,'.sgy'];
            Data = [];
            newMin = startMin;
            
            % make an exception for the last hour of reading data
            if(startHr == endHr)
                hourGlass = endMin+1;
            else
                hourGlass = 60;
            end
            for jj = startMin:1:(hourGlass-1)
                if(exist(fileNameTot,'file'))
                    [subData]=ReadSegyFast(fileNameTot);
                else
                    subData = zeros(fSamp*minuteGlass,1);
                end
                [Data] = [Data;subData];
                
                % disp(['new min is ',num2str(newMin)]);
                newMin = newMin+1;
                if(newMin<10)
                    newMinStr = ['0',num2str(newMin)];
                else
                    newMinStr = num2str(newMin);
                end
                fileNameTot = [fileName,startHrStr,newMinStr,startSecStr,'.sgy'];
            end
            fclose('all');
            %downSampledData = Data;%downsample(Data,5);
            AllData(1:lengthTotTrace,ii) = Data;
            
            AllData(:,ii) = AllData(:,ii)*(Vref/(Gain*2^23));
            
            % Remove mean from data
            meanMat = mean(AllData(:,ii));
            AllData(:,ii) = AllData(:,ii) - meanMat;
            % now decimate the data
            
            [AllData(:,ii)] = deconGF(AllData(:,ii),0.69,5,fSamp,S5Hz);
            outGData5HzDec(:,ii) = decimate(AllData(:,ii),DSRatio);
            disp(['Extracting data for Node ',num2str(nodeList(1,ii)),'...',...
                'Removing DC trend from data, Instrument Response Deconvolved']);
            
            % now time to filter the data
            outGData5HzDec(:,ii) = filtfilt(HdHP1Hz.Numerator,1,outGData5HzDec(:,ii));
            tVec = (1:1:(lengthTotTrace/DSRatio))/fDownSamp;
            % now do temporal normalization to data
            %[outGData5HzDec(:,ii)] = temporalNorm(outGData5HzDec(:,ii),...
                %tempNormWin,fDownSamp);
            [outGData5HzDec(:,ii)] = whitening(outGData5HzDec(:,ii),fDownSamp,[2,10]); 
            %             figure(2)
            %             subplot(2,1,ii)
            %             plot(tVec,outGData5HzDecimated(:,ii));
            %             title(['Node no = ',num2str(nodeList(ii))]);
            %
            %             figure(3)
            %             subplot(2,1,ii)
            %             plot(fVec,abs(whiteData));
            
            ii = ii + 1;
        end
        
        disp(['Data Extraction Complete for hours ',num2str(startHr),...
            ' to ', num2str(startHr+1)]);
        
        disp([' Day Number = ', num2str(currDate)]);
        
        %%
        % Window signals Now
        for nodeNo = 2:1:lNodeList
            disp('Windowing Signals and obtaining correlations');
            disp(['Window length is ',num2str(windowMin),' mins']);
            windowNo = 1;
            totalWinsHourly = (endMinWindows - startMin + 1)/(windowMin);
        
            psdCorrTempAvg = 0;
            while(windowNo <= totalWinsHourly)
                sInd = (windowNo-1)*windowMin*60*fDownSamp+1;
                eInd = windowNo*fDownSamp*windowMin*60;
                smallData = outGData5HzDec(sInd:eInd,:);
                [psdCorrTemp,lagVal] = xcorr(smallData(:,nodeNo),smallData(:,1),...
                    maxLag*fDownSamp,'coeff');
                psdCorrTempAvg =  psdCorrTemp + psdCorrTempAvg;
                windowNo = windowNo+1;
            end
            psdCorrAvg(:,(nodeNo-1)) = psdCorrAvg(:,(nodeNo-1))+...
                                        psdCorrTempAvg/totalWinsHourly;
            disp(['Computed correlation ',num2str(startHr),...
            ' to ', num2str(startHr+1),'!!!!!!!!!!']);
        end
        startHrNo = startHrNo + 1;
        startMin = 00;
        superCounter = superCounter+1;
    end
    % normalize the pwelch output
    %     pxxAvg = pxxAvg/(hourNo-1);
    startHr = 00;
    currDate = currDate +1;
end
psdCorrAvg = psdCorrAvg/(superCounter-1);
%now differentiate it wrt time
psdCorrAvg = diff(psdCorrAvg)*fDownSamp;
dummyElement = zeros(1,(lNodeList-1));
psdCorrAvg = [dummyElement;psdCorrAvg];
figure(2)
% plot(lagVal/fDownSamp,psdCorrAvg);
% hold on;
% plot(lagVal/fDownSamp,mean(psdCorrAvg,2),'r','LineWidth',1);
plotseis(psdCorrAvg,lagVal/fDownSamp,offset(2:end,1),1);
%surf(offset(2:end,1),lagVal/fDownSamp,psdCorrAvg);
%shading interp;
toc;