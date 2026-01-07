% this script was writte to asseble the frequency domain cross-correlations
% generated from the python code
% Cross-correlation with respect to a base station is assembled
% Assembling refers to averaging over certain number of days
% remember that the the good 125 stations were used for the analysis out of
% all the 183 stations
% the cross-correlation also includes self correlation which is unity
% throughout
clear; close all;
fPath = 'D:\FreqDomainCC\Nov';
refStn = 'Y8MAA';
load('D:\FreqDomainCC\Nov20\Nov20-Y8MAA-CC.mat');
stnList = string(stnNames);

stnSee = ["Y7L9A";"Y7L9B";"Y6MAA";"Y0MBA";"YYMLA";"ZSMIA";"YHLWA";"X4MBA";"XFMWA";"XON9A";...
    "XNOYA";"ZGPDA";"1OP6A";"1PO5A";"0NQPA";"WUJ2A"];

dayVec = [13,14,17,18,19,20];

ccXXAvg = zeros(length(fCC),length(ccXX(1,:)));
ccYYAvg = zeros(length(fCC),length(ccXX(1,:)));
ccXYAvg = zeros(length(fCC),length(ccXX(1,:)));
for i = 1:1:length(dayVec)
    if(dayVec(i)<10)
        dayStr = ['0',num2str(dayVec(i))];
    else
        dayStr = num2str(dayVec(i));
    end
    
    fPathFull = [fPath,dayStr,'\Nov',dayStr,'-',refStn,'-CC.mat'];
    
    load(fPathFull);
    ccXXAvg  = ccXXAvg + ccXX;
    ccYYAvg  = ccYYAvg + ccYY;
    ccXYAvg  = ccXYAvg + ccXY;
end

refStnInd = find(stnList==string(refStn),1,'first');
refStnLoc = nodeLocs(refStnInd,1:2);
figure(1)
hold on;

colDots = jet(length(stnSee));
plot(nodeLocs(:,1),nodeLocs(:,2),'bo','MarkerFaceColor','b','MarkerSize',6);
for i = 1:1:length(stnSee)
    stnSeeInd = find(stnList==stnSee(i),1,'first');
    nodeLocsNow = nodeLocs(stnSeeInd,1:2);
    plot(nodeLocsNow(1,1),nodeLocsNow(1,2),'ko','MarkerFaceColor',colDots(i,:),'MarkerSize',8);
    distNow(i,1) = sqrt((refStnLoc(1,1)-nodeLocsNow(1,1))^2 + ...
        (refStnLoc(1,2)-nodeLocsNow(1,2))^2);
    ccToPlot(:,i) = real(ccXYAvg(:,stnSeeInd)./sqrt(ccXXAvg(:,stnSeeInd).*ccYYAvg(:,stnSeeInd)));
end

nFigs = ceil(length(stnSee)/6);

ccNo = 1;
exitPlot = 0;
for i = 1:1:nFigs
    figure(10+i)
    
    for sPlotNo = 1:1:6
        subplot(2,3,sPlotNo)
        semilogx(fCC,ccToPlot(:,ccNo),'color',colDots(ccNo,:));
        title(['Dist = ',num2str(distNow(ccNo,1)),' m']);
        xlim([0.15,20]);
        ccNo = ccNo+1;
        if(ccNo>(length(stnSee)))
            exitPlot = 1;
            break;
        end
    end
    if(exitPlot==1)
        break;
    end
        
end

% now try to plot the peak correlation
fPeak = 2.08;
fPeakInd = find(fCC>=fPeak,1,'first');
Z = real(ccXYAvg./sqrt(ccXXAvg.*ccYYAvg));
for i = 1:1:length(stnList)
    ccPeak(i,1) = Z(fPeakInd,i);
end

%minCC = min(ccPeak)-0.1;
minCC = -0.1; maxCC = 0.1;
%maxCC = max(ccPeak)+0.1;
colLev = 100;
colMap = jet(colLev);
figure(4)
hold on;

for i = 1:1:length(stnList)
    if(ccPeak(i,1)>0.09)
        ccPeak(i,1)=0.09;
    end
    if(ccPeak(i,1)<-0.09)
        ccPeak(i,1)=-0.09;
    end
    colInd = floor((ccPeak(i,1)-minCC)/(maxCC-minCC)*colLev);
    
    if(~isnan(colInd))
        plot(nodeLocs(i,1),nodeLocs(i,2),'o','MarkerFaceColor',colMap(colInd,:),'MarkerSize',8);
    end
end
colorbar;
colormap('jet');
caxis([minCC,maxCC])
hold off;

    
