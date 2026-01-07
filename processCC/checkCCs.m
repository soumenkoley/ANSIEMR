% ok so I am writing this code to see all the cross correlation pairs in a
% subarray
% the idea is to select the ones with symmetric nature. Dont use others you
% will estimate the group velocity wong. Stick to the principles, dont be
% greedy :P

clear; close all;
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\subArray11NewCC.mat')
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\AllStnLoc.mat');
load('HdBp3_8Hz.mat');
fSamp = 25;
% now see the cross-correlations one by one

goodCount = 1;
for i = 1:1:length(ccStoreFinal(1,:))
    ccNow = filtfilt(HdBp3_8Hz.Numerator,1,ccStoreFinal(:,i));
    figure(28)
    plot(allStnLoc(:,1),allStnLoc(:,2),'bo','MarkerFaceColor','b','MarkerSize',6);
    hold on;
    stnALoc = rayAttribute(i,1:2); stnBLoc = rayAttribute(i,3:4);
    plot(stnALoc(1,1),stnALoc(1,2),'ro','MarkerFaceColor','r','MarkerSize',6)
    plot(stnBLoc(1,1),stnBLoc(1,2),'ro','MarkerFaceColor','r','MarkerSize',6)
    disp(['Station distance is ',num2str(rayAttribute(i,5)),' m']);
    hold off;
    
    figure(29)
    subplot(2,1,1);
    plot(tArray,ccNow);
    xlim([-10,10]);
    subplot(2,1,2)
    fVec = linspace(0,1,length(ccNow))*fSamp;;
    plot(fVec,abs(fft(ccNow)));
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlim([1,6]);
    
    stnAChar = convertStringsToChars(stnList(rayAttribute(i,7),1));
    stnBChar = convertStringsToChars(stnList(rayAttribute(i,8),1));
    disp(['Stn A is :',stnAChar,' Stn B is :', stnBChar]);
    prompt = 'Do you want to save this?';
    str = input(prompt,'s');
    
    if(str=='y')
        ccGood(:,goodCount) = ccStoreFinal(:,i);
        rayAttrGood(goodCount,:) = rayAttribute(i,:);
        goodCount = goodCount+1;
    end
end
