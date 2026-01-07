% this script was written as a test script for checking how the new
% cross-correlations look when compared to the old ones
% Basically it assembles all the cross-correlations

clear; %close all;
n1 = "XUKWA"; n2 = "0RPXA";

dayVec = 320; % this is Nov 15
yrVec = 2020;
hrVec = 0:1:23;

fPathBase = 'A:\Day';

fSamp = 25;

ccFull = zeros(1001,1);
for dayNo = 1:1:length(dayVec)
    for hrNo = 1:1:length(hrVec)
        fPath  = [fPathBase,num2str(dayVec(dayNo)),'\Hr',...
            num2str(hrVec(hrNo)),'\'];
        fNameA = n1;
        refStnChar = char(n1); endStnChar = char(n2);
        if(exist([fPath,refStnChar,'.mat']))
            load([fPath,refStnChar,'.mat']);
            stnEnd = string(stnEnd);
            endStnInd = find(stnEnd==n2);
            if(isempty(endStnInd))
                % then load otherstation file
                if(exist([fPath,endStnChar,'.mat']))
                    load([fPath,endStnChar,'.mat']);
                    stnEnd = string(stnEnd);
                    refStnInd = find(stnEnd==n1);
                    if(~isempty(refStnInd))
                        ccNow = flipud(ccStore(:,refStnInd));
                    else
                        disp('Stn pair not found');
                    end
                else
                    disp('Stn pair not found');
                end
            else
                ccNow = ccStore(:,endStnInd);
            end
        else
            if(exist([fPath,endStnChar,'.mat']))
                load([fPath,endStnChar,'.mat']);
                stnEnd = string(stnEnd);
                refStnInd = find(stnEnd==n1);
                if(~isempty(refStnInd))
                    ccNow = flipud(ccStore(:,refStnInd));
                else
                    disp('Stn pair not found');
                end
            else
                disp('Stn Pair not found!');
            end
        end
        ccFull = ccFull + ccNow;
    end
end

tArray = (-500:1:500)/fSamp;
figure(1)
subplot(2,1,1)
plot(tArray,ccFull);
xlim([-20,20])

fftCC = fft(ccFull);
fVec = linspace(0,1,length(ccFull(:,1)))*fSamp;
subplot(2,1,2)
loglog(fVec,abs(fftCC));
xlim([0.5,5]);

        