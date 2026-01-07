% thi script was written to pick the group travel times from a 
% FTAN plot
% don't use clear
figHand = figure(202);
fExt = fExt; % already a column vector
tExt = tOut';

disp('Pick the upper limit travel times!')
[uX,uY]= getpts(figHand);
uYIntp = interp1(uX,uY,fExt);

disp('Pick the lower limit travel times!')
[lX,lY]= getpts(figHand);
lYIntp = interp1(lX,lY,fExt);

% now pick the max in this time band
count = 1;
velPick = [];
for i = 1:1:length(fExt)
    if(~isnan(uYIntp(i,1)))
        if(~isnan(lYIntp(i,1)))
            uInd = find(tExt>=uYIntp(i,1),1,'first');
            lInd = find(tExt>=lYIntp(i,1),1,'first');
            [~,maxInd] = max(SAOut(uInd:lInd,i));
            velPick(count,1) = fExt(i,1);
            velPick(count,2) = tExt(uInd+maxInd-1,1);
            count = count +1;
        end
    end
end

figure(figHand);
hold on;
velPick(:,2) = smooth(velPick(:,2));
plot3(velPick(:,1),velPick(:,2),3000*ones(length(velPick),1),'k','LineWidth',2);

[tFinal(:,1),tFinal(:,2)] = getpts();
velFinal = [];
velFinal(:,1) = tFinal(:,1);
velFinal(:,2) = distOut./tFinal(:,2);
tFinal = [];
velFinal(:,2) = smooth(velFinal(:,2));
figure(403);
hold on;
plot(grpCurve(:,1),grpCurve(:,2),'k','LineWidth',2);
plot(velFinal(:,1),velFinal(:,2),'m','LineWidth',2);
hold off;


% final pick using getpts
figure(404);
hold on;
plot(velFinal(:,1),velFinal(:,2),'k--');
hold off;

% convert velocity to depth-velocity
depth = velFinal(:,2)./(3*velFinal(:,1));
figure(405);
plot(velFinal(:,2),depth)