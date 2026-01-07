% this script was written to plot the group velocities from the Stars in
% the arrat, these velocities are then used to create the group velocity
% band for cleaning the correlations of spurious energy
clear; close all;

starPath = 'A:\TestInver\StarsApril\';
starVal = 1:51;
nStar = length(starVal);

figure(1);
hold on;
for i = 1:1:nStar
    fPathFull = [starPath,'Star',num2str(starVal(i)),'\FundUpdApril.txt'];
    A = load(fPathFull);
    plot(A(:,1),A(:,2));
end

hold off;
