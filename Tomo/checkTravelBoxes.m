% set up grid first
clear all;
xExtent = 100; yExtent = 100;
dx = 20; dy = 20;

nodeA = [76,34]; 
nodeB = [16,87];

figure(1)
hold on;

for m = 1:1:(xExtent/dx)
    xStart = (m-1)*dx; yStart = 0;
    xEnd = (m-1)*dx; yEnd = 100;
    plot([xStart,xEnd],[yStart,yEnd],'k');
end

for n = 1:1:(yExtent/dy)
    yStart = (n-1)*dy; xStart = 0;
    yEnd = (n-1)*dx; xEnd = 100;
    plot([xStart,xEnd],[yStart,yEnd],'k');
end

plot([nodeA(1,1),nodeB(1,1)],[nodeA(1,2),nodeB(1,2)],'b');
hold off;

[propDetails] = findTravelBoxes(nodeA,nodeB,dx,dy,xExtent,yExtent);