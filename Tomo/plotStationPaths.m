clear;
load('delaysAllFirstTrial.mat');

load('ActiveNodeLocsVirgo.mat');
refNode = 143;
refNodeInd = find(ActiveNodeListVirgo(:,1)==refNode);
refNodeLoc = ActiveNodeListVirgo(refNodeInd,5:6);
nodeLocCartesian = zeros(length(ActiveNodeListVirgo),3);
figure(1)
hold on;
for i = 1:1:length(ActiveNodeListVirgo)
    nodeLocCartesian(i,1) = ActiveNodeListVirgo(i,1);
    [nodeLocCartesian(i,2),nodeLocCartesian(i,3)] = calculatedist(ActiveNodeListVirgo(i,5:6),...
                                                refNodeLoc,6378);
    plot(nodeLocCartesian(i,2),nodeLocCartesian(i,3),'bo','MarkerSize',8,'MarkerFaceColor','b');
end

for i = 1:1:333
    nodeA = delaysAll(i,3); nodeB = delaysAll(i,4);
    nodeAInd = find(nodeLocCartesian(:,1) == nodeA);
    nodeBInd = find(nodeLocCartesian(:,1) == nodeB);
    nodeALoc = nodeLocCartesian(nodeAInd,2:3);
    nodeBLoc = nodeLocCartesian(nodeBInd,2:3);
    plot([nodeALoc(1,1),nodeBLoc(1,1)],[nodeALoc(1,2),nodeBLoc(1,2)],'k');

end
hold off;




