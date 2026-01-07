% this script makes the plots of the group velocity maps
h1 = openfig('D:\LimburgSpatialResolution\Fig10.fig','reuse'); % open figure

ax1 = gca; % get handle to axes of figure

h2 = openfig('D:\LimburgSpatialResolution\Fig11.fig','reuse');

ax2 = gca;

h3 = openfig('D:\LimburgSpatialResolution\Fig12.fig','reuse');

ax3 = gca;


% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

h4 = figure; %create new figure

s1 = subplot(1,3,1); %create and get handle to the subplot axes

s2 = subplot(1,3,2);

s3 = subplot(1,3,3);

fig1 = get(ax1,'children'); %get handle to all the children in the figure

fig2 = get(ax2,'children');

fig3 = get(ax3,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes

copyobj(fig2,s2);

copyobj(fig3,s3);