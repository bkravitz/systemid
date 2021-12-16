function ax = setMapAGU()

% ax = axesm('eckert4', 'Frame', 'on', 'Grid', 'on');
ax = axesm('behrmann', 'Frame', 'on', 'Grid', 'on','Origin',[0 -150 0]);
coast = load('coast.mat');
plotm(coast.lat,coast.long,'k');
% landareas = shaperead('landareas.shp','UseGeoCoords',true);
% geoshow(landareas,'FaceColor',[1 1 0.5],'EdgeColor',[.6 .6 .6]);

end