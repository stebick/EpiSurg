% This script runs the brain shift correction commands on this page of the
% wiki:
%  http://episurg.pbworks.com/w/page/101851831/Electrode%20Localization


%% Test Yang-Wang brain shift correction code
makeIniLocTxtFile('PT001');
yangWangElecPjct('PT001');
% Note, 1st grid is 8x8, STP grid is 4 x 5 [1 5 20 16]


%% Test Dykstra brain shift correction code
dykstraElecPjct('PT001');


%% Test Plots: Mgrid on Pial Surf
plotMgridOnPial('PT001',1);


%% Test Plot: Mgrid on Slices
cfg=[]; cfg.printFigs=1;
plotMgridOnSlices('PT001',cfg);