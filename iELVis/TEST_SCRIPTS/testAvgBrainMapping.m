%% Test Wiki Example 1
groupAvgCoords=[];
groupLabels=[];
groupIsLeft=[];
subs={'PT001','PT002'};
for a=1:length(subs),
    fprintf('Working on Participant %s\n',subs{a});
    [avgCoords, elecNames, isLeft]=pial2AvgBrain(subs{a},[]);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; elecNames];
    groupIsLeft=[groupIsLeft; isLeft];
end

%% Test Wiki Example 2
cfg=[];
cfg.view='l';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title='PT001-2 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);


%% Same as above but on Inflated average brain
cfg=[];
cfg.view='l';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.surfType='inflated';
cfg.title='PT001-2 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);


%% Test on both hemispheres
groupAvgCoords=[];
groupLabels=[];
groupIsLeft=[];
subs={'TWH035','PT001'};
for a=1:length(subs),
    fprintf('Working on Participant %s\n',subs{a});
    [avgCoords, elecNames, isLeft]=pial2AvgBrain(subs{a},[]);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; elecNames];
    groupIsLeft=[groupIsLeft; isLeft];
end

%% Test omni plot of both hemispheres
cfg=[];
cfg.view='omni';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title='PT001, TWH035 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);

%% Test inflated plot 
cfg=[];
cfg.view='l';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title='PT001, TWH035 on Avg. Brain';
cfg.surfType='inflated';
cfgOut=plotPialSurf('fsaverage',cfg);

%% Test lomni plot 
cfg=[];
cfg.view='lomni';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title='PT001, TWH035 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);

%% Test romni plot 
cfg=[];
cfg.view='romni';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title='PT001, TWH035 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);


%%
disp('Script testAvgBrainMapping.m completed successfully.')