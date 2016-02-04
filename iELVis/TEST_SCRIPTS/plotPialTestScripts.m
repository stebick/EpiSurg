%%
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.overlayParcellation='D';
cfg.showLabels='y';
cfg.title='PT001: Destrieux Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);

%%
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.overlayParcellation='DK';
cfg.title='PT001: DK Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);

%%
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.overlayParcellation='Y7';
cfg.title='PT001: Yeo7 Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);

%%
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.overlayParcellation='Y7';
cfg.elecCoord='n';
cfg.title='PT001: Yeo7 Atlas'; 
cfgOut=plotPialSurf('fsaverage',cfg);
%cfgOut=plotPialSurf('PT001',cfg);