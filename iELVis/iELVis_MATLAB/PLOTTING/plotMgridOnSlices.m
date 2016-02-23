function plotMgridOnSlices(fsSub,cfg)
% function plotMgridOnSlices(fsSub,cfg)
%
% Creates a figure illustrating the location of each electrode in an mgrid
% file in a sagittal, coronal, and axial slice and indicates which part of
% the brain it is in.
%
% Required Inputs:
%  fsSub - Patient's freesurfer directory name
%
% Optional cfg parameters:
%  mgridFname - mgrid filename and path. If empty, name is assumed to be fsSub.mgrid. 
%  fullTitle - If 1, the mgrid and mri voxel coordinates are displayed in
%              the figure title along with the electrode name and anatomical 
%              location. {default: 0}
<<<<<<< HEAD
=======
%  markerSize - The size of the dot in each slice used to represent an
%              electrode's location. {default: 30}
>>>>>>> epiSurg/master
%  cntrst    - 0< number <=1 The lower this number the lower the brighter
%              the image (i.e., the lower the voxel value corresponding to 
%              white). {default: 0.5}
%  pauseOn   - If 1, Matlab pauses after each figure is made and waits for
%              a keypress. {default: 0}
%  printFigs - If 1, each figure is output to an eps file. {default: 0}

%
% Examples:
%  %Specify mgrid file and do NOT print
%  cfg=[];
%  cfg.mgridFname='/Applications/freesurfer/subjects/TWH001/elec_recon/TWH001_bi.mgrid';
%  plotMgridOnSlices('TWH001',cfg);
%
%  %Use FreeSurfer file structure and print
%  cfg=[];
%  cfg.printFigs=1;
%  plotMgridOnSlices('TWH011',cfg);
%
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto

% Future work:
% Add option for fsurf anatomy colors?

if ~isfield(cfg,'mgridFname'),    mgridFname=[];    else mgridFname=cfg.mgridFname; end
<<<<<<< HEAD
if ~isfield(cfg,'fullTitle'),    fullTitle=0;       else fullTitle=cfg.fullTitle; end
=======
if ~isfield(cfg,'fullTitle'),     fullTitle=0;      else fullTitle=cfg.fullTitle; end
if ~isfield(cfg,'markerSize'),    markerSize=30;    else markerSize=cfg.markerSize; end
>>>>>>> epiSurg/master
if ~isfield(cfg,'cntrst'),    cntrst=.5;          else cntrst=cfg.cntrst; end
if ~isfield(cfg,'pauseOn'),    pauseOn=0;          else pauseOn=cfg.pauseOn; end
if ~isfield(cfg,'printFigs'),    printFigs=0;          else printFigs=cfg.printFigs; end
checkCfg(cfg,'plotMgridOnSlices.m');

<<<<<<< HEAD
=======

>>>>>>> epiSurg/master
% FreeSurfer Subject Directory
fsdir=getFsurfSubDir();

mriFname=[fsdir '/' fsSub '/mri/brainmask.mgz'];
if ~exist(mriFname,'file')
   error('File %s not found.',mriFname); 
end
mri=MRIread(mriFname);
%mri.vol is ILA (i.e., S->I, R->L, P->A)
mx=max(max(max(mri.vol)))*cntrst;
mn=min(min(min(mri.vol)));
sVol=size(mri.vol);

% Load mgrid
% if strcmpi(mgridFname,'l') || strcmpi(mgridFname,'r')
%     [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub,mgridFname);
% else
%     [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(mgridFname); % mgrid coords are LIP
% end
if isempty(mgridFname)
    [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub);
end
nElec=length(elecLabels);
elecMatrix=round(elecMatrix);
xyz=zeros(size(elecMatrix));
xyz(:,1)=elecMatrix(:,2);
xyz(:,2)=elecMatrix(:,1);
xyz(:,3)=sVol(3)-elecMatrix(:,3);

depthElecs=zeros(nElec,1);
for a=1:nElec,
    if strcmpi(elecLabels{a}(2),'D')
        depthElecs(a)=1;
    end
end
    
for elecId=1:nElec,
    if depthElecs(elecId)
        figId=figure();
        set(figId,'position',[78 551 960 346],'paperpositionmode','auto');
        
        hm=zeros(1,3);
        figure(figId); clf;
        colormap gray;
        %subplot(131);
        wdth=.35;
        wDelt=.33;
        xStart=-.005;
        yStart=.03;
        ht=.9;
        axes('position',[xStart yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,xyz(elecId,2),:)),[mn mx]);
        axis square;
        set(gca,'xdir','reverse');
        hold on;
        hm(1)=plot(xyz(elecId,3),xyz(elecId,1),'r.');
<<<<<<< HEAD
        set(hm(1),'color',elecRgb(elecId,:));
=======
        set(hm(1),'color',elecRgb(elecId,:),'markersize',markerSize);
>>>>>>> epiSurg/master
        %find image limits
        mxX=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],2);
        mxY=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(1)/2),find(mxY==0)));
        limYb=min(intersect((sVol(1)/2:sVol(1)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax, 
            axis([tempMin tempMax tempMin tempMax]); 
        end
        set(gca,'xtick',[],'ytick',[]);
        
        %subplot(132);
        axes('position',[xStart+wDelt yStart wdth ht]);
        imagesc(squeeze(mri.vol(xyz(elecId,1),:,:)),[mn mx]);
        axis square;
        hold on;
        hm(2)=plot(xyz(elecId,3),xyz(elecId,2),'r.');
<<<<<<< HEAD
        set(hm(2),'color',elecRgb(elecId,:));
=======
        set(hm(2),'color',elecRgb(elecId,:),'markersize',markerSize);
>>>>>>> epiSurg/master
        %find image limits
        mxX=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],2);
        mxY=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[],'xdir','reverse');
<<<<<<< HEAD
        
=======
        set(hm(2),'color',elecRgb(elecId,:),'markersize',markerSize);
>>>>>>> epiSurg/master
        
        %subplot(133);
        axes('position',[xStart+wDelt*2 yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,:,xyz(elecId,3))),[mn mx]);
        axis square;
        hold on;
        hm(3)=plot(xyz(elecId,2),xyz(elecId,1),'r.');
<<<<<<< HEAD
        set(hm(3),'color',elecRgb(elecId,:));
=======
        set(hm(3),'color',elecRgb(elecId,:),'markersize',markerSize);
>>>>>>> epiSurg/master
        %find image limits
        mxX=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],2);
        mxY=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[]);
        
        anatLabel=vox2Seg(xyz(elecId,:),fsSub);
  
        % Remove first 3 characters that indicate hemisphere and electrode
        % type
        formattedLabel=elecLabels{elecId}(4:end);
        formattedLabel=rmChar(formattedLabel,'_'); % remove underscore between electrode stem and #
        
        if universalYes(fullTitle)
            ht=textsc2014([formattedLabel '; mgrid coords(' num2str(elecMatrix(elecId,:)-1) '); fsurf coords(' num2str(xyz(elecId,:)) '); ' anatLabel], ...
                'title');
            set(ht,'fontsize',14,'fontweight','bold');
        else
            ht=textsc2014([formattedLabel '; Anatomical Location: ' anatLabel], ...
                'title');
            set(ht,'fontsize',16,'fontweight','bold');
        end
        set(ht,'position',[.5 .97 0]);
        
        if universalYes(printFigs)
<<<<<<< HEAD
            for a=1:3,
                set(hm(a),'markersize',14);
            end
=======
>>>>>>> epiSurg/master
            % Make sure PICS directory exists
            erPath=fullfile(fsdir,fsSub,'elec_recon');
            outPath=fullfile(erPath,'PICS');
            if ~exist(outPath,'dir')
                dirSuccess=mkdir(outPath);
                if ~dirSuccess,
                    error('Could not create directory %s',dirSuccess);
                end
            end
            
<<<<<<< HEAD
=======
            drawnow;
>>>>>>> epiSurg/master
            figFname=fullfile(outPath,sprintf('%s_%sSlices',fsSub,elecLabels{elecId}));
            fprintf('Exporting figure to %s\n',figFname);
            %print(figId,figFname,'-depsc');
            print(figId,figFname,'-djpeg');
        end
        
        if universalYes(pauseOn)
            fprintf('Paused. Press any key for next electrode.\n');
            pause;
        end
    end
end
fprintf('Done showing all electrodes.\n');