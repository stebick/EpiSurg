function [avgCoords, elecNames, isLeft, avgVids, subVids]=pial2AvgBrain(subj,cfg)
%function [avgCoords, elecNames, isLeft, avgVids, subVids]=pial2AvgBrain(subj,cfg)
%
% This function takes RAS "pial" coordinates (snapped to pial surface)
% and maps it to the corresponding location on the pial surface of
% FreeSurfer's average brain.
%
% Inputs:
%   subj = FreeSurfer subject name
%
% Optional Inputs: passed as fields in a configuration structure
%   plotEm = 1 or 0.  If nonzero, a figure is created illustrating
%            electrode locations on subject and average pial surface.
%            Click on electrodes to see names. Depth electrodes are not
%            shown. {default: 1}
%   elecCoord = N-by-3 numeric array with RAS electrode coordinates. 
%               {default: not used; the function looks into the subject's 
%               Freesurfer folder for electrode coordinate file instead}
%   elecNames = cell array of strings with electrode names, corresponding
%               to the rows of elecCoord. This argument is required 
%               if elecCoord is used. {default: not used; the function
%               looks into the subject's Freesurfer folder for electrode
%               name file instead}
%   isLeft    = N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is on/in the right hemisphere.
%               1 indicates a left hemisphere location. This argument is
%               required if elecCoord is used. Otherwise this information
%               is read from the participant's electrodeNames file.
%   isSubdural= N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is a depth electrode. 1
%               indicates a subdural electrode. This argument is only used
%               if elecNames is used. Otherwise this information is read
%               from the participant's electrodeNames file {default:
%               all electrodes are assumed to be subdural}
%   rmDepths = 1 or 0. If nonzero, depth electrodes are ignored. {default: 1}
%
% Outputs:
%   avgCoords = Electrode coordinates on FreeSurfer avg brain pial surface
%                (RAS coordinates)
%   elecNames = Channel names with the participant's name appended to the
%               beginning (e.g., PT001-Gd1)
%   isLeft    = N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is on/in the right hemisphere.
%               1 indicates a left hemisphere location.
%   avgVids   = Index of subject pial surface vertices corresponding to each electrode
%   subVids   = Index of average pial surface vertices corresponding to each electrode
%
%
% Author:
% David Groppe
% Mehtalab
% March, 2012
%

% History
% 2015-6 Made compatible with Yang brain shift correction algorithm

% parse input parameters in cfg structure and set defaults
if  ~isfield(cfg,'plotEm'),         plotEm = 1;     else    plotEm = cfg.plotEm;            end
if  ~isfield(cfg,'elecCoord'),      elecCoord = []; else    elecCoord = cfg.elecCoord;      end
if  ~isfield(cfg,'elecNames'),      elecNames = []; else    elecNames = cfg.elecNames;      end
if  ~isfield(cfg,'isLeft'),        isLeft = [];   else    isLeft = cfg.isLeft;      end
if  ~isfield(cfg,'isSubdural'),     isSubdural = [];   else    isSubdural = cfg.isSubdural;      end
if  ~isfield(cfg,'rmDepths'),       rmDepths = 1;   else    rmDepths = cfg.rmDepths;      end
checkCfg(cfg,'pial2AvgBrain.m');


% FreeSurfer Subject Directory
fsDir=getFsurfSubDir();
avgDir=fullfile(fsDir,'fsaverage');
subDir=fullfile(fsDir,subj);

% Take care of electrode names, hemisphere, and type
if isempty(elecNames)
    % Import electrode names
    elecFname=fullfile(subDir,'elec_recon',[subj '.electrodeNames']);
    elecInfo=csv2Cell(elecFname,' ',2);
    elecNames=elecInfo(:,1);
    nElec=size(elecInfo,1);
    isLeft=zeros(nElec,1);
    isSubdural=zeros(nElec,1);
    for a=1:nElec,
        if ~strcmpi(elecInfo{a,2},'D')
            isSubdural(a)=1;
        end
        if strcmpi(elecInfo{a,3},'L')
            isLeft(a)=1;
        end
    end
else
    nElec=length(elecNames);
    if isempty(isLeft)
        error('You need to specify cfg.isLeft when using cfg.elecNames');
    else
        if length(isLeft)~=nElec,
            error('elecNames and isLeft do not have the same # of electrodes.');
        end
    end
    if isempty(isSubdural),
        % assume all electrodes are subdural
        isSubdural=ones(nElec,1);
    else
        if length(isSubdural)~=nElec,
            error('isSubdural and isLeft do not have the same # of electrodes.');
        end
    end
end


% Take care of electrode coordinates in participant space
if isempty(elecCoord) % no electrode coordinates have been passed in the function call:
    % Import electrode PIAL coordinates
    coordFname=fullfile(subDir,'elec_recon',[subj '.PIAL']);
    coordCsv=csv2Cell(coordFname,' ',2);
    elecCoord=zeros(nElec,3);
    for a=1:nElec,
        for b=1:3,
            elecCoord(a,b)=str2double(coordCsv{a,b});
        end
    end
else
    if size(elecCoord,1)~=nElec,
        error('Electrode coordinates need to have the same number of rows as electrode names.');
    end
end


% Remove depths (optional)
if universalYes(rmDepths),
    keepIds=find(isSubdural);
    isLeft=isLeft(keepIds);
    elecNames=elecNames(keepIds);
    elecCoord=elecCoord(keepIds,:);
    nElec=length(elecNames);
end


% Get electrode coordinates on average brain
avgCoords=zeros(nElec,3);
avgVids=zeros(nElec,1);
subVids=zeros(nElec,1);
plotCtOffset=0; % counts the # of electrodes that have been displayed
for hemLoop=1:2,
    if hemLoop==1
        % Do left hemisphere elecs
        hemElecIds=find(isLeft);
        hem='l';
    else
        % Do right hemisphere elecs
        hemElecIds=find(~isLeft);
        hem='r';
    end
    
    if ~isempty(hemElecIds)
        fprintf('Working on hemisphere: %s\n',hem);
        pialFname=fullfile(subDir,'surf',[ hem 'h.pial']);
        pial=readSurfHelper(pialFname);
        
        nHemElec=length(hemElecIds);
        n_pial_vert=size(pial.vert,1);
        for a=1:nHemElec,
            df=pial.vert-repmat(elecCoord(hemElecIds(a),:),n_pial_vert,1);
            dst=sum(df.^2,2);
            [~, subVids(hemElecIds(a))]=min(dst);
        end
        
        sphFnameSub=fullfile(subDir,'surf',[hem 'h.sphere.reg']);
        sph=readSurfHelper(sphFnameSub);
        
        sphFnameAvg=fullfile(avgDir,'surf',[hem 'h.sphere.reg']);
        avg_sph=readSurfHelper(sphFnameAvg);
        nAvgVert=length(avg_sph.vert);
        for a=1:nHemElec,
            df=avg_sph.vert-repmat(sph.vert(subVids(hemElecIds(a)),:),nAvgVert,1);
            dst=sum(df.^2,2);
            [~, avgVids(hemElecIds(a))]=min(dst);
        end
        
        avgPialFname=fullfile(avgDir,'surf',[ hem 'h.pial']);
        avgPial=readSurfHelper(avgPialFname);
        for a=1:nHemElec,
            avgCoords(hemElecIds(a),:)=avgPial.vert(avgVids(hemElecIds(a)),:);
        end
        
        
        % Plot results (optional)
        if universalYes(plotEm)
            % LATERAL VIEW
            hFig=figure;
            set(hFig,'position',[360 335 829 360]);
            
            % Plot Electrodes on Avg Brain
            subplot(1,2,1);
            map=[1 1 1]*.7;
            tripatchDG(avgPial,hFig,map);
            shading interp; lighting gouraud; material dull; axis off, hold on
            if strcmp(hem,'r')
                l=light('Position',[1 0 0]);
                view(90,0);
            else
                l=light('Position',[-1 0 0]);
                view(270,0);
            end
            for a=1:nHemElec,
                h=plot3(avgCoords(a+plotCtOffset,1),avgCoords(a+plotCtOffset,2),avgCoords(a+plotCtOffset,3),'r.');
                clickText(h,elecNames{a+plotCtOffset});
                set(h,'markersize',20);
            end
            rotate3d off;
            
            % Plot Electrodes on Participant Brain
            subplot(1,2,2);
            tripatchDG(pial,hFig,map);
            shading interp; lighting gouraud; material dull; axis off, hold on
            if strcmp(hem,'r')
                l=light('Position',[1 0 0]);
                view(90,0);
            else
                l=light('Position',[-1 0 0]);
                view(270,0);
            end
            for a=1:nHemElec
                d=subVids(a+plotCtOffset);
                h=plot3(pial.vert(d,1),pial.vert(d,2),pial.vert(d,3),'r.');
                clickText(h,elecNames{a+plotCtOffset});
                set(h,'markersize',20);
            end
            rotate3d off;
            set(gcf,'name',subj);
            
            % MEDIAL VIEW
            hFig=figure;
            set(hFig,'position',[360 325 829 360]);
            
            % Plot Electrodes on Avg Brain
            subplot(1,2,1);
            map=[1 1 1]*.7;
            tripatchDG(avgPial,hFig,map);
            shading interp; lighting gouraud; material dull; axis off, hold on
            if strcmp(hem,'r')
                l=light('Position',[-1 0 0]);
                view(270,0)
            else
                l=light('Position',[1 0 0]);
                view(90,0);
            end
            for a=1:nHemElec
                h=plot3(avgCoords(a+plotCtOffset,1),avgCoords(a+plotCtOffset,2), ...
                    avgCoords(a+plotCtOffset,3),'r.');
                %                 clickText(h,[elecNames{a} sprintf(' %.3f %.3f %.3f',avgCoords(a,1), ...
                %                     avgCoords(a,2),avgCoords(a,3))]);
                clickText(h,elecNames{a+plotCtOffset});
                set(h,'markersize',20);
            end
            rotate3d off;
            
            % Plot Electrodes on Participant Brain
            subplot(1,2,2);
            tripatchDG(pial,hFig,map);
            shading interp; lighting gouraud; material dull; axis off, hold on
            if strcmp(hem,'r')
                l=light('Position',[-1 0 0]);
                view(270,0)
            else
                l=light('Position',[1 0 0]);
                view(90,0);
            end
            for a=1:nHemElec
                d=subVids(a+plotCtOffset);
                h=plot3(pial.vert(d,1),pial.vert(d,2),pial.vert(d,3),'r.');
                clickText(h,elecNames{a+plotCtOffset});
                set(h,'markersize',20);
            end
            rotate3d off;
            set(gcf,'name',subj);
            
            drawnow;
            plotCtOffset=nHemElec+plotCtOffset;
        end
    end
end

%% Add subject name as prefix to electrode names:
for a=1:nElec,
   elecNames{a}=[subj '-' elecNames{a}];
end