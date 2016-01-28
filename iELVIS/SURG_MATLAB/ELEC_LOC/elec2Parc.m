%  elec2parc() - Returns the name of the nearest cortical area according
%                to FreeSurfer parcellations.
%
% Usage:
%  >>[elecParc, hem]=elec2parc(subj,format_names);
%
% Required Input:
%  subj      -Name of the subject's freesurfer directory (full path not
%             needed)
%
% Optional Input:
%  format_names -'y' or 'n'. If 'y' underscores are removed from electrode
%                names and "Grid" is converted to "G". {default: 'n'}
%
% Output:
%   elecParc - 2D cell array containing electrode names and their
%                 assigned cortical area
%   hem         - 1D cell array indicating the hemisphere each electrode
%                 is in
%
% Example:
% >> elecParc=elec2parc('AnRo');
%
% Author: David M. Groppe

%
% Future Work:
% -I should have a function called neat_labels.m that will format the DK
% atlas names to be more kind to the eye. It might help to incoporate that.

function elecParc=elec2Parc(subj,atlas)

if nargin<2,
    atlas='DK';
end

fsDir=getFsurfSubDir();

% ??
warning('Groppe says: modifiy this to use add_depth_aseg_labels2ecog.m to take care of depths');


% Folder with surface files
surfaceFolder=fullfile(fsDir,subj,'surf');
if ~isempty(surfaceFolder) && (surfaceFolder(end)~='/')
    surfaceFolder=[surfaceFolder '/'];
end

% Folder with cortical parcellation files
labelFolder=fullfile(fsDir,subj,'label');
if ~isempty(labelFolder) && (labelFolder(end)~='/')
    labelFolder=[labelFolder '/'];
end

% Import electrode locations
pvoxFname=fullfile(fsDir,subj,'elec_recon',sprintf('%s.PIAL',subj));
pvoxCoordStr=csv2Cell(pvoxFname,' ',2);
nElec=size(pvoxCoordStr,1);
pvoxCoord=zeros(nElec,3);
for a=1:nElec,
    for b=1:3,
        pvoxCoord(a,b)=str2num(pvoxCoordStr{a,b});
    end
end

% Import electrode labels
labelFname=fullfile(fsDir,subj,'elec_recon',sprintf('%s.electrodeNames',subj));
elecLabels=csv2Cell(labelFname,' ',2);

elecParc=cell(nElec,2);
hem=[];
for hemLoop=1:2,
    if hemLoop==1
        hem='L';
    else
        hem='R';
    end
    
    %% Are there any electrodes in this hemisphere?
    elecIdsThisHem=findStrInCell(hem,elecLabels(:,3));
    nElecThisHem=length(elecIdsThisHem);
    if nElecThisHem,
        
        %% READ SURFACE
        surfFname=fullfile(surfaceFolder,[lower(hem) 'h.pial']);
        [cort.vert, cort.tri]=read_surf(surfFname);
        nVertex=length(cort.vert);
        
        %% Get cortical parcellation
        switch upper(atlas)
            case 'DK'
                parcFname=fullfile(labelFolder,[lower(hem) 'h.aparc.annot']);
                [~, label, colortable]=read_annotation(parcFname);
                %[averts,label,colortable]=read_annotation(parcFname);
            case 'D'
                parcFname=fullfile(labelFolder,[lower(hem) 'h.aparc.a2009s.annot']);
                [~, label, colortable]=read_annotation(parcFname);
            case 'Y7'
                parcFname=fullfile(labelFolder,[lower(hem) 'h_Yeo2011_7Networks_N1000.mat']);
                load(parcFname);
            case 'Y17'
                parcFname=fullfile(labelFolder,[lower(hem) 'h_Yeo2011_17Networks_N1000.mat']);
                load(parcFname);
            otherwise
                error('Unrecognized value of atlas argument.')
        end
        
        for elecLoop=1:nElecThisHem,
            elecParc{elecIdsThisHem(elecLoop),1}=elecLabels{elecIdsThisHem(elecLoop),1};
            
            % Go through and set depth electrode assignments to depth:
            if elecLabels{elecIdsThisHem(elecLoop),2}=='D'
                elecParc{elecIdsThisHem(elecLoop),2}='Depth'; % ?? fix this later
            else
                % Find closest vertex
                [~, minId]=min(sum( (repmat(pvoxCoord(elecIdsThisHem(elecLoop),:),nVertex,1)-cort.vert).^2,2 ));
                
                % Grab parcellation label for that vertex
                elecParc{elecIdsThisHem(elecLoop),2}=colortable.struct_names{find(colortable.table(:,5)==label(minId))};
            end
        end
        
    end
end

