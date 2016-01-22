%  elec2parc() - Returns the name of the nearest cortical area according
%                to FreeSurfer parcellations.
%
% Usage:
%  >>elec_assign=elec2parc(subj,side);
%
% Required Inputs:
%  subj       -Name of the subject's freesurfer directory (full path not
%              needed)
%  side      -'r' or 'l', the hemisphere over which the electrodes are
%             implanted
%
% Output:
%   elec_assign - 2D cell array containing electrode names and their
%   assigned cortical area
%
% Example:
% >> elec_assign=elec2parc('AnRo','r');
%
% Author: David Groppe
% Mehtalab
% Sept. 2012 

function elec_assign=elec2parc(subj,side)

fs_dir=getenv('SUBJECTS_DIR');
surftype='pial';
verblevel=2;

% Folder with surface files
surfacefolder=fullfile(fs_dir,subj,'surf');
if ~isempty(surfacefolder) && (surfacefolder(end)~='/')
    surfacefolder=[surfacefolder '/'];
end

% Folder with cortical parcellation files
labelfolder=fullfile(fs_dir,subj,'label');
if ~isempty(labelfolder) && (labelfolder(end)~='/')
    labelfolder=[labelfolder '/'];
end

% Get side of brain to show
side=lower(side);
if side~='l' && side~='r'
    error('''side'' argument needs to be ''r'' or ''l''.');
end


%% READ SURFACE
if side == 'r'
    [cort.vert cort.tri]=read_surf([surfacefolder 'rh.' surftype]);
else
    [cort.vert cort.tri]=read_surf([surfacefolder 'lh.' surftype]);
end

%% Get cortical parcellation
annot_fname=[labelfolder side 'h.aparc.annot'];
[averts,albl,actbl]=read_annotation(annot_fname);
neat_labels=format_annot_labels(actbl.struct_names);
clear averts;
fvcdat=ones(length(albl),3);
ulvals=intersect(unique(albl),unique(actbl.table(:,5)));
for ac=1:length(ulvals),
    fvcdat(albl==ulvals(ac),1:3)=repmat(actbl.table((actbl.table(:,5)==ulvals(ac)),1:3),length(find(albl==ulvals(ac))),1)./255;
end


%% PLOT ELECTRODES (optional)
VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
coord='yes';
if universal_yes(coord), 
    display('...Getting electrode coordinates. Taking *.PIALVOX from elec_recon folder.');
    f=dir([fs_dir '/' subj '/elec_recon/*.PIALVOX']);
    files4bothsides=0;
    if length(f)>1,
        if side=='l',
            f=dir([fs_dir '/' subj '/elec_recon/*left.PIALVOX']);
            files4bothsides=1;
        else
            f=dir([fs_dir '/' subj '/elec_recon/*right.PIALVOX']);
            files4bothsides=1;
        end
    end
    if length(f)>1,
        error('To many possible PIALVOX files.  I do not know which to use.');
    end
    coord=[fs_dir '/' subj '/elec_recon/' f(1).name];
    
    if files4bothsides,
        if side=='l',
            temp=dir([fs_dir '/' subj '/elec_recon/*left.electrodeNames']);
        else
            temp=dir([fs_dir '/' subj '/elec_recon/*right.electrodeNames']);
        end
    else
        temp=dir([fs_dir '/' subj '/elec_recon/*.electrodeNames']);
    end
    if length(temp)>1,
        error('To many possible electrodeNames files.  I do not know which to use.');
    end
    enames_filename=[fs_dir '/' subj '/elec_recon/' temp(1).name];
    fprintf('Attempting to read electrode names from file %s\n', ...
        enames_filename);
    fid=fopen(enames_filename,'r');
    elec_names=textscan(fid,'%s');
    fclose(fid);
    elec_names=elec_names{1};
    elec_names=format_elec_names(elec_names);
end

if exist(coord,'file') && findstr(coord,'PIALVOX')
    VOX_coor=dlmread(coord);
    VOX_coor=VOX_coor(:,2:4);
    ispialvox=1;
elseif isnumeric(coord)
    if size(coord,2)==3
        display(['...Electrode input is matrix with coordinates.']);
        VOX_coor=coord;
    else
        error('...Electrode input is numeric but doesn''t have 3 coordinates');
    end
    ispialvox=0;
end
RAS_coor=(VOX2RAS*[VOX_coor'; ones(1, size(VOX_coor,1))])';

% exclude depth electrodes
if ispialvox==1
    c=1; depth_ind=[];
    for i=1:length(elec_names)
        if ~isempty(strfind(lower(elec_names{i}),'depth')) || strcmpi(elec_names{i}(1),'d')
            depth_ind(c)=i;
            c=c+1;
        end
    end
    if ~isempty(depth_ind),
        fprintf('%d depth electrodes removed.\n',length(depth_ind));
        RAS_coor(depth_ind,:)=[];
        temp_elecnames=cell(1,1);
        c=0;
        for i=1:length(elec_names),
            if ~ismember(i,depth_ind),
                c=c+1;
                temp_elecnames{c}=elec_names{i};
            end
        end
        elec_names=temp_elecnames;
        clear temp_elecnames;
    end
end
RAS_aft=RAS_coor(:,1:3);

n_vert=size(cort.vert,1);
elec_assign=cell(size(RAS_aft,1),2);
for j = 1:size(RAS_aft,1)
        abs_dist=sum(abs(repmat(RAS_aft(j,:),n_vert,1)-cort.vert),2);
        [min_dist, min_vert_id]=min(abs_dist);
        elec_rgb=fvcdat(min_vert_id,:)*255;
        elec_color_id=elec_rgb*[1 2^8 2^16]'; 
        table_id=find(actbl.table(:,5)==elec_color_id);
        if isempty(table_id),
            fprintf('Could not assign electrode %s to a cortical area, rgb: ',elec_names{j});
            disp(elec_rgb);
            fprintf('\n');
            elec_assign{j,2}='Not Cerebral Cortex';
        else
            elec_assign{j,2}=neat_labels{table_id};
            fprintf('%s: %s\n',elec_names{j},elec_assign{j,2});
        end
        elec_assign{j,1}=elec_names{j};
end


%%%% END OF MAIN FUNCTION %%%%


%%%% HELPER FUNCTIONS %%%%
function Y = get_loc_snap_mgh(electrodes,cortex,side,surfacefolder)
% by adykstra

%read smoothed surfaces into matlab
[cortex.lh.vert, cortex.lh.tri] = read_surf([surfacefolder '/lh.pial-outer-smoothed']);
[cortex.rh.vert, cortex.rh.tri] = read_surf([surfacefolder '/rh.pial-outer-smoothed']);

if side == 'l'
    vert_cart = cortex.lh.vert;
elseif side == 'r'
    vert_cart = cortex.rh.vert;
end

%loop over electrodes
for i = 1:max(size(electrodes))
    
    b_x = (vert_cart(:,1)-electrodes(i,1)).^2;
    b_y = (vert_cart(:,2)-electrodes(i,2)).^2;
    b_z = (vert_cart(:,3)-electrodes(i,3)).^2;
    
    temp = sqrt(b_x + b_y + b_z);
    [temp, index] = min(temp);
    clear temp*
    
    %Cartesian Coordiantes, assign electrode location to closest vertex
    x_new(i) = vert_cart(index,1);
    y_new(i) = vert_cart(index,2);
    z_new(i) = vert_cart(index,3);
end
Y = [x_new' y_new' z_new'];

