clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(scriptPath));

%% setup and options
%specify data folder that was used for segmentation (first round)
dataDir = fullfile(scriptPath,'data','230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17');

%type of input segmentation -'ilastik', 'cellpose', or 'combined'
%minor modifications should allow loading arbitrary segmentations from the
%user's preferred software
segmode = 'combined';

%options for nuclear and cytoplasmic masks
%size filtering thresholds for nuclear masks during consolidation
minarea = 200; %area threshold for getting rid of junk
maxarea = 3000; %area threshold for discarding fused clumps of nuclei

opts = struct(...
    'cytoSize',             4,... %width of the cytoplasmic annulus in pixels
    'cytoMargin',           2,... %margin from the nucleus to the cytoplasmic annulus
    'nucShrinkage',         0,... %number of pixels by which to shrink the nuclear mask
    'cytoplasmicLevels',    true); %boolean for whether to read out cytoplasmic levels

%options for linking 2D masks over z slices to make 3D masks
zopts = struct(...
    'IoU',                  0.5,... %intersection over union threshold for linking in z
    'maxZsize',             15,... %mazimum height for a nucleus in microns
    'minZsize',             3,... %minimum height for a nucleus in microns
    'maxCentroidDistance',  5); %max distance between centroids of adjacent nuclei

bare = 'stitched_p%.4d_w%.4d_t%.4d';
load(fullfile(dataDir,'meta.mat'),'meta')
nucChannel = meta.nucChannel;
npos = meta.nPositions;

%% combine ilastik, cellpose
overlaydir = fullfile(dataDir, 'SegOverlays');
if ~exist(overlaydir,'dir'), mkdir(overlaydir); end

%consolidate cpmasks from separate pngs into multipage tifs
if strcmp(segmode,'cellpose') || strcmp(segmode,'combined')
    consolidate_cp_masks(dataDir,npos,bare,nucChannel)
end

tic
for ii = 1:npos
    fprintf('position %d of %d\n',ii,npos)
    fname = sprintf(bare,ii-1,nucChannel,0);
    fname = fullfile(dataDir,fname);
    imname = [fname,'.tif'];
    t = imfinfo(imname);
    nz = length(t); m = t(1).Height; n = t(1).Width;

    cpname = [fname,'_cp_masks.tif'];
    ilastikname = [fname '_Simple Segmentation.h5'];

    if strcmp(segmode,'ilastik') || strcmp(segmode,'combined')
        ilastikseg = ilastikRead(ilastikname);
    else
        ilastikseg = false(m,n,nz);
    end

    sname = [fname,'_FinalSegmentation.tif'];
    for zi = 1:nz
        fprintf('.')
        if zi == 1
            mode = 'overwrite';
        else
            mode = 'append';
        end
        img = imread(imname,zi);
        I = imadjust(img,stitchedlim(img,0.01));
        
        if strcmp(segmode,'cellpose') || strcmp(segmode,'combined')
            cpmask = imread(cpname,zi);
        else
            cpmask = false(m,n);
        end
        
        ilastikmask = ilastikseg(:,:,zi);
        %difference of the masks
        mask = ilastikmask & ~(cpmask > 0);
        seg = imopen(mask,strel('disk',4));
        % combine the two masks
        L = uint16(labelmatrix(bwconncomp(seg)));
        L(L>0) = L(L>0) + max(cpmask,[],'all');

        finalseg = bitor(cpmask,L);
        finalseg = finalseg > 0 & ~boundarymask(finalseg);
        props = regionprops(finalseg,'PixelIdxList','Area');
        mask = [props.Area] < minarea | [props.Area] > maxarea;
        idxs = cell2mat({props(mask).PixelIdxList}');
        finalseg(idxs) = 0;

        finalseg = uint8(finalseg);
        
        imwrite(finalseg,sname,'WriteMode',mode)
        
        [~,name,~] = fileparts(fname);
        overlay = visualize_nuclei_v2(finalseg>0,I);
        savename = [name, sprintf('_z%.4d',zi-1), '_SegOverlay.jpg'];
        savename = fullfile(overlaydir,savename);
        imwrite(overlay,savename);
    end
    fprintf('\n')
end
toc

%% define nuclear and cytoplasmic masks for each cell and link across z slices

for pidx = 1:npos
    prefix = sprintf(bare,pidx-1,nucChannel,0);
    disp(prefix)
    segname = fullfile(dataDir,[prefix,'_FinalSegmentation.tif']);
    t = imfinfo(segname);
    nz = length(t);
    m = t(1).Height; n = t(1).Width;
    
    disp('reading segmentations')
    tic
    seg = false(m,n,nz);
    for zi = 1:nz
        fprintf('.')
        seg(:,:,zi) = imread(segname,zi) > 0;
    end
    fprintf('\n')
    toc
    
    % make nuclear and cytoplasmic masks for each cell in each frame
    disp('making nuclear and cytoplasmic masks')
    tic
    clear allmasks
    allmasks(nz) = struct; %#ok<SAGROW>
    for zi = 1:nz
        if zi == 1
            opts.suppressOutput = false;
        else
            opts.suppressOutput = true;
        end
        masks = makeNucCytMasks(seg(:,:,zi),[],opts);
        fields = fieldnames(masks);
        for fi = 1:length(fields)
            allmasks(zi).(fields{fi}) = masks.(fields{fi});
        end
        fprintf('.')
    end
    fprintf('\n')
    toc
    
    disp('consolidating in z')
    tic
    [masks, cellData, chain] = linkMasksInZ(allmasks, meta, zopts);
    toc

    bgmask = cell2mat(reshape({allmasks.bgmask},[1,1,length(allmasks)]));

    save(fullfile(dataDir,[prefix,'_masks.mat']),'masks','bgmask','cellData')
end


%% local functions
%don't run this block

function consolidate_cp_masks(dataDir,npos,bare,nucChannel)

disp('consolidating cellpose masks')

outputdir = fullfile(dataDir, 'cp_masks');
if ~exist(outputdir,'dir'), mkdir(outputdir); end

savebare = [bare,'_cp_masks.tif'];

%if already done, do not reprocess
prefix = sprintf(bare,0,nucChannel,0);
name = [prefix,sprintf('_z%.4d_cp_masks.png',0)];
readname = fullfile(dataDir,name);
writename = fullfile(outputdir,name);

if ~exist(readname,"file") && exist(writename,"file")
    doconsolidate = false;
else
    doconsolidate = true;
end

if doconsolidate
    for ii = 1:npos
        prefix = sprintf(bare,ii-1,nucChannel,0);
        listing = dir(fullfile(dataDir,[prefix,'*cp_masks.png']));
        nz = length(listing);
        disp(prefix)
        fprintf('nz = %d\n',nz)
        savename = fullfile(dataDir,sprintf(savebare,ii-1,nucChannel,0));
        for zi = 1:nz
            fprintf('.')
            if zi == 1
                mode = 'overwrite';
            else
                mode = 'append';
            end
            name = [prefix,sprintf('_z%.4d_cp_masks.png',zi-1)];
            readname = fullfile(dataDir,name);
            img = imread(readname);
            imwrite(img,savename,'WriteMode',mode)
    
            writename = fullfile(outputdir,name);
            movefile(readname,writename)
        end
        fprintf('\n')
    end
end

end




