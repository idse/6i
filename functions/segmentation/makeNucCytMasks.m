function masks = makeNucCytMasks(nucseg,fgseg,opts)

if nargin < 3
    opts = struct();
end

if ~isfield(opts, 'nucShrinkage')
    opts.nucShrinkage = 0;
end
if ~isfield(opts, 'cytoSize')
    opts.cytoSize = 5;
end
if ~isfield(opts, 'cytoMargin')
    opts.cytoMargin = 0;
end
if ~isfield(opts, 'bgMargin')
    opts.bgMargin = 2;
end
if isfield(opts,'junkmask')
    junkmask = opts.junkmask;
else
    junkmask = [];
end
if ~isfield(opts,'cytoplasmicLevels')
    opts.cytoplasmicLevels = false;
end

if ~isfield(opts,'suppressOutput')
    suppressOutput = false;
else
    suppressOutput = opts.suppressOutput;
end

if isempty(fgseg)
    if ~suppressOutput
        warning('using inverse dilated nuclear mask as background')
    end
    nz = 1;
    bgmask = ~imdilate(nucseg,strel('disk',2*opts.cytoSize));
    fgmask = [];
else
    nz = size(fgseg,3);
    % get background mask
    fgsegmip = sum(fgseg,3) > 0;
    %could alternatively also have a different bg mask for each z slice,
    %but common bgmask excluding every pixel that is foreground in at least
    %one z slice seems reasonable
    bgmask = imerode(~imclose(fgsegmip,strel('disk',5)),strel('disk', opts.bgMargin)) > 0;
    
    fgmask = zeros(size(fgseg),class(fgseg));
    for zi = 1:nz
        mask = imclose(fgseg(:,:,zi),strel('disk',5));
        fgmask(:,:,zi) = imerode(mask,strel('disk', opts.bgMargin));
    end
end

%exclude junk from the background mask
if ~isempty(junkmask)
    bgmask = bgmask & ~(sum(junkmask,3)>0);
end

nucmaskraw = nucseg;
nucmask = nucseg;

% make cytoplasmic mask
%----------------------
if opts.cytoplasmicLevels

    if opts.cytoMargin > 0
        nucmaskmarg = imdilate(nucmask,strel('disk',opts.cytoMargin));
    else
        nucmaskmarg = imerode(nucmask,strel('disk',-opts.cytoMargin));
    end

    % watershedding inside the dilation
    dilated = imdilate(nucmask, strel('disk',opts.cytoSize + opts.cytoMargin));
    basin = bwdist(dilated);
    basin = imimposemin(basin, nucmask);
    L = watershed(basin);

    stats = regionprops(L, 'PixelIdxList');
    npidxlist = cat(1,{stats.PixelIdxList});
    ncells = length(npidxlist);
    
    cytmask = repmat(npidxlist,nz,1);
    % replacement of exclusions above
    for cci = 1:ncells
        %starting cyt mask is the same for all z slices
        CCPIL = cytmask{1,cci};
        %remove things that are the same across z slices first
        CCPIL = CCPIL(dilated(CCPIL));    % exclude outside dilated nuclei
        if opts.cytoMargin ~= 0
            CCPIL = CCPIL(~nucmaskmarg(CCPIL)); % exclude margin around nuclei
        else
            CCPIL = CCPIL(~nucmaskraw(CCPIL));% exclude nuclei before cleanup
            CCPIL = CCPIL(~nucmask(CCPIL));% exclude nuclei 
        end
        %exclude junk if junkmask is segmented on the MIP
        if ~isempty(junkmask) && size(junkmask,3) == 1
            CCPIL = CCPIL(~junkmask(CCPIL));
        end
        %for each z slice, only keep the part of the cyt mask that is
        %segmented as foreground in that slice
        for zi = 1:nz
            if ~isempty(fgmask)
                mask = fgmask(:,:,zi);
                ccpil = CCPIL(mask(CCPIL));% exclude background
            else
                ccpil = CCPIL;
            end
            %if junkmask is segmented across z slices, exclude junk here
            if ~isempty(junkmask) && size(junkmask,3) == nz
                mask = junkmask(:,:,zi);
                ccpil = ccpil(mask(ccpil));
            end
            
            cytmask{zi,cci} = ccpil;
        end
        
    end
else
    cytmask = {};
end

% regionprops and indices from nuclear mask
% follows cyt so we can shrink without introducing new var
%-------------------------------------------

% shrink instead of erode to preserve number of connected
% components to not introduce mismatch between nuclear and
% cytoplasmic mask
if opts.nucShrinkage > 0
    nucmask = bwmorph(nucmask,'shrink',opts.nucShrinkage);
end

% initialize cellData
%---------------------------------------------------
% props = 'Area', 'Centroid','PixelIdxList','Orientation', 'MajorAxisLength','MinorAxisLength','Circularity';
stats = regionprops(nucmask,'PixelIdxList','Centroid',...
    'Orientation','MajorAxisLength','MinorAxisLength','Circularity');
nucmask = {stats.PixelIdxList};
XY = cell2mat({stats.Centroid}');
ncells = length(stats);

if opts.cytoplasmicLevels
    %shuffle cytoplasmic masks to be in the same order as
    %nuclear masks
    cytPixelLists = cell(nz,ncells);
    cytidxs = (1:ncells)'; nucidxs = NaN(ncells,1);
    for cellidx = 1:ncells
        nucidxs(cellidx) = unique(L(stats(cellidx).PixelIdxList));
        ind = unique(L(cytmask{1,cellidx}));
        if ~isempty(ind)
            if ind ~= cellidx
                error('cytCC indexing is wrong')
            end
        end

        cytidx = cytidxs == nucidxs(cellidx);
        for zi = 1:nz
            cytPixelLists{zi,cellidx} = cytmask{zi,cytidx};
        end
    end
    cytmask = cytPixelLists;
end


stats = rmfield(stats,{'PixelIdxList','Centroid'});

%rename fields to match what we usually use in cellData
sfields = {'Orientation','MajorAxisLength','MinorAxisLength','Circularity'}; %source field names
tfields = {'nucOrientation','nucMajorAxis','nucMinorAxis','nucCircularity'}; %target field names
if ncells > 0
    nucstats(ncells) = struct;
    for ii = 1:ncells
        for fi = 1:length(sfields)
            nucstats(ii).(tfields{fi}) = stats(ii).(sfields{fi});
        end
    end
else
    test = struct();
    for fi = 1:length(tfields)
        test.(tfields{fi}) = [];
    end
    nucstats = test([],1);
end

masks = struct(...
    'nucmask',      {nucmask},...
    'XY',           XY,...
    'cytmask',      {cytmask},...
    'bgmask',       bgmask,...
    'fgmask',       fgmask,...
    'nucmaskraw',   nucmaskraw,...
    'nucstats',     nucstats);



end