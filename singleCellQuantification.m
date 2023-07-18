clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(scriptPath));
%all data folders are in a subfolder called 'data' of the directory containing the script
baseDir = fullfile(scriptPath,'data');

%% setup

dirs = {...
    '230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17',...
    '230111_6i9rd_exp20_RD2_GATA3_OTX2_LEF1',...
    };
dataDirs = fullfile(baseDir,dirs);
nrounds = length(dataDirs);
r1 = 1; %first round; should be 1 in general

bare = 'stitched_p%.4d_w%.4d_t%.4d';
channelLabel = cell(1,nrounds);
nucChannels = NaN(1,nrounds);

%load metadata from each round
metas = cell(nrounds,1);
for ri = 1:nrounds
    meta = load(fullfile(dataDirs{ri},'meta.mat'),'meta');
    metas{ri} = meta.meta;
    channelLabel{ri} = metas{ri}.channelLabel;
    nucChannels(ri) = metas{ri}.nucChannel;
end
npos = metas{r1}.nPositions;
nucChannel = nucChannels(r1);

segDir = dataDirs{r1};

%load nuclear and cytoplasmic masks
masks = cell(1,npos); cellData = cell(1,npos); bgmasks = cell(1,npos);
for ii = 1:npos
    prefix = sprintf(bare,ii-1,nucChannel,0);
    segname = [prefix,'_masks.mat'];
    mask = load(fullfile(segDir,segname));
    masks{ii} = mask.masks; cellData{ii} = mask.cellData; bgmasks{ii} = mask.bgmask;
end


%% iterate over colonies, folders, and channels, and read out intensities

positions(npos) = Position();
nucLevels = cell(npos,nrounds); cytLevels = cell(npos,nrounds);
NCratios = cell(npos,nrounds); bgs = cell(npos,nrounds);
for ii = 1:npos    
    fprintf('colony %d of %d\n',ii,npos)
    positions(ii) = Position(metas{r1},ii);
    positions(ii).filenameFormat = [bare,'.tif'];
    positions(ii).cellData = cellData{ii};
    ncells = length(masks{ii});
    
    for ri = 1:nrounds
        fprintf('round %d of %d\n',ri,nrounds)
        nchannels = metas{ri}.nChannels;
        nucLevel = NaN(ncells,nchannels); cytLevel = NaN(ncells,nchannels);
        NCratio = NaN(ncells,nchannels); BG = NaN(1,nchannels);
        for ci = 1:nchannels
            img = positions(ii).loadImage(dataDirs{ri},ci-1,1);
            [nL,cL,ncR,bg] = readIntensityValues(img,masks{ii},bgmasks{ii});
            nucLevel(:,ci) = nL; cytLevel(:,ci) = cL; NCratio(:,ci) = ncR;
            BG(ci) = bg;
        end
        nucLevels{ii,ri} = nucLevel; cytLevels{ii,ri} = cytLevel;
        NCratios{ii,ri} = NCratio; bgs{ii,ri} = BG;
    end

    positions(ii).cellData.nucLevel = cell2mat(nucLevels(ii,:));
    positions(ii).cellData.cytLevel = cell2mat(cytLevels(ii,:));
    positions(ii).cellData.NCratio = cell2mat(NCratios(ii,:));
    positions(ii).cellData.background = cell2mat(bgs(ii,:));
end

meta = metas{r1};
meta.channelLabel = cat(2,channelLabel{:});

save(fullfile(baseDir,'positions.mat'),'positions')
save(fullfile(baseDir,'meta_combined.mat'),'meta')

%% load data if processing has already been done
load(fullfile(baseDir,'positions.mat'),'positions')
load(fullfile(baseDir,'meta_combined.mat'),'meta')

%% collect colony data and export data to a csv file

%radius of the colony in microns
radiusMicron = 350;

%colony indices for which to save data
pidxs = 1:npos;
np = length(pidxs);

%all channel names
channelLabels = meta.channelLabel;
channelLabels = renameDuplicateChannels(channelLabels);
%choose channels for which to save data
chans = 1:length(channelLabels);
nchan = length(chans);
cls = channelLabels(chans);

%find colony centers
cms = NaN(npos,2);
for ii = 1:npos
    cm = setCenter(positions, ii, radiusMicron, meta.xres);
    cms(ii,:) = cm;
end

%filename to which to save the csv
I = regexp(dirs{1},'_RD[1234567890]+_');
savename = [dirs{1}(1:I),'datamatrix.csv'];

excluded = {'XY','nucLevel','cytLevel','NCratio','background'};

fields = fieldnames(positions(1).cellData);
for fi = 1:length(excluded)
    fields = fields(~strcmp(fields,excluded{fi}));
end
nf = length(fields);

varnames = [fields',{'Colony','X','Y','CenterX','CenterY','RadialDist'}];

%collect data on nucleus geometry, cell position, etc.
geometricdata = cell(np,1);
for ii = 1:np
    pidx = pidxs(ii);
    ncell = size(positions(pidx).cellData.XY,1);
    A = NaN(ncell,nf+6);
    for fi = 1:nf
        A(:,fi) = positions(pidx).cellData.(fields{fi});
    end
    xys = positions(pidx).cellData.XY;
    d = meta.xres*pdist2(xys,cms(pidx,:));
    A(:,nf+1) = pidx;
    A(:,nf+2:nf+3) = xys;
    A(:,nf+4) = cms(pidx,1); A(:,nf+5) = cms(pidx,2);
    A(:,nf+6) = d;
    
    geometricdata{ii} = A;
end
geometricdata = cell2mat(geometricdata);

%nuclear, cytoplasmic, and N:C ratio readouts in each channel
fields = {'nucLevel','cytLevel','NCratio'};
prefs = {'','cyto','nc'};
intnames = [];
for fi = 1:length(fields)
    intnames = [intnames, strcat(prefs{fi},cls)]; %#ok<AGROW>
end

ncc = nchan*length(prefs);
intdata = cell(np,1);
for ii = 1:np
    pidx = pidxs(ii);
    ncell = size(positions(pidx).cellData.XY,1);
    A = NaN(ncell,ncc);
    for fi = 1:length(fields)
        A(:,(fi-1)*nchan + 1:fi*nchan) = positions(pidx).cellData.(fields{fi})(:,chans);
    end
    intdata{ii} = A;
end
intdata = cell2mat(intdata);

%make the collected data matrices into a table
A = array2table([intdata,geometricdata],'VariableNames', [intnames,varnames]);

%write the table to a csv file
writetable(A,fullfile(baseDir,savename))



%% local functions

function cm = setCenter(positions, pidx, radiusMicron, xres, margin)
% setCenter()
% setCenter(margin) % margin in microns

if ~exist('margin','var')
    margin = 20;
end

radiusPixel = radiusMicron/xres;

ntime = positions(pidx).nTime;

cm = zeros(ntime,2);
for ti = 1:ntime
    
    % extractData setting center based on mean cell centroid
    areas = positions(pidx).cellData(ti).nucArea;
    xys = positions(pidx).cellData(ti).XY;
    CM = mean(areas.*xys)/mean(areas);
    cm(ti,:) = CM;
    
    % exclude cells/junk outside colony
    d = sqrt((xys(:,1) - cm(ti,1)).^2 ...
            + (xys(:,2) - cm(ti,2)).^2);
        
    outside = d > radiusPixel + margin/xres;
    xys = xys(~outside,:);
    areas = areas(~outside);

    % recenter after removing junk outside
    CM = mean(areas.*xys)/mean(areas);
    cm(ti,:) = CM;
end

end




