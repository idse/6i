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

bare = 'stitched_p%.4d_w%.4d_t%.4d';
mipbare = 'stitched_MIP_p%.4d_w%.4d_t%.4d.jpg';
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
npos = metas{1}.nPositions;
meta = metas{1};

channelLabels = cat(2,channelLabel{:});

nchannel = length(channelLabels);
dirids = cumsum(strcmp(channelLabels,'DAPI'));
ndir = length(dataDirs);
cdi = NaN(1,nchannel);
for ri = 1:ndir
    cdi(dirids == ri) = 1:sum(dirids == ri);
end

channelLabels = renameDuplicateChannels(channelLabels);


%% choose channels for which to make overlays, load and overlay the images

%pick a colony/position
pidx = 1;
%pick channels for which to load images, in rgb order for display
channels = {'GATA3','LEF1','OTX2'};

%if the same channel is specified more than once, only load it once
[uchan, ~, ic] = unique(channels);
chans = getChannelInds(uchan,channelLabels);
nc = length(chans);

%load the images
ims = cell(1,1,nc);
for cii = 1:nc
    ci = chans(cii);
    fprintf("channel %d\n",cii)
    disp(channelLabels{ci})
    di = dirids(ci);
    disp(strcat(sprintf("channel %d in ",cdi(ci)), dirs{di}))
    fname = fullfile(dataDirs{di},'MIP',sprintf(mipbare,pidx-1,cdi(ci)-1,0));
    img = imread(fname);
    
    lim = seglim(img);
    ims{cii} = imadjust(img,lim);
end

rgb = cell2mat(ims(ic));

%display the combined color image
figure
imshow(rgb)
cleanSubplot

%add a label for which channels are in which colors
addColorLabel(channels,size(rgb),0.065)

%% cross section overlay
%this requires loading the full z stack in each channel so is generally
%much slower

%pick a colony/position
pidx = 1;
%pick channels for which to load images, in rgb order for display
channels = {'OTX2','LEF1','DAPI_1'};
[uchan, ~, ic] = unique(channels);
chans = getChannelInds(uchan,channelLabels);
nc = length(chans);

%choose whether to do a cross section in x or y and at what point in the
%image; if index is left empty ([]), use the middle of the image
crossmode = 'x';
index = [];

zres = meta.zres; xres = meta.xres;

ims = cell(1,1,nc);
for cii = 1:nc
    ci = chans(cii);
    fprintf("channel %d\n",cii)
    disp(channelLabels{ci})
    di = dirids(ci);
    disp(strcat(sprintf("channel %d in ",cdi(ci)), dirs{di}))
    fname = fullfile(dataDirs{di},[sprintf(bare,pidx-1,cdi(ci)-1,0),'.tif']);
    img = loadTiffStack(fname);
    %determine contrast limits using the mip
    mip = max(img,[],3);    
    lim = seglim(mip);

    if strcmp(crossmode,'x')
        if isempty(index)
            index = round(size(mip,2)/2);
        end
        im = transpose(squeeze(img(:,index,end:-1:1)));
    elseif strcmp(crossmode,'y')
        if isempty(index)
            index = round(size(mip,1)/2);
        end
        im = transpose(squeeze(img(index,:,end:-1:1)));
    end
    zsize = round(size(im,1)*zres/xres);
    im = imresize(im,[zsize,size(im,2)]);

    ims{cii} = imadjust(im,lim);
end

rgb = cell2mat(ims(ic));

figure
imshow(rgb,'Border','tight')
cleanSubplot

%add a label for which channels are in which colors
addColorLabel(channels,size(rgb),0.1);


%% local functions
%don't run this block

function chans = getChannelInds(channels,channelLabels)

chans = NaN(size(channels));
for ii = 1:length(chans)
    I = find(strcmp(channels{ii},channelLabels));
    if isempty(I)
        I = find(strcmp([channels{ii},'_1'],channelLabels));
        if isempty(I)
            error(strcat("could not find the channel ",channels{ii}))
        else
            chans(ii) = I;
        end
    else
        chans(ii) = I;
    end
end


end

function clabel = addColorLabel(channels,imsize,cfs)

m = imsize(1); n = imsize(2);

colors = {'red','green','blue'};
clabel = '';
for ii = 1:length(channels)
    clabel = strcat(clabel,"\color{",colors{ii},"}",channels{ii}," ");
end

ypos = m*(1-0.5*cfs); xpos = 0.01*n;
text(xpos, ypos, clabel,...
    'FontUnits','normalized','FontSize',cfs,'FontWeight','bold')

end








