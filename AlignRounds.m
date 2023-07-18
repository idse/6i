clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(scriptPath));
%all data folders are in a subfolder called 'data' of the directory containing the script
baseDir = fullfile(scriptPath,'data');

%% setup

%list folders containing image data
dirs = {...
    '230110_6i9rd_exp20_RD1_SMAD23_pAKT_SOX17',...
    '230111_6i9rd_exp20_RD2_GATA3_OTX2_LEF1',...
    };
dataDirs = fullfile(baseDir,dirs);
nrounds = length(dataDirs);

bare = 'stitched_p%.4d_w%.4d_t%.4d';
mipbare = 'stitched_MIP_p%.4d_w%.4d_t%.4d';

%load metadata from each round
channelLabel = cell(1,nrounds);
nucChannels = NaN(1,nrounds);
metas = cell(nrounds,1);
for ri = 1:nrounds
    meta = load(fullfile(dataDirs{ri},'meta.mat'),'meta');
    metas{ri} = meta.meta;
    channelLabel{ri} = metas{ri}.channelLabel;
    nucChannels(ri) = metas{ri}.nucChannel;
end
npos = metas{1}.nPositions;

%get file extension for mips
listing = dir(fullfile(dataDirs{1},'MIP','stitched_MIP_*'));
[~,~,mipext] = fileparts(fullfile(listing(1).folder,listing(1).name));


%% do alignment
%directories for which to do alignment
%each new folder can be aligned as it is acquired instead of all at once
rounds = 2:nrounds;

%iterate over positions
for pidx = 1:npos
    fprintf('position %d of %d\n',pidx,npos)
    %load DAPI image for the first channel
    fname = fullfile(dataDirs{1},[sprintf(bare,pidx-1,nucChannels(1),0),'.tif']);
    DAPI1 = loadTiffStack(fname);
    nz1 = size(DAPI1,3);
    mip1 = max(DAPI1,[],3);
    lim1 = stitchedlim(mip1);
    
    for ri = rounds
        fprintf('round %d of %d\n',ri,nrounds)
        QCdir = fullfile(dataDirs{ri},'overlays');
        if ~exist(QCdir,'dir'), mkdir(QCdir); end

        prefix = sprintf(bare,pidx-1,nucChannels(ri),0);
        fname = fullfile(dataDirs{ri},[prefix,'.tif']);
        DAPI2 = loadTiffStack(fname);
        mip2 = max(DAPI2,[],3);
        lim2 = stitchedlim(mip2);

        %xy alignment
        shiftyx = findImageShift(mip1,mip2,'automatic');
        mipa = alignImage(mip1,mip2,shiftyx);

        %save mip overlay
        overlay = makeMIPOverlay(mip1,mipa,lim1,lim2,0.075);
        savename = fullfile(QCdir,[prefix,'_MIPoverlay.png']);
        imwrite(overlay,savename);
        
        %align the entire DAPI stack in xy
        dapia = xyalignImageStack(DAPI2,shiftyx);

        %z alignment
        [shiftz, scalez] = zAlignImageStacks(DAPI1,dapia);
        %apply the z shift and scaling
        dapia = zShiftScale(dapia,shiftz,scalez,nz1);
        
        %save cross section overlay
        index = round(size(dapia,1)/2);
        crossSectionOverlay(DAPI1,dapia,metas{1},index,lim1,lim2,0.3)
        savename = fullfile(QCdir,[prefix,sprintf('_crossSectionX%d_overlay.png',index)]);
        saveas(gcf,savename)
        exportgraphics(gcf, savename,'Resolution',300)

        %iterate over channels, save aligned images
        nchannels = metas{ri}.nChannels;
        for ci = 1:nchannels
            fname = fullfile(dataDirs{ri},[sprintf(bare,pidx-1,ci-1,0),'.tif']);
            
            if ci - 1 == nucChannels(ri)
                img = dapia;
            else
                %yxz image stack in channel ci
                img = loadTiffStack(fname);
                %apply shift and warping to the image stack
                img = xyalignImageStack(img,shiftyx);
                %apply z shift + scaling
                img = zShiftScale(img,shiftz,scalez,nz1);
            end
            
            %write the aligned image stack
            writeTiffZStack(img, fname)
            %write the aligned MIP
            MIP = max(img,[],3);
            mipname = fullfile(dataDirs{ri},'MIP',[sprintf(mipbare,pidx-1,ci-1,0),mipext]);
            if strcmp(mipext,'.jpg')
                imwrite(im2double(MIP),mipname,'Quality',99)
            else
                imwrite(MIP,mipname)
            end
            
        end
    end
    close all
end


%% local functions

function IMa = xyalignImageStack(img,shiftyx)
m = size(img,1); n = size(img,2); nz = size(img,3); imclass = class(img);

IMa = zeros(m,n,nz,imclass);
for zi = 1:nz
    im = img(:,:,zi);
    ima = alignImage(im,im,shiftyx);    
    IMa(:,:,zi) = ima;
end

end

function imz = zShiftScale(img,shift,scale,nz1)
m = size(img,1); n = size(img,2); nz = size(img,3);
nznew = round(nz*scale);

[XX,YY,ZZ] = meshgrid(1:n,1:m,1:nz);
[Xnew,Ynew,Znew] = meshgrid(1:n,1:m,linspace(1 - shift,nz - shift,nznew));

imz = uint16(interp3(XX,YY,ZZ,single(img),Xnew,Ynew,Znew));

nznew = size(imz,3);
if nznew > nz1
    imz = imz(:,:,1:nz1);
elseif nznew < nz1
    imz = cat(3,imz,zeros(m,n,nz1-nznew,'uint16'));
end

end

function RGB = makeMIPOverlay(mip1,mip2,lim1,lim2,cfs)
m = size(mip1,1); n = size(mip1,2);

A = imadjust(mip1,lim1); B = imadjust(mip2,lim2);
overlay = cat(3,A,B,A);

%use arial font, but different matlab versions have different names for
%available fonts
fonts = listTrueTypeFonts;
font = 'Arial';
if sum(strcmp(font,fonts)) == 0
    font = 'Arial Unicode';
    if sum(strcmp(font,fonts)) == 0
        I = find(contains(fonts,'Arial'),1);
        font = fonts{I};
    end
end

%add label text
text = {'round 1', 'round 2'};
ypos = round([m*(1 - 0.75*cfs); m*(1 - 0.75*cfs)]);
xpos = round([0.005*n; 0.15*n]);
RGB = insertText(overlay,[xpos,ypos],text,'BoxColor','black',...
    'TextColor',{'m','g'},'BoxOpacity',0.3,'FontSize',96,...
    'Font',font);

end

function crossSectionOverlay(img1,img2,meta1,index,lim1,lim2,cfs)

xres = meta1.xres;
zres = meta1.zres;

cross1 = transpose(squeeze(img1(index,:,end:-1:1)));
zsize1 = round(size(cross1,1)*zres/xres);
cross1 = imresize(cross1,[zsize1,size(cross1,2)]);

cross2 = transpose(squeeze(img2(index,:,end:-1:1)));
zsize2 = round(size(cross2,1)*zres/xres);
cross2 = imresize(cross2,[zsize2,size(cross2,2)]);

if zsize1 > zsize2
    cross2 = [zeros(zsize1-zsize2,size(cross2,2),'uint16'); cross2];
elseif zsize2 > zsize1
    cross1 = [zeros(zsize2-zsize1,size(cross1,2),'uint16'); cross1];
end
zsize = size(cross1,1);
n = size(cross1,2);

C1 = imadjust(cross1,lim1); C2 = imadjust(cross2,lim2);
o1 = cat(3,C1,C2,C1);

titles = {'round 1','round 2',...
    "\color{magenta}round 1 \color{green} round 2"};
crosses = {C1,C2,o1};

figpos = figurePosition([n,2*zsize*length(crosses)]);
figure('Position',figpos)
for ii = 1:length(crosses)
    subplot_tight(length(crosses),1,ii)
    imshow(crosses{ii})
    cleanSubplot
    
    ypos = zsize*0.5*cfs; xpos = 0.01*n*cfs;
    text(xpos, ypos, titles{ii},...
        'Color','w','FontUnits','normalized','FontSize',cfs,...
        'FontWeight','bold')
end


end











