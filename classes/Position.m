classdef Position < handle
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        ncells              % number of cells

        cellData            % structure array indexed by time:
                            % - XY
                            % - area
                            % - nucLevel:       nCell x nChannels array
                            % - cytLevel
                            % - background:     nChannels long vector

        timeTraces          % structure: reorganizes cellData
                            % - nucLevelAvg:    vector indexed by time
                            % - cytLevelAvg
                            % - background
                            % NOT IMPLEMENTED YET: traces of tracked cells

        dataChannels        % we may not want to quantify all channels

        nChannels           % number of colors imaged
        nucChannel          % nuclear channel
        
        nTime               % number of time points
        filenameFormat
        tree                %stores a tree object with tracks (deprecated)
        G                   %digraph with links from tracking
    end

    properties (SetAccess = protected)
        ID                  % identifier of colony
    end

    properties (Dependent)

        filename
        barefilename;       % filename without extension
        %data                % cell by cell data for colony
                            % data cols: x, y, area, -1,
                            % (nuclear, cytoplasmic) mean for each channel
    end

    methods

        % constructor
        function this = Position(varargin)
            % Position(nChannels, filenameFormat, nTime, nucChannel)
            % Position(meta)
            % Position(meta, id)

            if nargin == 4
                nChannels = varargin{1};
                filenameFormat = varargin{2};
                nTime = varargin{3};
                nucChannel = varargin{4};

            elseif nargin == 1 || nargin == 2
                
                meta = varargin{1};
                nChannels = meta.nChannels;
                filenameFormat = meta.filenameFormat;
                nTime = meta.nTime;
                nucChannel = meta.nucChannel;
                
                if nargin == 2
                    this.ID = varargin{2};
                    conditionIdx = ceil(this.ID/meta.posPerCondition);
                    conditionPos = rem(this.ID-1, meta.posPerCondition)+1;
                    if isempty(filenameFormat) && ~isempty(meta.fileNames{1})
                        filenameFormat = meta.fileNames{conditionIdx}{conditionPos};
                    end
                end

            elseif nargin == 0 % matlab sucks, always call superclass constructor with no args
                return
            else
                error('wrong number of arguments');
            end

            if ~exist('nTime','var')
                nTime = 1;
                this.timeTraces = [];
            else
                this.timeTraces = struct();
            end

            if ~isnumeric(nChannels)
                error('first argument is nChannels, which is a number');
            end
            if ~isnumeric(nTime)
                error('third argument is nTime, which is a number');
            end

            this.nucChannel = nucChannel;
            this.nChannels = nChannels;
            this.nTime = nTime;

            % strip off path in case it was provided
            if ~isempty(filenameFormat)
                [~, name, ext] = fileparts(filenameFormat);
                this.filenameFormat = char(strcat(name,ext));
            end

            this.cellData = struct();
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channels, time)
            % load position image
            %
            % img = loadImage(dataDir, channels)
            % img = loadImage(dataDir, channels, time)
            %
            % dataDir:  main data directory
            % channels: desired channels to be loaded, leave empty for all
            % time:     time to load
            %
            % img:      loaded image

            if ~exist('channels','var')
                channels = 1:this.nChannels;
            end
            if ~exist('time','var')
                time = 1;
            end
            if ~exist('useMIP','var')
                useMIP = false;
            end

            [~,~,ext] = fileparts(this.filenameFormat);

            if strcmp(ext,'.tif') || strcmp(ext,'.btf')

                s = strsplit(this.filenameFormat,'_F[0-9]','DelimiterType','RegularExpression');
                if numel(s) > 1
                    fnameFormat = [s{1} '_MIP_p%.4d_w%.4d.tif'];
                else
                    fnameFormat = this.filenameFormat;
                end

                for cii = 1:numel(channels)

                    if contains(fnameFormat,'_t') && (contains(fnameFormat,'_c') || contains(fnameFormat,'_w'))
                        fname = sprintf(fnameFormat,this.ID - 1, channels(cii), time-1);
                        iminf = imfinfo(fullfile(dataDir,fname));
                        im = zeros([iminf(1).Height iminf(1).Width numel(iminf)],'uint16');
                        for zii = 1:numel(iminf)
                            im(:,:,zii) = imread(fullfile(dataDir, fname), zii);
                        end
                        img{cii} = im;
                        
                    elseif ~contains(fnameFormat,'_t') && (contains(fnameFormat,'_c') || contains(fnameFormat,'_w'))
                        % this doesn't handle z and t in the same tiff
                        fname = sprintf(fnameFormat,this.ID - 1, channels(cii));
                        img{cii} = imread(fullfile(dataDir, fname), time);
                        
                    else
                        % this doesn't handle z and t in the same tiff
                        fname = sprintf(fnameFormat,this.ID - 1);
                        img{cii} = imread(fullfile(dataDir, fname), time);
                    end
                end
                img = cat(3,img{:});

            % bioformats
            elseif strcmp(ext,'.vsi') || strcmp(ext, '.oif') ...
                    || strcmp(ext, '.oib') || strcmp(ext, '.ims') || strcmp(ext, '.nd2')

                fname = fullfile(dataDir, this.filename);

                r = bfGetReader(fname);
                img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ()], 'uint16');
                for cii = 1:numel(channels)
                    for zi = 1:r.getSizeZ()
                        img(:,:,cii,zi) = bfGetPlane(r, r.getIndex(zi-1,channels(cii),time-1)+1);
                    end
                end
                r.close();

                img = squeeze(img);
            else
                error('unknown format');
            end
            if time == 1
                disp(['loaded image ' fname]);
            end
        end

        function seg = loadSegmentation(this, dataDir, channel, segfileformat)
            % load segmentation
            %
            % seg = loadSegmentation(dataDir, channel)
            %
            % dataDir:  main data directory
            %
            % seg:      binary image of nuclei
            %
            % segmentation is expected to have same name as raw data
            % but with .tif -> '_Simple Segmentation.h5'
            %
            % by convention use green for foreground in Ilastik

            if ~exist('channel','var')
                error('please provide channel');
            end

            if exist('segfileformat','var')

                fname = sprintf(segfileformat, channel);
            else

                % guessing filename based on conventions
                %------------------------------------------

                % possible cases for filename vs segmentation:
                % 1) 1.a) Andor multiposition
                %    blabla_F00.ims -> blabla_MIP_p0000_w0000.tif 
                %                   -> blabla_MIP_p0000_w0000_Simple Segmentation.h5
                %    1.b) can use MIP instead of raw data for quantification
                %
                % 2) Snaps with well names for condition
                %    A1_blabla.vsi  -> A1_blabla_MIP_w0000.tif
                %                   -> A1_blabla_MIP_w0000_Simple Segmentation.h5
                
                barefname = {};

                % case 1:
                
                % split by _ <letter> <number>
                s = strsplit(this.filename,'_[Ffwptm][0-9]+','DelimiterType','RegularExpression');
                
                % get the position index from the raw data
                [si,ei] = regexp(this.filename,'_[Ffpm][0-9]+');
                pi = str2num(this.filename(si+2:ei));

                % add MIP if the 'raw' data filename is not a MIP (usual
                % case) because we always segment MIPS
                % 'raw' data is MIP in the case for the stitched grids
                if contains(s{1},'_MIP') % case 1.b)
                    barefname{1} = sprintf([s{1}  '_p%.4d'], pi);
                else  % case 1.a)
                    barefname{1} = sprintf([s{1}  '_MIP_p%.4d'], pi);
                end
                
                % case 2:
                barefname{2} = this.filename(1:end-4);
                
                % find the first barefname for which there is an h5 file
                i = 0;
                listing = [];
                while(isempty(listing))
                    i = i+1;
                    if i > numel(barefname)
                        error(['segmentation for ' this.filename ' for channel ' num2str(channel) ' not found in ' dataDir]);
                    else
                        listing = dir(fullfile(dataDir,[barefname{i} '*.h5']));
                    end
                end

                fname = [];
                for i = 1:numel(listing)
                    % DON'T COMMENT THE IF STATEMENT BELOW, IT WILL LOAD THE WRONG
                    % CHANNEL
                    if ~isempty(strfind(listing(i).name, sprintf('_w%.4d',channel)))...
                            || ~isempty(strfind(listing(i).name,sprintf('_c%d',channel)))

                        fname = fullfile(dataDir,listing(i).name);
                    end
                end
            end

            % actually reading the file
            %----------------------------

            if ~exist(fname,'file')

                warning(['segmentation for channel ' num2str(channel) ' not found in ' dataDir, ', may be naming convention problem']);
                seg = [];

            else
                seg = h5read(fname, '/exported_data');
                
                % Probabilities -> treshold to get binary
                % this is not set up to deal with more than 2 classes yet
                % (fg/bg)
                if size(seg,1) > 1
                    ssegsize = size(seg);
                    ssegsize(1) = 1;
                    simpleseg = false(ssegsize);
                    simpleseg(seg(1,:,:,:) >= 0.5) = true;
                    seg = squeeze(simpleseg);
                    
                % Simple Segmentation (could be more than 2 classes)
                else
                    seg = squeeze(seg);
                end
                % Ilastik output data has xy transposed
                seg = permute(seg, [2 1 3]);
                disp(['loaded segmentation ' fname]);
            end
        end

        function makeAvgTimeTraces(this, SNRcutoff)
            % make time traces of levels in cellData
            %
            % makeTimeTraces()
            % makeTimeTraces(SNRcutoff)
            %
            % populates Position.timeTraces
            %
            % SNRcutoff :   exclude cell if N/bg or C/bg < SNRcutoff,
            %               default 2

            if ~exist('SNRcutoff','var')
                SNRcutoff = zeros([numel(this.dataChannels) 1]);
            end

            nucLevelAvg = zeros([this.nTime numel(this.dataChannels)]);
            nucLevelMed = zeros([this.nTime numel(this.dataChannels)]);

            cytLevelAvg = zeros([this.nTime numel(this.dataChannels)]);
            cytLevelMed = zeros([this.nTime numel(this.dataChannels)]);
            ratioAvg = zeros([this.nTime numel(this.dataChannels)]);

            bg = zeros([this.nTime numel(this.dataChannels)]);

            for ti = 1:numel(this.cellData)

                nucLevel = this.cellData(ti).nucLevel;

                if isfield(this.cellData(ti),'cytLevel')
                    cytLevel = this.cellData(ti).cytLevel;
                else
                    cytLevel = [];
                end

                A = this.cellData(ti).nucArea;

                for ci = 1:numel(this.dataChannels)

                    bg(ti,ci) = this.cellData(ti).background(ci);

                    % implement SNR cutoff necessary for varying expression
                    % levels to exclude outliers from dim cells
                    if ~isempty(cytLevel)
                        tooDim =    cytLevel(:,ci)/bg(ci) < SNRcutoff(ci) |...
                                    nucLevel(:,ci)/bg(ci) < SNRcutoff(ci);
                        idx = ~tooDim & ~isnan(nucLevel(:,ci)) & ~isnan(cytLevel(:,ci));
                    else
                        idx = ~isnan(nucLevel(:,ci));
                    end

                    if sum(idx)/numel(idx) < 0.2
                        warning([num2str(100*(1-sum(idx)/numel(idx)),3) '% of cells excluded based on SNR in channel ' num2str(ci) ' at time ' num2str(ti) '!']);
                    end

                    if ~isempty(idx)

                        % weight by cell area - may or may not be sensible
                        % in context
%                         W = A(idx)/mean(A(idx));
                        W = ones(size(A(idx)));
                        nucLevelAvg(ti, ci) = mean(nucLevel(idx,ci).*W);
                        nucLevelMed(ti, ci) = median(nucLevel(idx,ci).*A(idx))/mean(A(idx));

                        if ~isempty(cytLevel)
                            cytLevelAvg(ti, ci) = mean(cytLevel(idx,ci).*W);
                            cytLevelMed(ti, ci) = median(cytLevel(idx,ci).*A(idx))/mean(A(idx));
                            ratioAvg(ti,ci) = mean((nucLevel(idx,ci) - bg(ti,ci))./(cytLevel(idx,ci) - bg(ti,ci)));
                        end
                    end
                end


            end

            this.timeTraces.nucLevelAvg = nucLevelAvg;
            this.timeTraces.nucLevelMed = nucLevelMed;

            this.timeTraces.cytLevelAvg = cytLevelAvg;
            this.timeTraces.cytLevelMed = cytLevelMed;
            this.timeTraces.ratioAvg = ratioAvg;

            this.timeTraces.background = bg;
        end

        function cellDataOverview(this, time)

            bg = this.cellData(time).background;
            nucl = mean(this.cellData(time).nucLevel);
            area = mean(this.cellData(time).nucArea);

            disp(['number of cells: ' num2str(this.ncells(time))]);
            disp(['mean nuclear area: ' num2str(round(area)) ' pixels']);
            disp(['background: ' num2str(round(bg))]);
            disp(['mean nuclear level: ' num2str(round(nucl))]);
            if isfield(this.cellData(time),'cytLevel')
                cytl = nanmean(this.cellData(time).cytLevel);
                R = nanmean(this.cellData(time).NCratio);
                disp(['mean cytoplasmic level: ' num2str(round(cytl))]);
                disp(['N:C ' num2str(R,2)]);
                disp(['C:N ' num2str(1./R,2)]);
            end
        end

        % getter for dependent properties
        %---------------------------------

        function filename = get.filename(this)
            % filename without extension

            filename = sprintf(this.filenameFormat, this.ID-1);
        end

        function barefname = get.barefilename(this)
            % filename without extension

            barefname = this.filenameFormat(1:end-4);
        end

        % setters
        %---------------------------------

        function setID(this, ID)

            this.ID = ID;
        end
    end
end
