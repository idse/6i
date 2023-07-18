function lims = seglim(varargin)

if nargin == 3
    im = varargin{1};
    seg = varargin{2};
    n = varargin{3};
elseif nargin == 2
    im = varargin{1};
    seg = varargin{2};
    n = 2;
elseif nargin == 1
    im = varargin{1};
    n = 2;
%     counts = imhist(im(im>0));
%     T = otsuthresh(counts);
%     seg = imbinarize(im,T);
    idxs = find(im==0); nzidxs = find(im>0); test = im;
    test(idxs) = test(nzidxs(1:length(idxs)));
    seg = imbinarize(im,adaptthresh(test));
    seg = imopen(imclose(seg,strel('disk',1)),strel('disk',3));
end

im = im2double(im);
fg = im(seg);
bg = im(~seg & im>0);

fgmean = mean(fg); fgstd = std(fg);
bgmean = mean(bg); bgstd = std(bg);

lims = [max(0,bgmean-n*bgstd), min(fgmean+n*fgstd,1)];

end