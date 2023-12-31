function [im, ncomps] = visualize_nuclei_v2(seg, img, s)

if ~exist('img','var')
    img = zeros(size(seg));
elseif ~isa(img,'double')
    img = im2double(img);
end

if ~exist('s','var')
    s = 0.5;
end

if size(seg,3) > 1
    seg = seg(:,:,1);
end

stats = regionprops(seg,'PixelIdxList');
idxlists = {stats.PixelIdxList};

% idxlists = {};
% vals = unique(seg);
% vals = vals(vals > 0);
% for vi = 1:length(vals)
%     CC = bwconncomp(seg == vals(vi),4);
%     idxlists = [idxlists, CC.PixelIdxList];
% end

ncomps = length(idxlists);
colors = zeros(ncomps,3);
for ci = 1:9000:ncomps
    idxs = ci:min(ncomps,ci+8999);
    colors(idxs,:) = distinguishable_colors(length(idxs),'k');
end
% colors = distinguishable_colors(ncomps,'k');
im1 = zeros(size(seg)); im2 = im1; im3 = im1;
for li = 1:ncomps
    im1(idxlists{li}) = colors(li,1);
    im2(idxlists{li}) = colors(li,2);
    im3(idxlists{li}) = colors(li,3);
end

im = s*cat(3,im1,im2,im3) + (1 - s)*repmat(img,1,1,3);

end

