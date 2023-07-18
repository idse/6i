function [shift, scale] = zAlignImageStacks(im1,im2)

if ~isa(im1,'double')
    im1 = double(im1);
end
if ~isa(im2,'double')
    im2 = double(im2);
end

nz1 = size(im1,3); nz2 = size(im2,3);
Rs = NaN(nz1,nz2);

tic
for z1 = 1:nz1
    fprintf('.')
    for z2 = 1:nz2
        test1 = im1(:,:,z1); test2 = im2(:,:,z2); mask = test1 > 0 & test2 > 0;
        R = corrcoef(test1(mask),test2(mask));
        Rs(z1,z2) = R(2);
    end
end
fprintf('\n')
toc

z1idxs = NaN(nz1,1);
corrs1 = NaN(nz1,1);
for zi = 1:nz1
    [mm,I] = max(Rs(zi,:));
    z1idxs(zi,1) = I;
    corrs1(zi,1) = mm;
end

zidxs = NaN(nz2,1);
zcorr = NaN(nz2,1);
corrs = NaN(nz2,1);
for zi = 1:nz2
    [mm,I] = max(Rs(:,zi));
    zidxs(zi,1) = I;
    corrs(zi,1) = mm;
    
    if z1idxs(I,1) == zi && abs(I - zi) < 25
        zcorr(zi,1) = I;
    end
end

% find a fit line
z2s = find(~isnan(zcorr(:,1)));
z1s = zcorr(z2s);
b = [ones(size(z2s)),z2s]\z1s;
shift = b(1); scale = b(2);
fprintf('shift = %g, scale = %g\n',shift,scale)



end