function img = loadTiffStack(fname)

t = imfinfo(fname);
nz = length(t);
m = t(1).Height; n = t(1).Width;

if t(1).BitDepth == 16
    imclass = 'uint16';
elseif t(1).BitDepth == 8
    imclass = 'uint8';
end

img = zeros(m,n,nz,imclass);
for zi = 1:nz
    img(:,:,zi) = imread(fname,zi);
end

end