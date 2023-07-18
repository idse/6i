function writeTiffZStack(img, fname)

nz = size(img,3);
for zi = 1:nz
    if zi == 1
        mode = 'overwrite';
    else
        mode = 'append';
    end
    imwrite(img(:,:,zi),fname,'WriteMode',mode);
end

end