function clabel = renameDuplicateChannels(channelLabels)

ulabel = unique(channelLabels);
clabel = channelLabels;

for ii = 1:length(ulabel)
    mask = strcmp(channelLabels,ulabel{ii});
    ids = find(mask);
    if sum(mask) > 1
        for jj = 1:length(ids)
            clabel{ids(jj)} = [ulabel{ii},'_',num2str(jj)];
        end
    end
end


end
