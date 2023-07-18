clear; close all; clc

%%
%periodically check for new image files and transfer the files as they are 
%written so processing can begin before imaging is completed

%specify the source directory
dir1 = '';
%specify the target directory
dir2 = '';
%create the target directory if it does not exist
if ~exist(dir2,'dir'), mkdir(dir2); end

%pause between finding and transferring files long enough to allow any
%position currently being imaged to complete imaging
nseconds = 120;
%total number of files we expect to see
totalfiles = (6*24*4);

total_transferred = 0; %running total of files transferred so far
filenames = cell(1,totalfiles); %list of all filenames, updated as we find them
transferred = false(1,totalfiles); %boolean to check whether each file has been transferred yet
while total_transferred < totalfiles
    t=datetime('now');
    disp("Starting next loop, start time:")
    disp(t)
    %find nd2 files in the source directory
    list1 = dir(fullfile(dir1,'*.nd2'));
    disp('paused')
    pause(nseconds)
    %iterate over source files to determine whether to transfer each one
    disp("Copying files")
    for fi = 1:length(list1)
        name = list1(fi).name;
        idx = find(strcmp(name,filenames)); %check previous filenames for matches with this one
	    %if this is a new file or has not yet been transferred
        if isempty(idx) || ~transferred(idx)
            if isempty(idx)
                %no matches; this is a new file
                idx = find(cellfun(@isempty,filenames),1,'first');
                filenames{idx} = name;
            end
            
            fprintf('on file %d\n',fi)
            copyfile(fullfile(dir1,name),fullfile(dir2,name))%,'f')
            transferred(idx) = true;
            total_transferred = sum(transferred);
                
        end
    end
    fprintf('%d files found, %d files transferred\n',sum(~cellfun(@isempty,filenames)),total_transferred)
end

fprintf('completed; %d files found, %d files transferred\n',sum(~cellfun(@isempty,filenames)),total_transferred)

