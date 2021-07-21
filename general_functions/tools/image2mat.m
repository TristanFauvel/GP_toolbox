%% Script to convert images to mat matrices in a fodler
folder = '/home/tristan/Desktop/optim_retina/Stimuli/landoltC';
files = dir(folder);
n= length(files);
k=0;
for j = 1:n
    if ~ files(j).isdir
        if ~ strcmp(files(j).name(end-2:end),'mat')
            k=k+1;
            image = imread([folder,'/',files(j).name]);
            image =rgb2gray(image); % Convert to grey
%             image =imresize(image, [ny, nx]);
            save([folder,'/',files(j).name(1:end-3),'mat'], 'image')
%             figure()
%             imagesc(image)
        end
    end
end
