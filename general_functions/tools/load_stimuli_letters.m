function So = load_stimuli_letters(experiment)
%image size corresponds to the size of the image in the implant coordinate
%frame

letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
add_directories;
Stimuli=  [Stimuli_folder, '/letters'];

stim_size = floor(0.9*experiment.display_size(1));
Si = load_Snellen(experiment.display_width, experiment.display_height, [Stimuli_folder, '/letters'], stim_size , letters);
So = NaN(experiment.ny*experiment.nx,numel(letters));
N = numel(letters);
for i = 1:N
    S = imresize(squeeze(Si(i,:,:)), [experiment.ny, experiment.nx], 'method', 'bilinear');
%     S = load_Snellen(experiment.display_width, experiment.display_height, Stimuli, stim_size, letters(i));

    So(:,i)=S(:);
    
    
    
end
% files = dir(stimuli_folder);
% S=[];
% k=0;
% n= length(files);
%
% im_ny = image_size(1);
% im_nx = image_size(2);
%
% padsize = floor([ny - im_ny,nx - im_nx]/2);
%
% for j = 1:n
%     if ~ files(j).isdir
%         if strcmp(files(j).name(end-2:end), 'mat')
%             k=k+1;
%             image = load([stimuli_folder,'/',files(j).name]);
%             image= struct2cell(image);
%             image = image{1};
%             extracted_image= double(imresize(image, [im_ny, im_nx],  'method', 'bilinear'));
%             s = padarray(extracted_image, padsize, 0, 'both');
%             s = imresize(s, [ny, nx], 'method','bilinear');
%             s= double(reshape(s, [],1));
%             S=[S, s];
%         end
%     end
% end
%

