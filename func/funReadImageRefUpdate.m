% ==============================================
% function ReadImage
% ==============================================
function [gridxy] = funReadImageRefUpdate(filename)

% fprintf('Choose method to load images:  \n')
% fprintf('     0: Select images folder;  \n')
% fprintf('     1: Use prefix of image names;  \n')
% fprintf('     2: Manually select images.  \n')
% prompt = 'Input here: ';
% LoadImgMethod = input(prompt);
 
% ==============================================
% Choose ZOI
fprintf('***** Update incremental reference image ***** \n');
disp('--- Define ROI corner points at the top-left and the bottom-right ---')
imshow( (imread(filename)));
title('Click top-left and the bottom-right corner points','fontweight','normal','fontsize',16);

gridx = zeros(1,2); gridy = zeros(1,2);
[gridx(1), gridy(1)] = ginput(1);
fprintf('Coordinates of top-left corner point are (%4.3f,%4.3f)\n',gridx(1), gridy(1))

[gridx(2), gridy(2)] = ginput(1);
fprintf('Coordinates of bottom-right corner point are (%4.3f,%4.3f)\n',gridx(2), gridy(2))

gridxy.gridx = round(gridx); gridxy.gridy = round(gridy);

% % Choose subset size
% prompt = '--- What is the subset size? --- Input here: ';
% winsize = input(prompt);
% 
% % Choose subset size
% prompt = '--- What is the subset step? --- Input here: ';
% winstepsize = input(prompt);
 


end