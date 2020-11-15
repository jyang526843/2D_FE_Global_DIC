% ==============================================
% function ReadImage
% ----------------------------------------------
%   This function is to load DIC images
%   Images can be loaded by selecting images folder, inputing image name
%   prefix, or directly selecting images manually.
%
% Author: Jin Yang.
% Last updated 02/2020.
% ==============================================
function [file_name,Img,DICpara] = ReadImage(varargin)

fprintf('Choose method to load images:  \n')
fprintf('     0: Select images folder;  \n')
fprintf('     1: Use prefix of image names;  \n')
fprintf('     2: Manually select images.  \n')
prompt = 'Input here: ';
LoadImgMethod = input(prompt);

switch LoadImgMethod 
    case 0
        % ==============================================
        imgfoldername = uigetdir(pwd,'Select images folder');
        addpath([imgfoldername,'\']);
        img1 = dir(fullfile(imgfoldername,'*.jpg'));
        img2 = dir(fullfile(imgfoldername,'*.jpeg'));
        img3 = dir(fullfile(imgfoldername,'*.tif'));
        img4 = dir(fullfile(imgfoldername,'*.tiff'));
        img5 = dir(fullfile(imgfoldername,'*.bmp'));
        img6 = dir(fullfile(imgfoldername,'*.png'));
        img7 = dir(fullfile(imgfoldername,'*.jp2'));
        file_name = [img1;img2;img3;img4;img5;img6;img7];
        file_name = struct2cell(file_name);
    case 1
        % ==============================================
        fprintf('What is prefix of DIC images? E.g. img_0*.tif.   \n')
        prompt = 'Input here: ';
        file_name = input(prompt,'s');
        [~,imgname,imgext] = fileparts(file_name);
        file_name = dir([imgname,imgext]);
        file_name = struct2cell(file_name);
    otherwise
        % ==============================================
        disp('--- Please load first image ---')
        file_name{1,1} = uigetfile('*.tif','Select reference Image (Deformed)');
        disp('--- Please load next image ---')
        file_name{1,2} = uigetfile('*.tif','Select deformed Image (Reference)');
        prompt = 'Do you want to load more deformed images? (0-Yes; 1-No)';
        DoYouWantToLoadMoreImages = input(prompt); imageNo = 2;
        while ( DoYouWantToLoadMoreImages == 0 )   
            imageNo = imageNo + 1;
            file_name{1,imageNo} = uigetfile('*.tif','Select Deformed Image');
            prompt = 'Do you want to load more deformed images? (0-Yes; 1-No)';
            DoYouWantToLoadMoreImages = input(prompt);
        end
end

% ==============================================
% The following codes only consider two images comparasion
numImages = size(file_name,2);
for i = 1:numImages
    Img{i} = imread(file_name{1,i});
    % Change color RGB images to grayscale images
    [~, ~, numberOfColorChannels] = size(Img{i});
    if (numberOfColorChannels==3)
        Img{i} = rgb2gray(Img{i});
    end
    Img{i} = double(Img{i})';
end


%f = imread(file_name{1});
%g = imread(file_name{2});

% % Decide rgb images or grayscale images
%if (numberOfColorChannels==3)
%    f = rgb2gray(f);
%    g = rgb2gray(g);
%end

% % Field of the undeformed and deformed greyscale;
%f = double(f'); g = double(g');

% images are:
% --------------------
% --  This is image --
% |                  |
% y                  |
% |                  |
% |  --> x direction |
% |                  |
% --------------------
% after transforming,  MatLab matrix direction:
% --  This is matrix in Matlab --
% |                             |
% x                             |
% |                             |
% |  --> y direction            |
% |                             |
% --------------------------------
 
% ==============================================
% Decide DIC subset parameters
% Choose ZOI
fprintf('\n');
disp('--- Define ROI corner points at the top-left and the bottom-right ---')
imshow( (imread(file_name{1}))); 
title('Click top-left and the bottom-right corner points','fontweight','normal','fontsize',16);

gridx = zeros(1,2); gridy = zeros(1,2);
[gridx(1), gridy(1)] = ginput(1);
fprintf('Coordinates of top-left corner point are (%4.3f,%4.3f)\n',gridx(1), gridy(1))

[gridx(2), gridy(2)] = ginput(1);
fprintf('Coordinates of bottom-right corner point are (%4.3f,%4.3f)\n',gridx(2), gridy(2))

gridxy.gridx = round(gridx); gridxy.gridy = round(gridy);

% Choose subset size
fprintf('\n');
fprintf('--- What is the finite element size? --- \n');
prompt = 'Input here: ';
winstepsize = input(prompt);

winsize = winstepsize + 10;
 
% ==============================================
% Deal with image sequence
NewFFTSearch = 1; % By default initialize parameters
if numImages > 2
    
    % ==============================================
    % DIC initial guess 
    NewFFTSearch = funParaInput('NewFFTSearch'); % Use last frame as init guess or not
    
    % ==============================================
    % Decide DIC as accumulative or incremental mode?
    fprintf('--- Choose accumulative or incremental mode ---  \n')
    fprintf('     0: Accumulative(By default);  \n')
    fprintf('     1: Incremental;  \n')
    prompt = 'Input here: ';
    DICIncOrNot = input(prompt);

    try
        switch DICIncOrNot
            case 0
                ImgSeqIncUnit = numImages+1;
                ImgSeqIncROIUpdateOrNot = 1;
            case 1
                fprintf('Incremental mode: How many frames to update reference image once? \n');
                prompt = 'Input here: ';
                ImgSeqIncUnit = input(prompt);
                fprintf('Update ROI at the same time of updating reference image? \n');
                fprintf('    0: Do not update ROI; \n'); 
                fprintf('    1: Manually(Recommended); \n'); 
                fprintf('    2: Automatically; \n'); 
                prompt = 'Input here: ';
                ImgSeqIncROIUpdateOrNot = input(prompt);
            otherwise
                ImgSeqIncUnit = numImages+1;
                ImgSeqIncROIUpdateOrNot = 1;
        end
         
    catch
        ImgSeqIncUnit = numImages+1; 
        ImgSeqIncROIUpdateOrNot = 1;
    end
    
    
    
% ================================    
else % Only two frames
    
    ImgSeqIncUnit = numImages+1; 
    ImgSeqIncROIUpdateOrNot = 1;
     
end

DICpara.winsize = winsize;
DICpara.winstepsize = winstepsize;
DICpara.gridxyROIRange = gridxy;
DICpara.LoadImgMethod = LoadImgMethod;
DICpara.ImgSeqIncUnit = ImgSeqIncUnit;
DICpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot;
DICpara.ImgSize = size(Img{1});

end