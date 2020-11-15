% ==============================================
% function funNormalizeImg
% ----------------------------------------------
%   This function is to normalize images with the formula: 
%       Normalized(Img) = (Img-Mean(Img))/sqrt(Std(Img))
%   
%   Normalized image grayscale value is between [-1,1]
%
% Author: Jin Yang.
% Last updated 02/2020.
% ==============================================
function [ImgNormalized,gridxyROIRange] = funNormalizeImg(Img,gridxyROIRange)

gridxROIRange = gridxyROIRange.gridx; gridyROIRange = gridxyROIRange.gridy; 

% ======= To garantuee ROI value is within image boundary =======
if gridxROIRange(1)< 1; gridxROIRange(1) = 1; end
if gridxROIRange(2)>size(Img{1},1); gridxROIRange(2) = size(Img{1},1);end
if gridyROIRange(1)< 1; gridyROIRange(1) = 1; end
if gridyROIRange(2)>size(Img{1},2); gridyROIRange(2) = size(Img{1},2);end 

gridxyROIRange.gridx = gridxROIRange; gridxyROIRange.gridy = gridyROIRange;
 
% ============== Normalize or compress images ==============
% ============== Normalize f and g ==============
ImgNormalized = cell(size(Img));

for i = 1:length(Img)
   tempf = Img{i}(gridxROIRange(1):gridxROIRange(2), gridyROIRange(1):gridyROIRange(2));
   fAvg = mean(tempf(:)); fstd = std(tempf(:));
   ImgNormalized{i} = (Img{i}-fAvg)/fstd;
end


% % Normalize f and g
% tempf = f(gridxROIRange(1):gridxROIRange(2), gridyROIRange(1):gridyROIRange(2));
% fAvg = mean(tempf(:)); fstd = std(tempf(:));
% 
% tempg = g(gridxROIRange(1):gridxROIRange(2), gridyROIRange(1):gridyROIRange(2));
% gAvg = mean(tempg(:)); gstd = std(tempg(:));
% 
% fNormalized = (f-fAvg)/fstd;
% gNormalized = (g-gAvg)/gstd;

% ========= Don't Normalized images and use original images =========
% fNormalized = f; gNormalized = g;

% ========= Do wavelets compression =========
% [fNormalized] = func_compress_dw2d(fNormalized0,'sym4',0.1);
% [gNormalized] = func_compress_dw2d(gNormalized0,'sym4',0.1);