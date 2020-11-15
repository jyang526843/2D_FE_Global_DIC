% ==============================================
% function funImgGradient
% ==============================================

function [Df] = funImgGradient(fNormalized,gNormalized,varargin)
%[imgfNormalizedbc,imggNormalizedbc,imgSize,DfAxis] = funImgGradient(fNormalized,gNormalized,varargin)

imgSize = size(gNormalized);
disp('--- Start to compute image gradients ---');

if nargin > 2
    method = varargin{2};
else
    method = [];
end

if strcmp(method,'Splines_interp')
    
    imgfNormalizedbc = Spline2D('bicubic',[1:1:size(fNormalized,1)],[1:1:size(fNormalized,2)],fNormalized);
    imggNormalizedbc = Spline2D('bicubic',[1:1:size(gNormalized,1)],[1:1:size(gNormalized,2)],gNormalized);
    % [XX,YY] = ndgrid([1:1:size(fNormalized,1)],[1:1:size(fNormalized,2)]);
    % DfDxNormalizedbc = imgfNormalizedbc.eval_Dx(XX,YY);
    % DfDyNormalizedbc = imgfNormalizedbc.eval_Dy(XX,YY);
    DfAxis = [0,size(fNormalized,1)-1,0,size(fNormalized,2)-1];
    
else
    
    imgradientMatrix = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]';
    DfDxStartx = 4; % max([4,coordinatesFEM(1,1)-0.5*winsize-1]);
    DfDxStarty = 4; % max([4,coordinatesFEM(1,2)-0.5*winsize-1]);
    DfDxEndx = size(fNormalized,1)-3; % min([coordinatesFEM(M*N,1)+0.5*winsize,size(fNormalized,1)-3]);
    DfDxEndy = size(fNormalized,2)-3; % min([coordinatesFEM(M*N,2)+0.5*winsize,size(fNormalized,2)-3]);
    % DfDxStartx = coordinates(1,1)-1; DfDxStarty = coordinates(1,2)-1;
    % I = fNormalized( coordinates(1,1)-3:coordinates((M+1)*(N+1),1)+3, coordinates(1,2)-3:coordinates((M+1)*(N+1),2)+3);
    I = fNormalized( DfDxStartx-3:DfDxEndx+3, DfDxStarty-3:DfDxEndy+3 );
    
    DfDxNormalizedtemp = imfilter(I, imgradientMatrix);
    DfDyNormalizedtemp = imfilter(I, imgradientMatrix');
    DfDxNormalized = DfDxNormalizedtemp(4:end-3, 4:end-3);
    DfDyNormalized = DfDyNormalizedtemp(4:end-3, 4:end-3);
    
    DfAxis = [DfDxStartx,DfDxEndx,DfDxStarty,DfDxEndy]-1;
    Df.DfAxis = DfAxis; Df.imgSize = imgSize;
    Df.DfDx = DfDxNormalized;
    Df.DfDy = DfDyNormalized;
    
end

disp('--- Computing image gradients done ---');