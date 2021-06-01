% Compute image grayscale value difference (in pixels) between the 
% reference and deformed images
% 
% Author: Jin Yang, jyang526@wisc.edu; aldicdvc@gmail.com
% Date: 2020.11
%%

function PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized)

prompt = 'Do you want to compute image differences using current solved deformations? (0-yes; 1-no)';
ComputeImgDiffOrNot = input(prompt);

if ComputeImgDiffOrNot == 0
    
    x1=x0(1):1:x0(end); y1=y0(1):1:y0(end);
    %[y1Grid,x1Grid] = meshgrid(y1,x1);
    uFine = gridfit(x0,y0,u',x1,y1,'interp','bilinear'); uFine = uFine';
    vFine = gridfit(x0,y0,v',x1,y1,'interp','bilinear'); vFine = vFine';
    
    DispMask = zeros(length(x1),length(y1),2);
    DispMask(:,:,2) = uFine; DispMask(:,:,1) = vFine;
    f1 = fNormalized(x1(1):1:x1(end),y1(1):1:y1(end));
    % g1 = gNormalized(x1(1)-3:x1(end)+3,y1(1)-3:y1(end)+3);
    g3 = gNormalized(x1(1):1:x1(end),y1(1):1:y1(end));
    g2 = imwarp(g3,DispMask,'cubic');
     
    % g2 = ba_interp2(g1,vFine+y1Grid+3,uFine+x1Grid+3,'cubic');
    % pause;
    
    figure; surf(flip((f1-g2)',1),'edgecolor','none');view(2);colorbar;set(gca,'fontsize',18); % caxis([-1,1]); 
    colormap(gray); axis equal; axis tight; set(gcf,'color','w'); a = gca; a.TickLabelInterpreter = 'latex'; 
    title('Image grayscale difference before DIC','interpreter','latex');

    figure; surf(flip((f1-g3)',1),'edgecolor','none');view(2);colorbar;axis tight;set(gca,'fontsize',18); % caxis([-1,1]); 
    colormap(gray); axis equal; axis tight; set(gcf,'color','w'); a = gca; a.TickLabelInterpreter = 'latex'; 
    title('Image grayscale difference after DIC','interpreter','latex');

    sum(sum((f1-g2).^2)); 
    sum(sum((f1-g3).^2))
    
end

end