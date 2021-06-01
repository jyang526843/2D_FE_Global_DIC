function [stress_sxx, stress_sxy, stress_syy, stress_principal_max_xyplane, ...
    stress_principal_min_xyplane, stress_maxshear_xyplane, ...
    stress_maxshear_xyz3d, stress_vonMises] = Plotstress0(DICpara,ResultStrain,sizeOfImg)
%PLOTSTRESS0: to compute and plot DIC solved stress fields 
%   [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
%    stress_principal_min_xyplane, stress_maxshear_xyplane, ...
%    stress_maxshear_xyz3d, stress_vonMises]     = Plotstress0(DICpara,ResultStrain,sizeOfImg)
%
%   INPUT: DICpara          DIC para in the ALDIC code
%          ResultStrain     ALDIC computed strain field result
%          SizeOfImg        Size of the DIC raw image
%           
%   OUTPUT: stress_sxx      Cauchy stress xx-compoent
%           stress_sxy      Cauchy stress xy-compoent s_xy = s_yx
%           stress_syy      Cauchy stress yy-compoent
%           stress_principal_max_xyplane    max principal stress on the xy-plane
%           stress_principal_min_xyplane    min principal stress on the xy-plane
%           stress_maxshear_xyplane         max shear stress on the xy-plane
%           stress_maxshear_xyz3d           max shear stress on the xyz-three dimensional space
%           stress_vonMises                 von Mises stress
%
%   Plots:       
%       1) stress sxx
%       2) stress sxy
%       3) stress syy
%       4) max principal stress on the xy-plane 
%       5) min principal stress on the xy-plane
%       6) max shear stress on the xy-plane
%       7) max shear stress on the xyz-three dimensional space
%       8) equivalent von Mises stress
%
%   TODO: users could change caxis range based on their own choices.
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2020.12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)
   

%% Load computed strain fields
x2 = ResultStrain.strainxCoord; 
y2 = ResultStrain.strainyCoord;
[M,N] = size(x2);
dudx = ResultStrain.dudx; dvdx = ResultStrain.dvdx; 
dudy = ResultStrain.dudy; dvdy = ResultStrain.dvdy; 
strain_exx = dudx; 
strain_exy = 0.5*(dvdx + dudy); 
strain_eyy = dvdy;


%% Load displacement components to deform the reference image
disp_u = ResultStrain.dispu; 
disp_v = ResultStrain.dispv;


%% Compute stress components

% ------ For each material model ------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DICpara.MaterialModel == 1 
    
    % Linear elasticity -- Plane stress
    % stress_xyz = [sxx  sxy  0;      strain_xyz = [exx  exy  0;
    %               sxy  syy  0;                    exy  eyy  0;
    %               0    0    0];                   0    0    ezz];
    E = DICpara.MaterialModelPara.YoungsModulus;
    nu = DICpara.MaterialModelPara.PoissonsRatio;
    stress_sxx = E/(1-nu^2)*(strain_exx + nu*strain_eyy);
    stress_syy = E/(1-nu^2)*(nu*strain_exx + strain_eyy);
    stress_sxy = E/(1+nu)*strain_exy;
     
    % Principal stress
    stress_maxshear_xyplane = sqrt((0.5*(stress_sxx-stress_syy)).^2 + stress_sxy.^2);
    stress_principal_max_xyplane = 0.5*(stress_sxx+stress_syy) + stress_maxshear_xyplane;
    stress_principal_min_xyplane = 0.5*(stress_sxx+stress_syy) - stress_maxshear_xyplane;
    stress_maxshear_xyz3d = max(stress_maxshear_xyplane, 0.5*abs(stress_principal_max_xyplane));
    
    % von Mises stress
    stress_vonMises = sqrt(0.5*( (stress_principal_max_xyplane-stress_principal_min_xyplane).^2 + ...
        (stress_principal_max_xyplane).^2 + (stress_principal_min_xyplane).^2 ));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif DICpara.MaterialModel == 2 
    
    % Linear elasticity -- Plane strain
    % stress_xyz = [sxx  sxy  0;      strain_xyz = [exx  exy  0;
    %               sxy  syy  0;                    exy  eyy  0;
    %               0    0    szz];                 0    0    0]; 
    E = DICpara.MaterialModelPara.YoungsModulus;
    nu = DICpara.MaterialModelPara.PoissonsRatio;
    stress_sxx = E*(1-nu)/(1+nu)/(1-2*nu)*(strain_exx + nu/(1-nu)*strain_eyy);
    stress_syy = E*(1-nu)/(1+nu)/(1-2*nu)*(strain_eyy + nu/(1-nu)*strain_exx);
    stress_sxy = E/(1+nu)*strain_exy;
    stress_szz = nu*(stress_sxx + stress_syy);

    % Principal stress
    stress_maxshear_xyplane = sqrt((0.5*(stress_sxx-stress_syy)).^2 + stress_sxy.^2);
    stress_principal_max_xyplane = 0.5*(stress_sxx+stress_syy) + stress_maxshear_xyplane;
    stress_principal_min_xyplane = 0.5*(stress_sxx+stress_syy) - stress_maxshear_xyplane;
    stress_maxshear_xyz3d = reshape(  max( [ stress_maxshear_xyplane(:), 0.5*abs(stress_principal_max_xyplane(:)-stress_szz(:)), ...
                                0.5*abs(stress_principal_min_xyplane(:)-stress_szz(:)) ], [], 2 ),  size(stress_maxshear_xyplane) ) ;
                            
    % von Mises stress
    stress_vonMises = sqrt(0.5*( (stress_principal_max_xyplane-stress_principal_min_xyplane).^2 + ...
        (stress_principal_max_xyplane-stress_szz).^2 + (stress_principal_min_xyplane-stress_szz).^2 ));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif DICpara.MaterialModel == 3 % Neo-Hookean or other material models
    
    disp('User needs to modify the code by yourself.'); pause; 
    % TODO: Please modify by yourself
    
end

 
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Stress sxx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2, um2px*(sizeOfImg(2)+1)-y2,stress_sxx,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Stress $s_{xx}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Strain sxy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_sxy,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Stress $s_{xy}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain syy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_syy,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Stress $s_{yy}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 4) Strain stress_principal_max_xyplane ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_principal_max_xyplane,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('$xy$-plane principal stress $s_{\max}$','FontWeight','Normal','Interpreter','latex'); 
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 5) Strain stress_principal_min_xyplane ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_principal_min_xyplane,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('$xy$-plane principal stress $s_{\min}$','FontWeight','Normal','Interpreter','latex'); 
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 6) Strain stress_maxshear_xyplane ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_maxshear_xyplane,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('$xy$-plane max shear stress','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 7) Strain stress_maxshear_xyz3d ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_maxshear_xyz3d,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('$xyz$-3D max shear stress','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 8) von Mises stress ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,stress_vonMises,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('von Mises equivalent stress','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
 
 