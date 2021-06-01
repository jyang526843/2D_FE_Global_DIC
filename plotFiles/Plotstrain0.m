function [x2,y2,disp_u,disp_v,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy, ...
    strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = Plotstrain0( ...
    U,F,Rad,x0,y0,sizeOfImg,DICpara)
%PLOTSTRAIN0: to compute and plot DIC solved strain fields
%   [x2,y2,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy, ...
%   strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = Plotstrain0(U,F,Rad,x0,y0,sizeOfImg,DICpara)
%
%   INPUT: U                DIC solved displacement vector
%          F                DIC solved deformation gradient tensor
%          x,y              x and y coordinates of each points on the image domain
%          SizeOfImg        Size of the DIC raw image
%          DICpara          DIC chosen parameters
%
%   OUTPUT: x2,y2                   x- and y-coordinates of points whose strain values are computed
%           disp_u,disp_v           Interpolated dispu and dispv at points {x2,y2}
%           dudx,dvdx,dudy,dvdy     E.g., dudx = d(disp_u)/dx at points {x2,y2}
%           strain_exx              strain xx-compoent
%           strain_exy              strain xy-compoent
%           strain_eyy              strain yy-compoent
%           strain_principal_max    max principal strain on the xy-plane
%           strain_principal_min    min principal strain on the xy-plane
%           strain_maxshear         max shear strain on the xy-plane
%           strain_vonMises         equivalent von Mises strain
%
%   Plots:       
%       1) strain sxx
%       2) strain sxy
%       3) strain syy
%       4) max principal strain on the xy-plane 
%       5) min principal strain on the xy-plane
%       6) max shear strain on the xy-plane
%       7) equivalent von Mises strain
%
%   TODO: users could change caxis range based on their own choices.
%
%
% Author: Jin Yang  (jyang526@wisc.edu)
% Last date modified: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%% Compute strain components

x = x0(1+Rad:end-Rad,1+Rad:end-Rad); 
y = y0(1+Rad:end-Rad,1+Rad:end-Rad);

M = size(x,1); N = size(x,2);
u_x = F(1:4:end); v_x = F(2:4:end);
u_y = F(3:4:end); v_y = F(4:4:end);
 
u_x = reshape(u_x,M,N); v_x = reshape(v_x,M,N);
u_y = reshape(u_y,M,N); v_y = reshape(v_y,M,N);
 
u = U(1:2:end); v = U(2:2:end);
u0 = reshape(u,M+2*Rad,N+2*Rad); v0 = reshape(v,M+2*Rad,N+2*Rad);
u = u0(1+Rad:end-Rad,1+Rad:end-Rad); v = v0(1+Rad:end-Rad,1+Rad:end-Rad);

% imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
if M < 9, x2 = x(:,1)'; else x2 = linspace(x(1,1),x(end,1),4*(length(x(:,1))-1)+1); end
if N < 9, y2 = y(1,:);  else y2 = linspace(y(1,1),y(1,end),4*(length(y(1,:))-1)+1); end
     

%% Compute displacement components to manipulate the reference image
disp_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);
disp_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);


%% Compute strain components
dudx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);
dvdx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_x,M*N,1),x2,y2);
dudy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_y,M*N,1),x2,y2);
dvdy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);

strain_exx = dudx;
strain_exy = 0.5*(dvdx + dudy);
strain_eyy = dvdy;

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);
         
% Please don't delete this line, to deal with the image and physical world coordinates       
[x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ====== 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,strain_exx,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Strain exy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,strain_exy,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,strain_eyy,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 4) Strain e_principal_max ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,strain_principal_max,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Principal strain $e_{\max}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 5) Strain e_principal_min ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2,strain_principal_min,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Principal strain $e_{\min}$','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 6) Strain e_max_shear ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2, strain_maxshear ,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Max shear strain','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 7) von Mises equivalent strain ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2,um2px*(sizeOfImg(2)+1)-y2, strain_vonMises ,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('von Mises equivalent strain','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x(M,N)<um2px*200, set(gca,'XTick',[]); end
if y(M,N)<um2px*200, set(gca,'YTick',[]); end
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
