function [x2,y2,disp_u,disp_v,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy, ...
    strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = Plotstrain( ...
    U,F,Rad,x0,y0,sizeOfImg,CurrentImg,DICpara)
%PLOTSTRAIN: to compute and plot DIC solved strain fields on the original DIC images
%   [x2,y2,disp_u,disp_v,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy, ...
%   strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = Plotstrain( ...
%   U,Rad,F,x0,y0,sizeOfImg,CurrentImg,DICpara)
%
%   INPUT: U                DIC solved displacement fields
%          F                DIC solved deformation gradient tensor
%          Rad              Parameter used to compute the final F. If we use the direct
%                           output from the ALDIC code, Rad=0. If we apply a finite  
%                           difference or plane fitting of the U to compute the F, Rad>0.
%          x0,y0            x and y coordinates of each points on the image domain
%          SizeOfImg        Size of the DIC raw image
%          CurrentImg       File name of current deformed image
%          DICpara          DIC para in the ALDIC code
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
% Author: Jin Yang  (jyang526@wisc.edu)
% Last date modified: 2020.12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

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
% 2D equivalent von Mises strain
strain_vonMises = 2/3 * sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min );
% There was a mistake before in computing "strain_vonMises":   previous "+ 3*strain_maxshear.^2" term should not be included.
% Please be careful on computing von Mises strain, there can be a different formula to use 
% Some people define it as: " sqrt(strain_principal_max.^2 + strain_principal_min.^2 ...
%                                   - strain_principal_max.*strain_principal_min );"
% Link: https://help.csiamerica.com/help/sap2000/26/26.0.0/SAP2000/WebHelp/Menus/Display/Show_Forces_Stresses/Von_Mises_Stress.htm

% Please don't delete this line, to deal with the image and physical world coordinates       
[x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_exx, 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); %caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Strain exy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_exy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.008,0.008]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_eyy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.002,0.014]);
% caxis([-0.55,0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 4) Strain e_principal_max ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_principal_max,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.002,0.014]) % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 5) Strain e_principal_min ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_principal_min,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 6) Strain e_max_shear ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_maxshear,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 7) von Mises equivalent strain ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, strain_vonMises,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto;  
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;  colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.022]) % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 

end
 
 
