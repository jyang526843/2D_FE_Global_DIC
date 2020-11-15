%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot DIC solved strain components on the original DIC images
%   1) strain exx
%   2) strain exy
%   3) strain eyy
%   4) principal strain max
%   5) principal strain min
%   6) max shear
%   7) equivalent von Mises strain
%
% Author: Jin Yang  
% Last date modified: 2020.11.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function Plotstrain(U,Rad,FSubpb3,x,y,sizeOfImg,CurrentImg,OrigDICImgTransparency)

warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

M = size(x,1); N = size(x,2);
u_x = FSubpb3(1:4:end); v_x = FSubpb3(2:4:end);
u_y = FSubpb3(3:4:end); v_y = FSubpb3(4:4:end);
 
u_x = reshape(u_x,M,N); v_x = reshape(v_x,M,N);
u_y = reshape(u_y,M,N); v_y = reshape(v_y,M,N);

u = U(1:2:end); v = U(2:2:end);
u0 = reshape(u,M+2*Rad,N+2*Rad); v0 = reshape(v,M+2*Rad,N+2*Rad);
u = u0(1+Rad:end-Rad,1+Rad:end-Rad); v = v0(1+Rad:end-Rad,1+Rad:end-Rad);

% imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
if M < 9,x2 = x(:,1)'; else x2 = interp(x(:,1)',4); end
if N < 9,y2 = y(1,:);  else y2 = interp(y(1,:),4);  end


disp_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);
disp_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);


strain_exx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);
strain_exy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(0.5*(v_x+u_y),M*N,1),x2,y2);
strain_eyy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);

[x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';


%% ====== 1) Strain exx ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_exx,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% ====== 2) Strain exy ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_exy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Strain $e_{xy}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% ====== 3) Strain eyy ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_eyy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Strain $e_{yy}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');

 

%% ====== 4) Strain e_principal_max ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_principal_max,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$Principal strain e_{\max}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% ====== 5) Strain e_principal_min ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_principal_min,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$Principal strain e_{\min}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% ====== 6) Strain e_max_shear ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_maxshear,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Max shear strain','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');

 
%% ====== 7) von Mises equivalent strain ======
fig1=figure; ax1=axes; 
try
    h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch
    h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; h2=surf(x2+disp_u,sizeOfImg(2)+1-(y2-disp_v),strain_vonMises,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;
colormap(coolwarm(32));  colormap jet; colormap(cMap);
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray'); % colormap(ax2,coolwarmmap);  colormap( coolwarmmap);

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('von Mises equivalent strain','FontWeight','Normal','Interpreter','latex');
 




end
 
 