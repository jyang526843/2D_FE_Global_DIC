%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot DIC solved displacement components on the original DIC images
%   1) dispx 
%   2) dispy
%
% Author: Jin Yang  
% Last date modified: 2020.11.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function Plotdisp(U,x,y,sizeOfImg,CurrentImg,OrigDICImgTransparency)
  
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
  
M = size(x,1); N = size(x,2);
u = U(1:2:end); v = U(2:2:end);
u = reshape(u,M,N); v = reshape(v,M,N);
 
if M < 9, x2 = x(:,1)'; else x2 = linspace(x(1,1),x(end,1),4*(length(x(:,1))-1)+1); x2=x2(:)'; end  
if N < 9, y2 = y(1,:); else y2 = linspace(y(1,1),y(1,end),4*(length(y(1,:))-1)+1); y2=y2(:)'; end  


disp_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);
disp_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);


[x2,y2]=ndgrid(x2,y2);x2=x2'; y2=y2';


%% ====== 1) dispx u ======
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u ,sizeOfImg(2)+1-(y2-disp_v),disp_u,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(coolwarm(32)); colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  % Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');


 
%% ====== 2) dispy v ======
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+disp_u ,sizeOfImg(2)+1-(y2-disp_v),disp_v,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(coolwarm(32)); colormap jet; colormap(cMap);

linkaxes([ax1,ax2]);  % Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



end
