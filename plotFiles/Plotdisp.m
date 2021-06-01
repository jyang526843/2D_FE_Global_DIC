function Plotdisp(U,x,y,sizeOfImg,CurrentImg,varargin)
%FUNCTION Plotdisp(U,x,y,sizeOfImg,CurrentImg,DICpara)
% To plot DIC solved displacement components
% ----------------------------------------------
%
%   INPUT: U                 Displacement vector: U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          x,y               FE mesh x- and y-coordinates
%          sizeOfImg         DIC deformed image size
%          CurrentImg        Current deformed image
%          DICpara           varargin{1}: DIC paramters
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
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
% Last time updated: 11/2020.
% ==============================================


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DICpara] = parseargs(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

M = size(x,1); N = size(x,2);
u = U(1:2:end); v = U(2:2:end);
u = reshape(u,M,N); v = reshape(v,M,N);
 
if M < 9, x2 = x(:,1)'; else x2 = linspace(x(1,1),x(end,1),4*(length(x(:,1))-1)+1); end
if N < 9, y2 = y(1,:);  else y2 = linspace(y(1,1),y(1,end),4*(length(y(1,:))-1)+1); end


%% Compute displacement components to manipulate the reference image
disp_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);
disp_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);
 
% Please don't delete this line, to deal with the image and physical world coordinates 
[x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) dispx u ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1 = imshow(  flipud( imread(CurrentImg) ),'InitialMagnification','fit');
catch h1 = surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end
 
axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2 = surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    (sizeOfImg(2)+1)-((y2-Image2PlotResults*disp_v))/um2px, disp_u, 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-1.7,.2]); % caxis([-0.025,0.025]); 
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

% xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) dispx v ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2 = surf( (x2+Image2PlotResults*disp_u)/um2px, (sizeOfImg(2)+1)-((y2-Image2PlotResults*disp_v))/um2px, disp_v, 'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([5,12]); % caxis([-0.025,0.025]); 
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

% xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');


end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [DICpara] = parseargs(vargin)
DICpara=[];  
 
try 
    DICpara=vargin{1};
     
catch
    
end

end

