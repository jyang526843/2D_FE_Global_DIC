%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot DIC solved displacement components
%   1) dispx 
%   2) dispy
%
% Author: Jin Yang  
% Last date modified: 2020.11.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plotdisp_show(U0,coordinatesFEM,elementsFEM)

warning off; 

%% ====== 1) dispx u ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,U0(1:2:end)); 
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar; % view([90 -90])

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; box on;


%% ====== 2) dispy v ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,U0(2:2:end)); 
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; box on;

end

