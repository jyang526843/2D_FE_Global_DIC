function Plotuv(USubpb2,x,y )
 
M = size(x,1); N = size(x,2);
u = USubpb2(1:2:end); v = USubpb2(2:2:end);
u = reshape(u,M,N); v = reshape(v,M,N);
figure; surf(x,y, u);  set(gca,'fontSize',18); 
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
axis tight; % set(gca,'XTick',[] );
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
box on; colormap jet;
 

figure; surf(x,y,  v);  set(gca,'fontSize',18);  
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
axis tight; % set(gca,'XTick',[] );
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
box on; colormap jet;

% print -painters -dpng -r600  fig_Sample14_Local_u.png


% figure; mesh(x,y,reshape( Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact x displacement');
% axis(l);
% % print -painters -dpng -r600  fig_Sample14_Local_uExact.png
% 
% 
% figure; mesh(x,y, v);  set(gca,'fontSize',20); view(-5, 30);title('Computed y displacement'); axis tight
% l = axis;
% % print -painters -dpng -r600  fig_Sample14_Local_v.png
% 
% 
% figure; mesh(x,y,reshape(0* Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact y displacement');
% axis(l);
% % print -painters -dpng -r600  fig_Sample14_Local_vExact.png