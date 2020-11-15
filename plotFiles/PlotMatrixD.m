%% Plot D matrix
close all; 
tempD = DerivativeOp(10,10,1);
[tempy,tempx] = meshgrid(1:1:size(tempD,1),1:1:size(tempD,2));
figure; imshow(int8(full(tempD+1)*0.5*128)); axis tight;
b=colorbar; caxis([0,128])
b.Ticks = [0,64,128]
b.TickLabels = {'-1','0','1'}

set(gca,'fontsize',40);

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 9];
print('Fig_MatrixD_10by10','-dpdf','-r600')

