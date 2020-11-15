


U0World = U0; U0World(2:2:end) = -U0(2:2:end); 
% DICmesh.y0World = (size(fNormalized,2)+1-y0); DICmesh.coordinatesFEMWorld = [DICmesh.coordinatesFEM(:,1),size(fNormalized,2)+1-DICmesh.coordinatesFEM(:,2)];

close all; Plotuv(U0World,DICmesh.x0,DICmesh.y0World); % Plotdisp_show(U,coordinatesFEM,elementsFEM); % Plotuv(UExact,x0,y0);
% figure(1); view([-140,60]);% colormap(coolwarm(32)); % caxis([-1.0,1.0]); % caxis([-0.45,0.45]);
% figure(2); view([-140,60]);% colormap(coolwarm(32)); % caxis([-0.085, 0.005]); % caxis([-0.085,0.015 ]); 
Plotdisp_show(U0World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM); % Plot initial values
figure(3); axis equal;% colormap(coolwarm(32)); colorbar; view(2); axis tight; set(gca,'fontsize',18); box on; % caxis([-1.0,1.0]); % caxis([-0.45,0.45]);
figure(4); axis equal; %colormap(coolwarm(32)); colorbar; view(2); axis tight; set(gca,'fontsize',18); box on; % caxis([-0.085, 0.005]); % caxis([-0.085,0.015 ]); 
