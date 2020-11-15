function [xAxis,avgStrainxx,stdStrainxx,ExactStrainxx]= PlotStrainErr(coordinatesFEM,FSubpb2,x,y,M,N, Exact)
 
ErrStrain = zeros(size(coordinatesFEM,1),1); avgErrStrain = 0;
for i = 1:size(coordinatesFEM,1)
   ErrStrain(i) = sqrt((FSubpb2(4*i-3) + 1e-6*Exact(coordinatesFEM(i,1),7))^2  + (FSubpb2(4*i-2) - 0)^2 + (FSubpb2(4*i-1) - 0)^2 + (FSubpb2(4*i) - 0)^2 );
   avgErrStrain = avgErrStrain + ErrStrain(i)^2; 
end
avgErrStrain = sqrt(avgErrStrain/size(coordinatesFEM,1))
% figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
% figure; mesh(x,y,reshape(ErrStrain,M,N)); axis tight; set(gca,'fontSize',20); view(-20, 50); title('Strain Absolute Error')
  
  Strainxx = zeros(size(coordinatesFEM,1),1); avgStrainxx = 0; stdStrainxx = 0; xAxis = 0; ExactStrainxx = 0;
  for i = 1:size(coordinatesFEM,1)
  Strainxx(i) =  FSubpb2(4*i-3)  ;
  end
  Strainxx = reshape(Strainxx,M,N); 
  for i = 1:M
      avgStrainxx(i) = sum(Strainxx(i,:))/N;
      stdStrainxx(i) = std(Strainxx(i,:));
      ExactStrainxx(i) = -1e-6*Exact(x(i,1),7) ;
  end
  avgStrainxx=avgStrainxx';
  stdStrainxx=stdStrainxx';
    ExactStrainxx=ExactStrainxx';
    xAxis = x(:,1);  
  %figure; errorbar(xAxis, avgStrainxx  , stdStrainxx );