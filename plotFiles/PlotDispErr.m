function [xAxis,avgDispxx,stdDispxx,ExactDispxx]=PlotDispErr(coordinatesFEM,USubpb2,x,y,M,N, Exact)
  
ErrDisp = zeros(size(coordinatesFEM,1),1); avgDispErr = 0;
for i = 1:size(coordinatesFEM,1)
   ErrDisp(i) = sqrt( (USubpb2(2*i-1) + Exact(coordinatesFEM(i,1),6))^2  + (USubpb2(2*i) - 0)^2 );
   avgDispErr = avgDispErr + ErrDisp(i)^2; 
end
avgDispErr = sqrt(avgDispErr/size(coordinatesFEM,1))
% figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
% figure; mesh(x,y,reshape(ErrDisp,M,N)); axis tight; set(gca,'fontSize',18); view(-20, 50); title('Displacement Absolute Error')


Dispxx = zeros(size(coordinatesFEM,1),1); avgDispxx = 0; stdDispxx = 0; xAxis = 0; ExactDispxx = 0;
for i = 1:size(coordinatesFEM,1)
    Dispxx(i) =  USubpb2(2*i-1)  ;
end
Dispxx = reshape(Dispxx,M,N);
for i = 1:M
    avgDispxx(i) = sum(Dispxx(i,:))/N;
    stdDispxx(i) = std(Dispxx(i,:));
    ExactDispxx(i) = - Exact(x(i,1),6) ;
end
avgDispxx=avgDispxx';
stdDispxx=stdDispxx';
ExactDispxx=ExactDispxx';
xAxis = x(:,1);