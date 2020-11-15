%close all;

clear u_ dm m2vx
L = 9;
u_ = cell(1);
temp1 = reshape(ULocal(1:2:end),M,N); temp1 = repmat(temp1,1,1,L); temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{1}{1} = temp10;
temp1 = reshape(ULocal(2:2:end),M,N); temp1 = repmat(temp1,1,1,L); temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{2}{1} = temp10;
 
u_{3}{1} = 0*temp10;

temp1 = reshape(sqrt((ULocal(1:2:end).^2 + ULocal(2:2:end).^2  )),M,N );  temp1 = repmat(temp1,1,1,L); temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{4}{1} = temp10;
 
dm  = [winstepsize];

m2vx  = [1,1,1];

[Fij, J] = calculateFij(u_,dm,m2vx,'optimal9'); 

%% Calculate Lagrangian strain (Eij), eulerian strain (eij), Strain energy density (U), and stress (Sij)
[Eij, eij] = calculateEij(Fij);              
%[Sij, U] = calculateSij(Fij, J, eij{1,1}, Eij{1,1}, [E v], mechModel); 
%save('mechanicsvariables_dyedcell_xy5.mat', 'Fij','J','Eij','eij','U','Sij');

%%
FStrain = zeros(M*N*4,1);
FStrain(1:4:end) =  reshape(Eij{1}{5}(:,:,1),M*N ,1);
FStrain(4:4:end) = reshape(Eij{1}{1}(:,:,1),M*N ,1);
% FSubpb3(9:9:end) = reshape(Fij{1}{9},M*N*L,1)-1;
FStrain(2:4:end) = reshape(Eij{1}{4}(:,:,1),M*N ,1);
FStrain(3:4:end) = reshape(Eij{1}{2}(:,:,1),M*N ,1);
% FSubpb3(3:9:end) = reshape(Fij{1}{6},M*N*L,1);
% FSubpb3(7:9:end) = reshape(Fij{1}{8},M*N*L,1);
% FSubpb3(6:9:end) = reshape(Fij{1}{3},M*N*L,1);
% FSubpb3(8:9:end) = reshape(Fij{1}{7},M*N*L,1);
%for index=[1:9]
%FSubpb3(index:9:end) = reshape(Fij{1}{index},M*N*L,1);
%end 
%Plotstrain_show3(FSubpb3,ResultcoordinatesFEM{1}.coordinatesFEM,ResultelementsFEM{1}.elementsFEM);


%%
Rad = 1; % FilterSizeInput-2;
% Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
temp = 1:1:size(coordinatesFEM,1); temp = temp';
temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);

temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
for i = 1:(M-2*Rad)*(N-2*Rad)
    temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
end
FStraintemp = FStrain(temp3);


