% ==================================
% Compute strain
% ----------------------------------
% switch MethodToComputeStrain
%   case 0: finite element Gauss points;
%   case 1: central finite difference;
%   case 2: plane fitting method;
% ==================================

switch DICpara.MethodToComputeStrain
        
    case 1
        % Compute strain method I: Use Finite difference operator or FEM solver
        %if Subpb2FDOrFEM == 1 %FD
			D = funDerivativeOp(M,N,DICpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
			FStrain = D*reshape(ULocal,length(ULocal),1); % JY!!! change intrisic coords to world coords
        %else %FEM
		%	GaussPtOrder = 2; [FStrain,~,~] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder);
        %end
        % Plotstrain(FLocal,x,y,f,g);
        % End of computing strain method I.
        
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
         
    case 2
        D = funDerivativeOp(M,N,DICpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
        FStrain = D*reshape(ULocal,length(ULocal),1); % JY!!! change intrisic coords to world coords
        
        % Compute strain method II: Use Plane Fitting method
        prompt = 'What is your half window size: ';
        Rad = input(prompt);
        
        [Uytemp, Uxtemp, UNewtemp] = PlaneFit2(reshape(ULocal(1:2:end),M,N), DICpara.winstepsize, DICpara.winstepsize,Rad);
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Uytemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Uxtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),UNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        [Vytemp, Vxtemp, VNewtemp] = PlaneFit2(reshape(ULocal(2:2:end),M,N), DICpara.winstepsize, DICpara.winstepsize,Rad);
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Vytemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),Vxtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        % figure; mesh(x((Rad+1):M-Rad,(Rad+1):N-Rad),y((Rad+1):M-Rad,(Rad+1):N-Rad),VNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad))
        
        FStraintemp = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        FStraintemp(1:4:end) = reshape(Uxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(2:4:end) = reshape(Vxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(3:4:end) = reshape(Uytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        FStraintemp(4:4:end) = reshape(Vytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        
        % FStraintemp = -FStraintemp; % JY!!! change intrisic coords to world coords
        
        % ULocaltemp = zeros(2*(M-2*Rad)*(N-2*Rad),1);
        % ULocaltemp(1:2:end) = reshape(UNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        % ULocaltemp(2:2:end) = reshape(VNewtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
        
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
        temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
        
        % Plotstrain0(FStraintemp,x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad),f,g);
        
        temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
        for i = 1:(M-2*Rad)*(N-2*Rad)
            temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
        end
        
        FStrain(temp3) = FStraintemp;
        
    case 0
        
        % Compute strain method I: Use FEM solver
        GaussPtOrder = 2; [FStrain] = funGlobal_NodalStrainAvg(DICmesh,USubpb2,GaussPtOrder);
         
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
        
    otherwise
        disp('Wrong Input to compute strain field!')
        
end


%% Update infinitesimal strain to large deformation gradient tensor
% FStrainSmall = FStrain;
FStrainFinite = FStrain;
for tempi = 1:4:length(FStrain)
    
    % Obtain each component of def grad tensor
    dudx = FStrain(tempi);
    dvdx = FStrain(tempi+1);
    dudy = FStrain(tempi+2);
    dvdy = FStrain(tempi+3);
    
    switch DICpara.StrainType
        case 0 % Infinitesimal stran
            % Do nothing
        case 1 % Eluerian strain
            FStrainFinite(tempi) = 1/(1-dudx)-1;
            FStrainFinite(tempi+3) = 1/(1-dvdy)-1;
            FStrainFinite(tempi+2) = dudy/(1-dvdy);
            FStrainFinite(tempi+1) = dvdx/(1-dudx);
        case 2 % Green-Lagrangian strain: E=(C-I)/2
            FStrainFinite(tempi) = 0.5*(dudx*2-dudx^2-dvdx^2);
            FStrainFinite(tempi+3) = 0.5*(dvdy*2-dudy^2-dvdy^2);
            FStrainFinite(tempi+2) = 0.5*(dudy+dvdx-dudx*dudy-dvdx*dvdy);
            FStrainFinite(tempi+1) = 0.5*(dvdx+dudy-dudy*dudx-dvdy*dvdx);
        case 3  
            disp('Press "Ctrl+C" to modify by yourself.'); pause;
        otherwise
            disp('Wrong strain type!');
    end

end

FStraintemp = FStrainFinite(temp3);
