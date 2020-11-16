% =========================================================
% function funGlobalICGN to solve FE iterations for global DIC
% ---------------------------------------------------------
%   INPUT: 
%       DIC mesh, DIC image pair, 
%       displacement initial guess, regularizer coefficient
%
%   OUTPUT:
%       U: Solved displacement field;
%       normOfW: FE-based global DIC iteration update norm;
%       TimeICGN: Time cost for each FE-based global DIC iteration;
%
% Author: Jin Yang, jyang526@wisc.edu or aldicdvc@gmail.com
% Date: 2020.10
% =========================================================

function [U, normOfW, TimeICGN] = funGlobalICGN(DICmesh,Df,Img1,Img2,U,alpha,tol)


coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;
NodesPerEle = 4; % Num of nodes for each finite element
DIM = 2; % Problem dimension
winsize = (coordinatesFEM(2,1)-coordinatesFEM(1,1))*ones(1,DIM); % or: DICpara.winsize;
FEMSize = DIM*size(coordinatesFEM,1); % FEM problem size
DfDx = Df.DfDx; DfDy = Df.DfDy; DfAxis = Df.DfAxis; DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3); 

 
%% ============================================================== 
% Optional: construct Navier-Lame elasticity regularizer 
% Ignore this section if you don't apply elasticity regularization

% MatrixGrad = [-1 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 -1; 0 -1 0 0 1 0 0 0; 0 0 1 0 -1 0 0 0;
%                0 0 1 0 0 -1 0 0; 0 0 0 -1 0 1 0 0; 0 0 0 1 0 0 -1 0; -1 0 0 0 0 0 1 0];
% MatrixGradUpdate = MatrixGrad(:,1:4);
% MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,7) + 0.5*MatrixGrad(:,8);
% MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,8) + 0.5*MatrixGrad(:,5);
% MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,5) + 0.5*MatrixGrad(:,6);
% MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,6) + 0.5*MatrixGrad(:,7);
% 
% Matrix1 = MatrixGradUpdate'*MatrixGradUpdate;
% Matrix1Update = zeros(8,8); Matrix1Update(1:2:end,1:2:end)=Matrix1; Matrix1Update(2:2:end,2:2:end)=Matrix1;
%  
% % Lame elasticity constants 
% mu = alpha*1; lamda = mu;
     

%% %%%%%%%%%%%%%% Start FE ICGN iteration %%%%%%%%%%%%%%% 
for stepwithinwhile = 1:100 % Max iteration number is set to be 100 by default
     
    tic;  
    
    % ====== Initialize stiffness matrix at the first iteration ======
    if (stepwithinwhile==1)
       disp(['--- Global IC-GN iterations ---']);
       INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = []; % A = sparse(FEMSize,FEMSize);
    end
    INDEXBI = []; INDEXBVAL = []; % clear b; b = sparse(FEMsize,1); % Update external force vector
      
    
    % ------ Define ksi and eta list ------
    ksiList = -1:2/winsize(1):1; etaList = -1:2/winsize(2):1;
    [ksiMat,etaMat] = ndgrid(ksiList,etaList);
    
    NMat = cell(NodesPerEle,1); % Shape functions
    NMat{1} = 0.25*(1-ksiMat).*(1-etaMat); 
    NMat{2} = 0.25*(1+ksiMat).*(1-etaMat); 
    NMat{3} = 0.25*(1+ksiMat).*(1+etaMat); 
    NMat{4} = 0.25*(1-ksiMat).*(1+etaMat); 
         
    hbar = waitbar(0,['Global ICGN iteration step ',num2str(stepwithinwhile)]);
    % ============= Each element, assemble stiffness matrix ============
    for indEle = 1: size(elementsFEM,1) % indEle is the element index
        
        waitbar(indEle/size(elementsFEM,1));  
        
        tempA = zeros(DIM*NodesPerEle, DIM*NodesPerEle); tempb = tempA(:,1);
         
        % ------ Find four corner points in Q4 FE element ------
        pt1x = coordinatesFEM(elementsFEM(indEle,1),1); pt1y = coordinatesFEM(elementsFEM(indEle,1),2);
        pt2x = coordinatesFEM(elementsFEM(indEle,2),1); pt2y = coordinatesFEM(elementsFEM(indEle,2),2);
        pt3x = coordinatesFEM(elementsFEM(indEle,3),1); pt3y = coordinatesFEM(elementsFEM(indEle,3),2);
        pt4x = coordinatesFEM(elementsFEM(indEle,4),1); pt4y = coordinatesFEM(elementsFEM(indEle,4),2);
         
        % ------ Find the element nodal indices ------
        tp = ones(1,DIM); tempIndexU = DIM*elementsFEM(indEle,[tp,2*tp,3*tp,4*tp]);
        for tempDIM = 1:DIM-1
            tempIndexU(tempDIM:DIM:end) = tempIndexU(tempDIM:DIM:end)-(DIM-tempDIM);
        end % size of tempIndexU: 1*8
        % or using the following lines
        % tempIndexU = [2*elementsFEM(indEle,1)-1 2*elementsFEM(indEle,1) 2*elementsFEM(indEle,2)-1 2*elementsFEM(indEle,2)...
        %               2*elementsFEM(indEle,3)-1 2*elementsFEM(indEle,3) 2*elementsFEM(indEle,4)-1 2*elementsFEM(indEle,4)];
        
        [ptOfxAll, ptOfyAll] = ndgrid(pt1x:pt3x, pt1y:pt3y); % To compute each pixels
        
        tempUMat = zeros(winsize+ones(1,2)); tempVMat = tempUMat;
        for tempk = 1:NodesPerEle % Compute affine deformations for all the pixels within the finite element
            tempUMat = tempUMat + (U(tempIndexU(DIM*tempk-1))*ones(winsize+ones(1,DIM))).*NMat{tempk};
            tempVMat = tempVMat + (U(tempIndexU(DIM*tempk-0))*ones(winsize+ones(1,DIM))).*NMat{tempk};
        end
        
        tempg = ba_interp2(Img2, ptOfyAll+tempVMat, ptOfxAll+tempUMat, 'cubic'); % Deformed g(x+u)
        
        
        ptOfxAll = ptOfxAll(:); ptOfyAll = ptOfyAll(:); % Garantuee ptOfxAll and ptOfyAll are column vectors
        for tempjj = 1:length(ptOfxAll) % Write into one for-loop instead of two for-loops
        % for ptOfx = pt1x:pt3x
        %    for ptOfy = pt1y:pt3y
            ptOfx = ptOfxAll(tempjj); ptOfy = ptOfyAll(tempjj);
              
                % ------ Calculate ksi and eta ------
                ksi = ksiMat(tempjj); eta = etaMat(tempjj);
                % ksi = l(1)*pointOfx*pointOfy + l(2)*pointOfx + l(3)*pointOfy + l(4) ;
                % eta = m(1)*pointOfx*pointOfy + m(2)*pointOfx + m(3)*pointOfy + m(4) ;
                
                % ------ Calculate N matrix ------
                N1 = NMat{1}(tempjj); N2 = NMat{2}(tempjj); 
                N3 = NMat{3}(tempjj); N4 = NMat{4}(tempjj);
                % N1 = (1-ksi)*(1-eta)*0.25; N2 = (1+ksi)*(1-eta)*0.25;
                % N3 = (1+ksi)*(1+eta)*0.25; N4 = (1-ksi)*(1+eta)*0.25;
                
                % ------ Generate [N] shape function matrix ------
                N = [N1 0 N2 0 N3 0 N4 0;
                     0 N1 0 N2 0 N3 0 N4];
                 
                % ------ Build J matrix ------
                J = [funDN1Dksi(ksi,eta),funDN2Dksi(ksi,eta),funDN3Dksi(ksi,eta),funDN4Dksi(ksi,eta);
                    funDN1Deta(ksi,eta),funDN2Deta(ksi,eta),funDN3Deta(ksi,eta),funDN4Deta(ksi,eta)] * ...
                    [pt1x,pt1y;pt2x,pt2y;pt3x,pt3y;pt4x,pt4y];
                % J11 = funDN1Dksi(ksi,eta)*pt1x + funDN2Dksi(ksi,eta)*pt2x + ...
                %       funDN3Dksi(ksi,eta)*pt3x + funDN4Dksi(ksi,eta)*pt4x;
                % J12 = funDN1Dksi(ksi,eta)*pt1y + funDN2Dksi(ksi,eta)*pt2y + ...
                %       funDN3Dksi(ksi,eta)*pt3y + funDN4Dksi(ksi,eta)*pt4y;
                % J21 = funDN1Deta(ksi,eta)*pt1x + funDN2Deta(ksi,eta)*pt2x + ...
                %       funDN3Deta(ksi,eta)*pt3x + funDN4Deta(ksi,eta)*pt4x;
                % J22 = funDN1Deta(ksi,eta)*pt1y + funDN2Deta(ksi,eta)*pt2y + ...
                %       funDN3Deta(ksi,eta)*pt3y + funDN4Deta(ksi,eta)*pt4y;
                % J = [J11 J12; J21 J22];
                
                Jacobian = det(J);
                InvJ = 1/Jacobian*[J(2,2) -J(1,2); -J(2,1) J(1,1)]; % Inverse of J matrix
                
                
                % ------ Compute DN matrix ------
                DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                    [funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta) 0;
                     funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta) 0;
                     0 funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta);
                     0 funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta)];
                
                 
                % ------ Here approximate Dg(x+u)=Df(x) ------
                DfEle = [DfDx(ptOfx-DfDxStartx, ptOfy-DfDxStarty); 
                         DfDy(ptOfx-DfDxStartx, ptOfy-DfDxStarty)];
                
                % ------ Only assemble stiffness in the first step ------
                if (stepwithinwhile==1)
                    %A(tempIndexU,tempIndexU) =  A(tempIndexU,tempIndexU) + (N'*Df)*(N'*Df)' + alpha*(DN')*DN ;
                    tempA = tempA + (N'*DfEle)*(N'*DfEle)' + alpha*(DN')*DN;
                end
                
                % ------ Construct b vector ------
                %temp1 = [pointOfx;pointOfy] + N*U(temp);
                %tempb = ((f(pointOfx ,pointOfy) - fungInterpolation_g(temp1(1), temp1(2), g(floor(temp1(1))-1:floor(temp1(1))+2, floor(temp1(2))-1:floor(temp1(2))+2))) * (N'*Df));
                %b(tempIndexU) = b(tempIndexU) + tempb - (alpha*(DN')*DN)*U(temp);
                temp2 = (Img1(ptOfx,ptOfy) - tempg(tempjj)) * (N'*DfEle);
                tempb = tempb + temp2 - (alpha*(DN')*DN) * U(tempIndexU);
                
                % ====== Optional ======
                % --- To use Navier-Lame elasticity regularization instead ---
                % b(temp) = b(temp) + tempb + (mu*Matrix1Update + (mu+lamda)*Matrix2Update)*U(temp);
                % ------------------------------------------------------------
                 
            %end % for ptOfx = pt1x:pt3x
        % end % for ptOfy = pt1y:pt3y
        end
          
        % --- To store A_ele for each element ---
        if (stepwithinwhile==1)
            [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
            INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)];
        end
        INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
        
    end
    
    close(hbar);
    
    if (stepwithinwhile==1)
        A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
    end
    b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
    
    % ========= Solve FEM problem ===========
    W = A\b;
    
    normW = norm(W)/sqrt(size(W,1)); 
    normOfW(stepwithinwhile) = normW;
    TimeICGN(stepwithinwhile) = toc; 
    U = reshape(U,length(U),1); W = reshape(W,length(W),1);
    
    disp(['normW = ',num2str(normW),' at iter ',num2str(stepwithinwhile),'; time cost = ',num2str(toc),'s']);
     
    if stepwithinwhile == 1
        normWOld = normW*10;
    else
        normWOld = normOfW(stepwithinwhile-1);
    end
    
    if (normW < tol) % || ((normW/normWOld > 0.9) && (normW/normWOld < 1))
        U = U + W; 
        break;
    elseif (normW >= tol && normW < (0.1/tol)) % || ((normW/normWOld >= 1) && (normW/normWOld < 100)))
        U = U + W;
    else
        warning('Get diverged in Global_ICGN!!!')
        break;
    end
    
    
end 
   
TotalTimeICGN = sum(TimeICGN);
disp(['Elapsed time is ',num2str(TotalTimeICGN),' seconds.']);

end


%% ========= subroutines for  FEM Q4 shape function derivatives ========
function DN1x=funDN1Dksi(ksi,eta)
DN1x = -(1-eta)/4 ;
end
function DN1y=funDN1Deta(ksi,eta)
DN1y =  -(1-ksi)/4 ;
end
function DN2x=funDN2Dksi(ksi,eta)
DN2x =  (1-eta)/4 ;
end
function DN2y=funDN2Deta(ksi,eta)
DN2y =  -(1+ksi)/4 ;
end
function DN3x=funDN3Dksi(ksi,eta)
DN3x = (1+eta)/4 ;
end
function DN3y=funDN3Deta(ksi,eta)
DN3y =  (1+ksi)/4 ;
end
function DN4x=funDN4Dksi(ksi,eta)
DN4x = -(1+eta)/4 ;
end
function DN4y=funDN4Deta(ksi,eta)
DN4y = (1-ksi)/4 ;
end


