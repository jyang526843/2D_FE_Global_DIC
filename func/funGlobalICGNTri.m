%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Triangulation FE-based Global DVC ICGN-code     %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2019.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, normOfW, TimeICGN] = funGlobalICGNTri(coordinatesFEM,elementsFEM,Df,f,g,U,alpha,tol,clusterNo)

DIM = 2; NodesPerEle = 3; % Using cubic elements
FEMSize = DIM*size(coordinatesFEM,1);

DfDx = Df.DfDx; DfDy = Df.DfDy;
DfAxis = Df.DfAxis; DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3);

ImgPydUnit = 1; dirichlet = []; neumann = [];

for stepwithinwhile = 1:100
    tic
    disp(['--- Global IC-GN iteration step',num2str(stepwithinwhile),' ---']);
    
    if (stepwithinwhile==1)
        INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = [];
        %A = sparse(2*size(coordinatesFEM,1),2*size(coordinatesFEM,1));
    end
    INDEXBI = []; INDEXBVAL = []; %clear b; b = sparse(2*size(coordinatesFEM,1),1);
    
    % ================= Navier-Lame elasticity regularization ============
    % MatrixGrad = [-1 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 -1; 0 -1 0 0 1 0 0 0; 0 0 1 0 -1 0 0 0;
    %                0 0 1 0 0 -1 0 0; 0 0 0 -1 0 1 0 0; 0 0 0 1 0 0 -1 0; -1 0 0 0 0 0 1 0];
    %
    % MatrixGradUpdate = MatrixGrad(:,1:4);
    % MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,7) + 0.5*MatrixGrad(:,8);
    % MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,8) + 0.5*MatrixGrad(:,5);
    % MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,5) + 0.5*MatrixGrad(:,6);
    % MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,6) + 0.5*MatrixGrad(:,7);
    %
    % Matrix1 = MatrixGradUpdate'*MatrixGradUpdate;
    % Matrix1Update = zeros(8,8);
    % Matrix1Update(1:2:end,1:2:end) = Matrix1;
    % Matrix1Update(2:2:end,2:2:end) = Matrix1;
    %
    % Matrix2 = zeros(16,8);
    % Matrix2(1:2:end,1:2:end) = MatrixGradUpdate;
    % Matrix2(2:2:end,2:2:end) = MatrixGradUpdate;
    % Matrix2Update = 0.25*Matrix2'*diag([1 0 1 0 0 1 0 1 1 0 1 0 0 1 0 1])*Matrix2;
    %
    % % ------- Lame elasticity constants -------
    % mu = alpha*1; lamda = mu;
    % ====================================================================
    
    hbar = waitbar(0,['Global ICGN iteartion step: ',num2str(stepwithinwhile)]);
    % ============= Each element, assemble stiffness matrix ============
    for j = 1: size(elementsFEM,1) % j is the element index
        
        waitbar(j/size(elementsFEM,1));
        
        tempA = zeros(DIM*NodesPerEle,DIM*NodesPerEle); tempb = tempA(:,1);
        
        % ------ Find three corner pts ------
        pt1x = coordinatesFEM(elementsFEM(j,1),1); pt1y = coordinatesFEM(elementsFEM(j,1),2);
        pt2x = coordinatesFEM(elementsFEM(j,2),1); pt2y = coordinatesFEM(elementsFEM(j,2),2);
        pt3x = coordinatesFEM(elementsFEM(j,3),1); pt3y = coordinatesFEM(elementsFEM(j,3),2);
        
        % ------ Compute triangle area --------
        TriArea = det([1 pt1x pt1y; 1 pt2x pt2y; 1 pt3x pt3y]);
        
        % ------ Calculate DN Matrix for CST ------
        funDN1x = 1/(2*TriArea)*(pt2y-pt3y); funDN1y = 1/(2*TriArea)*(pt3x-pt2x);
        funDN2x = 1/(2*TriArea)*(pt3y-pt1y); funDN2y = 1/(2*TriArea)*(pt1x-pt3x);
        funDN3x = 1/(2*TriArea)*(pt1y-pt2y); funDN3y = 1/(2*TriArea)*(pt2x-pt1x);
        
        DN = [ funDN1x 0 funDN2x 0 funDN3x 0  ;
            funDN1y 0 funDN2y 0 funDN3y 0  ;
            0 funDN1x 0 funDN2x 0 funDN3x  ;
            0 funDN1y 0 funDN2y 0 funDN3y  ];
        
        % ------ Find the element nodal indices ------
        tempIndexU = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2) 2*elementsFEM(j,3)-1 2*elementsFEM(j,3) ];
        
        [ptOfxAll,ptOfyAll] = ndgrid(min([pt1x,pt2x,pt3x]):max([pt1x,pt2x,pt3x]), ...
                                     min([pt1y,pt2y,pt3y]):max([pt1y,pt2y,pt3y]));
        ptOfxAll = ptOfxAll(:); ptOfyAll = ptOfyAll(:); 
        
        for tempjj = 1:length(ptOfxAll)
        
        %for ptOfx = min([pt1x,pt2x,pt3x]):max([pt1x,pt2x,pt3x])
        %    for ptOfy = min([pt1y,pt2y,pt3y]):max([pt1y,pt2y,pt3y])
                
                ptOfx = ptOfxAll(tempjj); ptOfy = ptOfyAll(tempjj); 
                
                % Judge pt is inside triangle or not
                ptInTriangleOrNot = funPointInTriangleCheck(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,ptOfx,ptOfy);
                
                if ptInTriangleOrNot == 1
                    
                    % ------ Calculate N ------
                    N1 = det([1 ptOfx ptOfy; 1 pt2x pt2y; 1 pt3x pt3y])/TriArea;
                    N2 = det([1 ptOfx ptOfy; 1 pt3x pt3y; 1 pt1x pt1y])/TriArea;
                    N3 = det([1 ptOfx ptOfy; 1 pt1x pt1y; 1 pt2x pt2y])/TriArea;
                    
                    N = [N1 0 N2 0 N3 0;
                        0 N1 0 N2 0 N3];
                    
                    % ------ Here approximate Dg(x+u)=Df(x) ------
                    DfEle = [DfDx(ptOfx-DfDxStartx, ptOfy-DfDxStarty);
                        DfDy(ptOfx-DfDxStartx, ptOfy-DfDxStarty)];
                    
                    % ------ Only assemble stiffness in the first step ------
                    if (stepwithinwhile==1)
                        %A(temp,temp) =  A(temp,temp) + (N'*DfEle)*((N'*DfEle)') + alpha*(DN')*DN ;
                        tempA = tempA + (N'*DfEle)*((N'*DfEle)') + alpha*(DN')*DN ;
                    end
                    
                    temp1 = [ptOfx;ptOfy] + N*U(tempIndexU);
                    temp2 = ((f(ptOfx,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), g(floor(temp1(1))-1*ImgPydUnit:ImgPydUnit:floor(temp1(1))+2*ImgPydUnit, floor(temp1(2))-1*ImgPydUnit:ImgPydUnit:floor(temp1(2))+2*ImgPydUnit))) * (N'*DfEle));
                    
                    %b(temp) = b(temp) + tempb - (alpha*(DN')*DN)*U(tempIndexU);
                    tempb = tempb + temp2 - (alpha*(DN')*DN)*U(tempIndexU);
                    
                end
                
         %   end
        end
        
        
        % --- To store A_ele for each element ---
        if (stepwithinwhile==1)
            %  A_ele(j,1:36)=reshape(A(temp,temp),1,36);
            [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
            INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)]; %INDEXAREG = [INDEXAREG;tempAreg(:)];
        end
        INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
        
    end
    
    close(hbar);
    
    if (stepwithinwhile==1)
        % A = sparse(DIM*FEMSize+NodesPerEle*DIM, DIM*FEMSize+NodesPerEle*DIM);
        % b = sparse(DIM*FEMSize+NodesPerEle*DIM,1);
        A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
        A = A+1e-3*max(diag(A))*speye(FEMSize);
       %Areg = sparse( INDEXAI, INDEXAJ, INDEXAREG ,FEMSize,FEMSize);
    end
    b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
    
    
    
    coordsIndexInvolved = unique(elementsFEM);
    if coordsIndexInvolved(1) == 0
        UIndexInvolved = [coordsIndexInvolved(2:end);coordsIndexInvolved(2:end)];
        % Not including the first 0-th entry
        for tempi = 1:(size(coordsIndexInvolved,1)-1)
            UIndexInvolved(2*tempi-1:2*tempi) = [2*coordsIndexInvolved(tempi+1)-1; 2*coordsIndexInvolved(tempi+1)];
        end
    else
        UIndexInvolved = [2*coordsIndexInvolved-1;2*coordsIndexInvolved];
        
    end
    
    W = sparse(2*size(coordinatesFEM,1),1);
    W(2*unique(dirichlet)) = 0;
    W(2*unique(dirichlet)-1) = 0;
    
    dirichlettemp = [2*dirichlet; 2*dirichlet-1];
    FreeNodes = setdiff(UIndexInvolved,unique(dirichlettemp));
    
    W(FreeNodes) = A(FreeNodes,FreeNodes)\b(FreeNodes);
    
    normW = norm(full(W))/sqrt(size(W,1))
    normOfW(stepwithinwhile) = normW;
    TimeICGN(stepwithinwhile) = toc;
    toc
    U = reshape(U,length(U),1); W = reshape(W,length(W),1);
    
    if stepwithinwhile == 1
        normWOld = normW*10;
    else
        normWOld = normOfW(stepwithinwhile-1);
    end
    
    if (normW < tol) || ((normW/normWOld > 1-tol) && (normW/normWOld < 1))
        U = U + W;
        break;
    elseif (normW >= tol && normW < 1/tol)
        U = U + W;
    else
        warning('Get diverged in Global_ICGN!!!')
        break;
    end
    
    
end

end




