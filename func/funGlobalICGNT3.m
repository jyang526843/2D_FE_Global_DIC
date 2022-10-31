%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Triangulation FE-based Global DVC ICGN-code     %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2019.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, normOfW, TimeICGN] = funGlobalICGNT3(DICmesh,Df,Img1,Img2,U,alpha,tol,maxIter)

coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM; 
clusterNo = 8;

DIM = 2; NodesPerEle = 3; % Used elements
FEMSize = DIM*size(coordinatesFEM,1);

DfDx = Df.DfDx; DfDy = Df.DfDy;
DfAxis = Df.DfAxis; DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3);

ImgPydUnit = 1; dirichlet = []; neumann = [];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for stepwithinwhile = 1:maxIter
    tic
    disp(['--- Global IC-GN iteration step',num2str(stepwithinwhile),' ---']);
     
%     if clusterNo==0 || clusterNo==1
%     if (stepwithinwhile==1)
%         INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = [];
%         %A = sparse(2*size(coordinatesFEM,1),2*size(coordinatesFEM,1));
%     end
%     INDEXBI = []; INDEXBVAL = []; %clear b; b = sparse(2*size(coordinatesFEM,1),1);
%     end
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
    if clusterNo== 0 || clusterNo== 1
        hbar = waitbar(0, ['Global ICGN iteartion step: ',num2str(stepwithinwhile)]);
    else
        hbar=parfor_progressbar(size(elementsFEM,1),['Global ICGN iteartion step: ',num2str(stepwithinwhile)]);
    end
    
    % ============= Each element, assemble stiffness matrix ============
    parfor j = 1: size(elementsFEM,1) % j is the element index
        
        if clusterNo== 0 || clusterNo== 1
            waitbar(j/size(elementsFEM,1));
        else
            hbar.iterate(1);
        end
        
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
        
        DN = [  funDN1x 0 funDN2x 0 funDN3x 0  ;
                funDN1y 0 funDN2y 0 funDN3y 0  ;
                0 funDN1x 0 funDN2x 0 funDN3x  ;
                0 funDN1y 0 funDN2y 0 funDN3y  ];
        
        % ------ Find the element nodal indices ------
        tempIndexU = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2) ...
                                                              2*elementsFEM(j,3)-1 2*elementsFEM(j,3) ];
         
        
       % ------ Combine all the points -------
       [ptOfxAll,ptOfyAll] = ndgrid(floor(min([pt1x,pt2x,pt3x])):ceil(max([pt1x,pt2x,pt3x])), ...
                                    floor(min([pt1y,pt2y,pt3y])):ceil(max([pt1y,pt2y,pt3y])));
                                
       % ------ Shape function N matrix ------
       NMat = cell(3,1); % Three nodes for each triangle element
       NMat{1} = (pt2x*pt3y + ptOfxAll*pt2y + ptOfyAll*pt3x - pt2x*ptOfyAll - pt2y*pt3x - pt3y*ptOfxAll)/TriArea;
       NMat{2} = (pt3x*pt1y + ptOfxAll*pt3y + ptOfyAll*pt1x - pt3x*ptOfyAll - pt3y*pt1x - pt1y*ptOfxAll)/TriArea;
       NMat{3} = (pt1x*pt2y + ptOfxAll*pt1y + ptOfyAll*pt2x - pt1x*ptOfyAll - pt1y*pt2x - pt2y*ptOfxAll)/TriArea;
       
       tempUMat = zeros(size(ptOfxAll)); tempVMat = tempUMat;
       for tempk = 1:NodesPerEle
           tempUMat = tempUMat + (U(tempIndexU(2*tempk-1))*ones(size(ptOfxAll))).*NMat{tempk};
           tempVMat = tempVMat + (U(tempIndexU(2*tempk-0))*ones(size(ptOfxAll))).*NMat{tempk};
       end
       
       tempg = ba_interp2(Img2, ptOfyAll+tempVMat, ptOfxAll+tempUMat, 'cubic');
       
       
       ptOfxAll = ptOfxAll(:); ptOfyAll = ptOfyAll(:);
       ptInTriOrNotMat = funPtInTriCheck([pt1x,pt1y;pt2x,pt2y;pt3x,pt3y],[ptOfxAll,ptOfyAll]);
       
       
       
       % ====== Assemble stiffness matrix and force vector ======
      %  for ptOfx = min([pt1x,pt2x,pt3x]):ImgPydUnit:max([pt1x,pt2x,pt3x])
      %      for ptOfy = min([pt1y,pt2y,pt3y]):ImgPydUnit:max([pt1y,pt2y,pt3y])
      for tempjj = 1:length(ptOfxAll)          
           
                % Judge pt is inside triangle or not
                % ptInTriangleOrNot = ptInTriangleCheck(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,ptOfx,ptOfy);
                
                if ptInTriOrNotMat(tempjj) == 1 % 0-Yes, inside the triangle; 1-No, outside the triangle;
                    
                    ptOfx = ptOfxAll(tempjj); ptOfy = ptOfyAll(tempjj);
                     
                    % ------ Calculate N ------
                    %N1 = det([1 ptOfx ptOfy; 1 pt2x pt2y; 1 pt3x pt3y])/TriArea;
                    %N2 = det([1 ptOfx ptOfy; 1 pt3x pt3y; 1 pt1x pt1y])/TriArea;
                    %N3 = det([1 ptOfx ptOfy; 1 pt1x pt1y; 1 pt2x pt2y])/TriArea;
                    N1 = NMat{1}(tempjj); N2 = NMat{2}(tempjj); N3 = NMat{3}(tempjj);
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
                    
                    % ------ Construct b vector ------
                    %temp1 = [ptOfx;ptOfy] + N*U(tempIndexU);
                    %temp2 = ((Img1(ptOfx,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), Img2(floor(temp1(1))-1*ImgPydUnit:ImgPydUnit:floor(temp1(1))+2*ImgPydUnit, floor(temp1(2))-1*ImgPydUnit:ImgPydUnit:floor(temp1(2))+2*ImgPydUnit))) * (N'*DfEle));
                    temp2 = ((Img1(ptOfx,ptOfy) - tempg(tempjj)) * (N'*DfEle));
                    
                    %b(temp) = b(temp) + tempb - (alpha*(DN')*DN)*U(tempIndexU);
                    tempb = tempb + temp2 - (alpha*(DN')*DN)*U(tempIndexU);
                     
                end
                
            %end
        %end
       end
        
       % --- To store A_ele for each element ---
       if (stepwithinwhile==1)
           %  A_ele(j,1:36)=reshape(A(temp,temp),1,36);
           [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
           % INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)]; %INDEXAREG = [INDEXAREG;tempAreg(:)];
           INDEXAIpar{j}=IndexAXX(:); INDEXAJpar{j}=IndexAYY(:); INDEXAVALpar{j}=tempA(:);
       end
       % INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
       INDEXBIpar{j}=tempIndexU(:); INDEXBVALpar{j}=tempb(:);
        
    end
    
    close(hbar);
    
    if (stepwithinwhile==1)
        A = sparse(FEMSize, FEMSize);
        for eleInd = 1:size(elementsFEM,1)
           A = A + sparse(INDEXAIpar{eleInd}, INDEXAJpar{eleInd}, INDEXAVALpar{eleInd}, FEMSize,FEMSize) ;
        end
        % A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
    end
    
    %b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
    b = sparse(FEMSize,1);
    for eleInd = 1:size(elementsFEM,1)
        b = b + sparse(double(INDEXBIpar{eleInd}),ones(length(INDEXBIpar{eleInd}),1),INDEXBVALpar{eleInd}, FEMSize,1) ;
    end
     
    
    
    % ====== Find involved coordiantes index ======
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
    
    normW = norm(W)/sqrt(size(W,1))
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




