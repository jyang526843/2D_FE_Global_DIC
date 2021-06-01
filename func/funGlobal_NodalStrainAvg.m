function [StrainNodalPt] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,U,GaussPtOrder)
%FUNCTION [StrainNodalPt] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,U,GaussPtOrder))                 
% Object: to compute strain fields by the FE-method
% ----------------------------------------------
%
%	INPUT: coordinatesFEM      DIC FE Q4 mesh coordinates 
%          elementsFEM         DIC FE Q4 mesh elements 
%          U                   Disp vector: U = [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%          GaussPtOrder        Gauss point order used in FE Q4 element,  GaussPtOrder = 2 BY DEFAULT
%            
%   OUTPUT: StrainNodalPt      Solved strains at the nodal points,
%   
%            
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2018.03, 2020.12 
% ==============================================


%% Initialization
FStrainAvgTimes = zeros(4*size(coordinatesFEM,1),1);
FStrain = zeros(4*size(coordinatesFEM,1),1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for eleInd = 1:size(elementsFEM,1) % eleInd is the element index
    
    StrainWithinEachElementGausspoint = zeros(4,4);
    
    % ----- Find four corner points -------
    pt1x = coordinatesFEM(elementsFEM(eleInd,1),1);
    pt1y = coordinatesFEM(elementsFEM(eleInd,1),2);
    pt2x = coordinatesFEM(elementsFEM(eleInd,2),1);
    pt2y = coordinatesFEM(elementsFEM(eleInd,2),2);
    pt3x = coordinatesFEM(elementsFEM(eleInd,3),1);
    pt3y = coordinatesFEM(elementsFEM(eleInd,3),2);
    pt4x = coordinatesFEM(elementsFEM(eleInd,4),1);
    pt4y = coordinatesFEM(elementsFEM(eleInd,4),2);
    
    % ------ Find the element nodal indices ------
    tempIndexU = 2*elementsFEM(eleInd,[1,1,2,2,3,3,4,4]);
    tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;
    
    % ------ Four Gauss points ------
    for tempi = 1:2
        for tempj = 1:2
            
            %ksi = ( 2*tempj-3)/sqrt(3); eta = (2*tempi-3)/sqrt(3);
            ksi = 2*tempj-3; eta = 2*tempi-3;
            if (tempi == 1) && (tempj == 1)
                ksi = -1/sqrt(3); eta = -1/sqrt(3);
            elseif (tempi == 1) && (tempj == 2)
                ksi = 1/sqrt(3); eta = -1/sqrt(3);
            elseif (tempi == 2) && (tempj == 1)
                ksi = 1/sqrt(3); eta = 1/sqrt(3);
            elseif (tempi == 2) && (tempj == 2)
                ksi = -1/sqrt(3); eta = 1/sqrt(3);
            end
            
            
            % ------ Calculate N ------
            N1 = (1-ksi)*(1-eta)*0.25;
            N2 = (1+ksi)*(1-eta)*0.25;
            N3 = (1+ksi)*(1+eta)*0.25;
            N4 = (1-ksi)*(1+eta)*0.25;
            N = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            
            
            % ------ Build J matrix ------
            J11 = funDN1Dksi(ksi,eta)*pt1x + funDN2Dksi(ksi,eta)*pt2x + ...
                funDN3Dksi(ksi,eta)*pt3x + funDN4Dksi(ksi,eta)*pt4x;
            J12 = funDN1Dksi(ksi,eta)*pt1y + funDN2Dksi(ksi,eta)*pt2y + ...
                funDN3Dksi(ksi,eta)*pt3y + funDN4Dksi(ksi,eta)*pt4y;
            J21 = funDN1Deta(ksi,eta)*pt1x + funDN2Deta(ksi,eta)*pt2x + ...
                funDN3Deta(ksi,eta)*pt3x + funDN4Deta(ksi,eta)*pt4x;
            J22 = funDN1Deta(ksi,eta)*pt1y + funDN2Deta(ksi,eta)*pt2y + ...
                funDN3Deta(ksi,eta)*pt3y + funDN4Deta(ksi,eta)*pt4y;
            J = [J11 J12; J21 J22];
            
            Jacobian = det(J);
            InvJ = 1/Jacobian*[J22 -J12; -J21 J11];
            
            % ------ Compute DN matrix ------
            DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                [funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta) 0;
                funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta) 0;
                0 funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta);
                0 funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta)];
            
            StrainWithinEachElementGausspoint(2*(tempi-1)+tempj,1:4) = DN*U(tempIndexU);
            
            
        end
    end
    
    
    % ------ Extrapolation matrix ------
    MatrixExtrapolation = [ 1+0.5*sqrt(3),   -0.5,            1-0.5*sqrt(3),   -0.5;
                            -0.5,            1+0.5*sqrt(3),   -0.5,            1-0.5*sqrt(3);
                            1-0.5*sqrt(3),   -0.5,            1+0.5*sqrt(3),   -0.5;
                            -0.5,            1-0.5*sqrt(3),   -0.5,            1+0.5*sqrt(3)];
    
    % ------ Nodal points strain extrapolation using Gauss points -----
    StrainExxWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(1:4,1);
    StrainExyWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(1:4,2); 
    StrainEyxWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(1:4,3);
    StrainEyyWithinEachElementNodalpoint = MatrixExtrapolation * StrainWithinEachElementGausspoint(1:4,4);
             
    StrainWithinEachElementGausspoint(1:4,1) = reshape(StrainExxWithinEachElementNodalpoint(1:4),4,1);
    StrainWithinEachElementGausspoint(1:4,2) = reshape(StrainExyWithinEachElementNodalpoint(1:4),4,1);
    StrainWithinEachElementGausspoint(1:4,3) = reshape(StrainEyxWithinEachElementNodalpoint(1:4),4,1);
    StrainWithinEachElementGausspoint(1:4,4) = reshape(StrainEyyWithinEachElementNodalpoint(1:4),4,1);
     
    
    % ------ Find the element nodal indices for strain ------
    tempStrainIndex = [ 4*elementsFEM(eleInd,1)-3 4*elementsFEM(eleInd,1)-2  4*elementsFEM(eleInd,1)-1 4*elementsFEM(eleInd,1)  ...
                        4*elementsFEM(eleInd,2)-3 4*elementsFEM(eleInd,2)-2  4*elementsFEM(eleInd,2)-1 4*elementsFEM(eleInd,2) ...
                        4*elementsFEM(eleInd,3)-3 4*elementsFEM(eleInd,3)-2  4*elementsFEM(eleInd,3)-1 4*elementsFEM(eleInd,3) ...
                        4*elementsFEM(eleInd,4)-3 4*elementsFEM(eleInd,4)-2  4*elementsFEM(eleInd,4)-1 4*elementsFEM(eleInd,4)];
    
    
    FStrain(tempStrainIndex) = FStrain(tempStrainIndex)+[StrainWithinEachElementGausspoint(1,1);
        StrainWithinEachElementGausspoint(1,2); StrainWithinEachElementGausspoint(1,3);
        StrainWithinEachElementGausspoint(1,4); StrainWithinEachElementGausspoint(2,1);
        StrainWithinEachElementGausspoint(2,2); StrainWithinEachElementGausspoint(2,3);
        StrainWithinEachElementGausspoint(2,4); StrainWithinEachElementGausspoint(3,1);
        StrainWithinEachElementGausspoint(3,2); StrainWithinEachElementGausspoint(3,3);
        StrainWithinEachElementGausspoint(3,4); StrainWithinEachElementGausspoint(4,1);
        StrainWithinEachElementGausspoint(4,2); StrainWithinEachElementGausspoint(4,3);
        StrainWithinEachElementGausspoint(4,4)];
    
    FStrainAvgTimes(tempStrainIndex) = FStrainAvgTimes(tempStrainIndex) + ones(16,1);
    
     
end
 
StrainNodalPt = FStrain./FStrainAvgTimes;


end
            


%% ========= subroutines for FEM Q4 shape function derivatives ========
function DN1Dksi=funDN1Dksi(ksi,eta)
    DN1Dksi = -(1-eta)/4 ;
end
function DN1Deta=funDN1Deta(ksi,eta)
    DN1Deta =  -(1-ksi)/4 ;
end
function DN2Dksi=funDN2Dksi(ksi,eta)
    DN2Dksi =  (1-eta)/4 ;
end
function DN2Deta=funDN2Deta(ksi,eta)
    DN2Deta =  -(1+ksi)/4 ;
end
function DN3Dksi=funDN3Dksi(ksi,eta)
    DN3Dksi = (1+eta)/4 ;
end
function DN3Deta=funDN3Deta(ksi,eta)
    DN3Deta =  (1+ksi)/4 ;
end
function DN4Dksi=funDN4Dksi(ksi,eta)
    DN4Dksi = -(1+eta)/4 ;
end
function DN4Deta=funDN4Deta(ksi,eta)
    DN4Deta = (1-ksi)/4 ;
end

