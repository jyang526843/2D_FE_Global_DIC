%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute strain field in the FE-based global DIC method   %
% Author: Jin Yang                                         %
% Last date modified: 2019.03; 2020.10                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F] = funGlobal_NodalStrainAvg(DICmesh,U,GaussPtOrder)

coordinatesFEM = DICmesh.coordinatesFEM; 
elementsFEM = DICmesh.elementsFEM;
DIM = 2; NodesPerEle = 4; GaussPtOrder = 2;

FStrainAvgTimes = zeros(NodesPerEle*size(coordinatesFEM,1),1);
FStrain = zeros(NodesPerEle*size(coordinatesFEM,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for j = 1:size(elementsFEM,1) % j is the element index
    
    StrainWithinEachElementGausspoint = zeros(4,4);
    
    % ----- Find four corner points -------
    pt1x = coordinatesFEM(elementsFEM(j,1),1);
    pt1y = coordinatesFEM(elementsFEM(j,1),2);
    pt2x = coordinatesFEM(elementsFEM(j,2),1);
    pt2y = coordinatesFEM(elementsFEM(j,2),2);
    pt3x = coordinatesFEM(elementsFEM(j,3),1);
    pt3y = coordinatesFEM(elementsFEM(j,3),2);
    pt4x = coordinatesFEM(elementsFEM(j,4),1);
    pt4y = coordinatesFEM(elementsFEM(j,4),2);
    
    % ------ Find element nodal indices ------
    temp = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2)...
        2*elementsFEM(j,3)-1 2*elementsFEM(j,3) 2*elementsFEM(j,4)-1 2*elementsFEM(j,4)];
    
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
            
            StrainWithinEachElementGausspoint(2*(tempi-1)+tempj,1:4) = DN*U(temp);
            
            
        end
    end
    
    %% Extrapolation to strains at nodal points using Gauss points
    
    MatrixExtrapolation = [1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;
                            -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);
                            1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;
                            -0.5 1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3)];
                        
    % ------ Nodal points strain extrapolation using Gauss points -----
    StrainExxWithinEachElementNodalpoint =  MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,1);
        StrainWithinEachElementGausspoint(2,1);
        StrainWithinEachElementGausspoint(3,1);
        StrainWithinEachElementGausspoint(4,1)];
    
    StrainExyWithinEachElementNodalpoint = MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,2);
        StrainWithinEachElementGausspoint(2,2);
        StrainWithinEachElementGausspoint(3,2);
        StrainWithinEachElementGausspoint(4,2)];
    
    StrainEyxWithinEachElementNodalpoint = MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,3);
        StrainWithinEachElementGausspoint(2,3);
        StrainWithinEachElementGausspoint(3,3);
        StrainWithinEachElementGausspoint(4,3)];
    
    StrainEyyWithinEachElementNodalpoint = MatrixExtrapolation * [StrainWithinEachElementGausspoint(1,4);
        StrainWithinEachElementGausspoint(2,4);
        StrainWithinEachElementGausspoint(3,4);
        StrainWithinEachElementGausspoint(4,4)];
    
    StrainWithinEachElementGausspoint(1,1) = StrainExxWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,1) = StrainExxWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,1) = StrainExxWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,1) = StrainExxWithinEachElementNodalpoint(4);
    
    StrainWithinEachElementGausspoint(1,2) = StrainExyWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,2) = StrainExyWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,2) = StrainExyWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,2) = StrainExyWithinEachElementNodalpoint(4);
    
    StrainWithinEachElementGausspoint(1,3) = StrainEyxWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,3) = StrainEyxWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,3) = StrainEyxWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,3) = StrainEyxWithinEachElementNodalpoint(4);
    
    StrainWithinEachElementGausspoint(1,4) = StrainEyyWithinEachElementNodalpoint(1);
    StrainWithinEachElementGausspoint(2,4) = StrainEyyWithinEachElementNodalpoint(2);
    StrainWithinEachElementGausspoint(3,4) = StrainEyyWithinEachElementNodalpoint(3);
    StrainWithinEachElementGausspoint(4,4) = StrainEyyWithinEachElementNodalpoint(4);
    
    
    % ------ Find the element nodal indices for strain ------
    tempStrainIndex = [4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-2  4*elementsFEM(j,1)-1 4*elementsFEM(j,1)  ...
        4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-2  4*elementsFEM(j,2)-1 4*elementsFEM(j,2) ...
        4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-2  4*elementsFEM(j,3)-1 4*elementsFEM(j,3) ...
        4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-2  4*elementsFEM(j,4)-1 4*elementsFEM(j,4)];
    
    
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


F = FStrain./FStrainAvgTimes;
                       
end


% ========= subroutines for  FEM Q4 shape function derivatives ========

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


