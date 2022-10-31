function [Strain,StrainGaussPt,CoordGaussPt] = funGlobal_NodalStrainT3(coordinates,elements,U)

FStrainAvgTimes = zeros(4*size(coordinates,1),1); FStrain = zeros(4*size(coordinates,1),1);
StrainGaussPt = zeros(size(elements,1),4);  CoordGaussPt = zeros(size(elements,1),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(elements,1)
     
    % ------ Find three corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    
    % ------ Compute triangle area --------
    TriArea =  det([1 point1x point1y; 1 point2x point2y; 1 point3x point3y]);
    
    % ------ Calculate DN Matrix for CST ------
    funDN1x = 1/(2*TriArea)*(point2y-point3y);
    funDN1y = 1/(2*TriArea)*(point3x-point2x);
    funDN2x = 1/(2*TriArea)*(point3y-point1y);
    funDN2y = 1/(2*TriArea)*(point1x-point3x);
    funDN3x = 1/(2*TriArea)*(point1y-point2y);
    funDN3y = 1/(2*TriArea)*(point2x-point1x);
    
    DN = [ funDN1x 0 funDN2x 0 funDN3x 0  ;
        funDN1y 0 funDN2y 0 funDN3y 0  ;
        0 funDN1x 0 funDN2x 0 funDN3x  ;
        0 funDN1y 0 funDN2y 0 funDN3y  ];
      
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) 2*elements(j,3)-1 2*elements(j,3) ];
     
    % ------ Compute the strain of the Guass point ------
    StrainGaussPt(j,1:4) = (DN*U(temp))';
    CoordGaussPt(j,1:2) = 1/3*[(point1x+point2x+point3x), (point1y+point2y+point3y)];
    
    % ------ Find the element nodal indices for strain fields ------
    tempStrainIndex = [4*elements(j,1)-3, 4*elements(j,1)-2, 4*elements(j,1)-1, 4*elements(j,1), ...
                       4*elements(j,2)-3, 4*elements(j,2)-2, 4*elements(j,2)-1, 4*elements(j,2), ...
                       4*elements(j,3)-3, 4*elements(j,3)-2, 4*elements(j,3)-1, 4*elements(j,3)];
                    
    FStrain(tempStrainIndex) = FStrain(tempStrainIndex) + ...
                        [StrainGaussPt(j,1);StrainGaussPt(j,3);StrainGaussPt(j,2);StrainGaussPt(j,4);
                        StrainGaussPt(j,1);StrainGaussPt(j,3);StrainGaussPt(j,2);StrainGaussPt(j,4);
                        StrainGaussPt(j,1);StrainGaussPt(j,3);StrainGaussPt(j,2);StrainGaussPt(j,4)];
                    
    FStrainAvgTimes(tempStrainIndex) = FStrainAvgTimes(tempStrainIndex) + ones(12,1);
    
    
end


Strain = FStrain./FStrainAvgTimes;