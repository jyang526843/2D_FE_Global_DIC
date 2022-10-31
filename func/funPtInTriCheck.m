function ptInTriOrNot = funPtInTriCheck(TriNodeCoord,PtToCheckCoord)
%FUNPTINTRICHECK: To check points inside or outside a given 2D triangle 
%   ptInTriOrNot = funPtInTriCheck(TriNodeCoord,PtToCheckCoord)
%   
%   INPUT: 
%       Triangle node coordinates:  TriNodeCoord: 3*2 matrix
%                                   [coordx_triangle_pt1, coordy_triangle_pt1;
%                                    coordx_triangle_pt2, coordy_triangle_pt2;
%                                    coordx_triangle_pt3, coordy_triangle_pt3];
%
%       Coordinates of the checking points:  PtToCheckCoord: N*2 matrix with N points to check
%                                   [coordx_pti  coordy_pti] = PtToCheckCoord(i,1:2)
%
%   OUTPUT:
%       Logical variable: ptInTriOrNot: N*1 vector
%                         1 --> Yes, point is inside the given triangle
%                         0 --> No,  point is outside the given triangle
%
% -----------------------------------------------
% Author: Jin Yang (jyang526@wisc.edu)
% Date: 07-07-2020
%
% Reference
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point1x = TriNodeCoord(1,1); point1y = TriNodeCoord(1,2); 
point2x = TriNodeCoord(2,1); point2y = TriNodeCoord(2,2); 
point3x = TriNodeCoord(3,1); point3y = TriNodeCoord(3,2); 

pointOfx = PtToCheckCoord(:,1); pointOfy = PtToCheckCoord(:,2);


%%
p1x = pointOfx; p1y = pointOfy; p2x = point1x; p2y = point1y; p3x = point2x; p3y = point2y;
b1 = sign( (p1x-p3x)*(p2y-p3y) - (p2x-p3x)*(p1y-p3y) );

p1x = pointOfx; p1y = pointOfy; p2x = point2x; p2y = point2y; p3x = point3x; p3y = point3y;
b2 = sign( (p1x-p3x)*(p2y-p3y) - (p2x-p3x)*(p1y-p3y) );

p1x = pointOfx; p1y = pointOfy; p2x = point3x; p2y = point3y; p3x = point1x; p3y = point1y;
b3 = sign( (p1x-p3x)*(p2y-p3y) - (p2x-p3x)*(p1y-p3y) );


ptInTriOrNot = 1-logical(sign(abs(b1-b2)+abs(b2-b3)));
% Comment: The above line is to run the following codes
% if (b1==b2) && (b2==b3)
%    pointInTriangleOrNot = 0; % point is inside of the triangle
% else
%    pointInTriangleOrNot = 1; % point is outside of the triangle
% end



