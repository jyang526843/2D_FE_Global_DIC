function show(elements3,elements4,coordinates,u)

trisurf(elements3,coordinates(:,1),coordinates(:,2),u',...
'facecolor','interp' )
hold on
% trisurf(elements4,coordinates(:,1),coordinates(:,2),u',...
% 'facecolor','interp' )

Sqx = zeros(4,size(elements4,1)); Sqy = zeros(4,size(elements4,1)); Sqc = zeros(4,size(elements4,1));
for j = 1:size(elements4,1)
    Sqx(1:4,j) = coordinates(elements4(j,1:4),1);
    Sqy(1:4,j) = coordinates(elements4(j,1:4),2);
    Sqc(1:4,j) = u(elements4(j,1:4));
end
patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none'); 

hold off
view(10,40); 
% title('Solution of the Problem')
box on