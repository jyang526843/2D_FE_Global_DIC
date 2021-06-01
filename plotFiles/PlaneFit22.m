function [UNew,UX,UY,iGrid,jGrid] = PlaneFit22(U,xstep,ystep,Rad)
%PLANEFIT22: 2D least square fittings
%	[UY,UX,UNew] = PLANEFIT2(U,xstep,ystep,Rad)
%
% Input: U        a gridded data with size (M,N)
%        xstep    x spacing between gridded data
%        ystep    y spacing between gridded data
%        Rad      plane fitting radius
%
% Output: UY      Fitted dudy
%         UX      Fitted dudx
%         UNew    Fitted U
%
% Formula to fit:  u(x,y) = UNew + UX*(x-x0) + UY*(y-y0)
%
% Author: Jin Yang, aldicdvc@gmail.com
% Date: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
[M,N] = size(U); % size of gridded data
 
[iGrid,jGrid] = ndgrid(Rad+1:M-Rad, Rad+1:N-Rad); % Generate computation points
UNew = 0*iGrid; UX = UNew; UY = UNew; % Initalization {UNew,UY,UX}
iGrid = iGrid(:); jGrid = jGrid(:); % Transform into long vectors
[xGrid,yGrid] = ndgrid(-Rad:Rad, -Rad:Rad); 

%%
% hbar = parfor_progressbar(length(iGrid) ,'Please wait for PlaneFit2!');
hbar = waitbar(0,'wait for PlaneFit2!');
for ptInd = 1:length(iGrid) % for each inner point
    
    % hbar.iterate(1);
	 waitbar(ptInd/length(iGrid)  );

    ii = iGrid(ptInd); jj = jGrid(ptInd);
    
    if isnan(U(ii,jj)) == 1
        
        UNew(ptInd) = nan;  UX(ptInd) = nan;  UY(ptInd) = nan;
        
    else % Not nan
        
        LSMatrix = [ ones((2*Rad+1)^2,1), xGrid(:)*xstep, yGrid(:)*ystep ]; % Left hand side matrix
        LSb = U( sub2ind([M,N], xGrid(:)+ii, yGrid(:)+jj ) ); % Right hand side vector
         
        [row,col] = find(isnan(LSb)==0); % Remove nans of the surrounding points
        
        tempVector = LSMatrix(row,:)\LSb(row,:); % Least square fitting
        
        UNew(ptInd) = tempVector(1);
        UX(ptInd) = tempVector(2);
        UY(ptInd) = tempVector(3);
         
    end
    
end
close(hbar);

end
