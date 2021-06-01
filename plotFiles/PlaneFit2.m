function [UY,UX,UNew] = PlaneFit2(U,winsizeOfx,winsizeOfy, Rad)

[M,N] = size(U); UY = U; UX = U; UNew = U;

% [tempjj,tempii] = meshgrid(Rad+1:N-Rad, Rad+1:M-Rad); tempii = tempii(:); tempjj = tempjj(:);
 
%hbar = parfor_progressbar(length(tempii),'Please wait for PlaneFit2!');
h=waitbar(0,'wait for PlaneFit2!');

 
 for ii = (Rad+1):M-Rad
  for jj = (Rad+1):N-Rad
        
        LSMatrix = ones((2*Rad+1)^2, 3);
        for tempj = -Rad:Rad
            for tempi = -Rad:Rad
                LSMatrix((2*Rad+1)*(tempj+Rad)+tempi+Rad+1,:) = ...
                    [1 tempi*winsizeOfx tempj*winsizeOfy];
            end
        end

        LSb = zeros((2*Rad+1)^2, 1);
        for tempj = -Rad:Rad
            for tempi = -Rad:Rad
                LSb((2*Rad+1)*(tempj+Rad)+tempi+Rad+1,:) = ...
                    U(ii+tempi, jj+tempj);
            end
        end
         
        tempVector = (LSMatrix'*LSMatrix)\(LSMatrix'*LSb);
        UNew(ii,jj) = tempVector(1);
        UX(ii,jj) = tempVector(2);
        UY(ii,jj) = tempVector(3);
         
% end
% hbar.iterate(1);
tempij = (ii -Rad-1)*(N-2*Rad)+jj-Rad ;
waitbar(tempij/((M-2*Rad)*(N-2*Rad))  );
  end
end




