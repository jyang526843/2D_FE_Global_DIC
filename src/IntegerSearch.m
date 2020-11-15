% ==============================================
% function AdjSearchDomain
% ==============================================
function [SizeOfFFTSearchRegion,x0,y0,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara)

gridxROIRange = DICpara.gridxyROIRange.gridx;
gridyROIRange = DICpara.gridxyROIRange.gridy;
winsize = DICpara.winsize;
winstepsize = DICpara.winstepsize;




% Input initial integer search zone size
InitialGuessSatisfied = 1; % gridxBackup = gridxROIRange; gridyBackup = gridyROIRange;

% ========= Start FFT to find initial guess with integer accuracy =======
while InitialGuessSatisfied == 1
    
    fprintf('--- Whole field initial guess (0) OR Several seeds (1) initial guess? ---  \n')
    prompt = 'Input here: ';
    tempNoOfInitPt = input(prompt);
    
    fprintf('--- What is your initial guess search zone size (pixels)? ---  \n')
    prompt = 'Input here: ';
    tempSizeOfSearchRegion = input(prompt);
    
    if length(tempSizeOfSearchRegion) == 1, tempSizeOfSearchRegion = tempSizeOfSearchRegion*[1,1]; end
    
    
    if tempNoOfInitPt ~= 1  
       
        tempNoOfInitPt = 0;
        % tempNoOfInitPt == 0, whole field for initial guess, 
        
        [x0,y0,u,v,cc] = funIntegerSearch(fNormalized,gNormalized,tempSizeOfSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,tempNoOfInitPt,winstepsize);

    else % tempNoOfInitPt ~= 0, several local seeds for initial guess
        
        % input local seeds coordinates
        figure; imshow( (imread(file_name{1}))); % surf(fNormalized,'EdgeColor','none','LineStyle','none'); view(2);
        [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
        
        % initial FFT search for local seeds
        [x0,y0,u,v,cc] = funIntegerSearch(fNormalized,gNormalized,tempSizeOfSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,tempNoOfInitPt,[row,col]);

        % plot-view local seeds initial results
         
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Have a look at integer search
    % --------------------------------------
    close all;
    figure; surf(u); colorbar;
    title('Displacement u','fontweight','normal')
    set(gca,'fontSize',18);
    title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
    axis tight; % set(gca,'XTick',[] );
    xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w');
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    box on; colormap jet;


    figure; surf(v); colorbar;
    title('Displacement v','fontweight','normal')
    set(gca,'fontSize',18);
    title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
    axis tight; % set(gca,'XTick',[] );
    xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w');
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    box on; colormap jet;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('--- Are you satisfied with initial guess with current search region? (0-yes; 1-no)? ---  \n')
    prompt = 'Input here: ';
    InitialGuessSatisfied = input(prompt);
    
end

% ======== Find some bad inital guess points ========
cc.ccThreshold = 1.25; % bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)
qDICOrNot = 1; Thr0 = 100; [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0);

%[u,v] = funRemoveOutliers(u,v);

% ======== Output final search region radius ========
SizeOfFFTSearchRegion = tempSizeOfSearchRegion;


