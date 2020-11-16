% ---------------------------------------------
% Finite-element-based global DIC (FE-Global-DIC)
% Author: Jin Yang, Postdoc @UW-Madison;  PhD @Caltech 19'; 
% Contact: aldicdvc@gmail.com; jyang526@wisc.edu
% 2015.04,06,07; 2016.03,04; 2020.11
% ---------------------------------------------

%% Section 1: Clear MATLAB environment & mex set up Spline interpolation  
close all; clear; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
% % cd("./Splines_interp/lib_matlab"); CompileLib; cd("../../");  % % mex bi-cubic spline interpolations
% % addpath("./Splines_interp/lib_matlab"); % dbstop if error % % Old version codes.
mex -O ba_interp2.cpp; 
addpath("./func"); addpath("./src"); addpath("./plotFiles/");  
% addpath("./YOUR IMAGE FOLDER"); 
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ======
[file_name,Img,DICpara] = ReadImage; close all; 
% %%%%%% Uncomment the line below to change the DIC computing ROI manually %%%%%%
%gridxROIRange = [gridxROIRange1,gridxROIRange2]; gridyROIRange = [Val1, Val2];
%gridxROIRange = [224,918]; gridyROIRange = [787,1162];
% ====== Normalize images ======
[ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(Img,DICpara.gridxyROIRange);  
% ====== Initialize variable storage ======
ResultDisp = cell(length(ImgNormalized)-1,1); 
ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrain = cell(length(ImgNormalized)-1,1);
ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/DICpara.ImgSeqIncUnit),1); % For incremental DIC mode
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fNormalized = ImgNormalized{1}; 
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); fprintf('\n');
    
    %% Section 3:
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update initial guess for ALDIC
    % The key idea is to either to use FFT peak fitting, or use last frame
    % results for the next new frame;
    % Particularly in incremental mode, the reference image can also be updated.
    % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
    gNormalized = ImgNormalized{ImgSeqNum}; NewFFTSearchCheck = 0; DICpara.NewFFTSearch = 0;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1
        % ====== Integer Search ======
        % Old version: search a subset of img1 in a larger region of img2 
        [DICpara.SizeOfFFTSearchRegion,x0temp,y0temp,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara);
        % New version: adaptive search initial guess
        % [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);
        % ====== FEM mesh set up ======
        [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0); %PlotuvInit; [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);
        % ====== Deal with incremental mode ======
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
       
    elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % TO update ref image in incremental mode
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1,  fNormalizedNewIndex = fNormalizedNewIndex-1; end
        fNormalized = ImgNormalized{fNormalizedNewIndex}; % Update reference
        [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
        U0 = zeros(2*size(DICmesh.coordinatesFEM,1),1); % PlotuvInit;
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    %Plotdisp_show(U0,[coordinatesFEM(:,1),size(fNormalized,2)+1-coordinatesFEM(:,2)],elementsFEM); % Plot initial values
    % ====== Spline interpolation images ======
    %[imgfNormalizedbc,imggNormalizedbc,imgSize,DfAxis] = funImgGradient(fNormalized,gNormalized);
    Df = funImgGradient(fNormalized,gNormalized); % % using finite difference;
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    fprintf('------------ Section 3 Done ------------ \n \n')

    
    %% Section 4
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finite element based global DIC iterations
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DICpara.tol = 1e-4; % iteration stopping threshold 
    DICpara.maxIter = 100; % Maximum IC-GN iterations in IC-GN iterations
    DICpara.alpha = 1e2; % Set regularization coefficient, alpha, as 10
    alphaList = DICpara.alpha;
    
    % ====== Tune regularization coefficient ======
    % If you don't know the best alpha (coefficient), please run the following
    % codes to tune the best value of the coefficient of the regularizer |grad u|^2):
    % 
    % %%%%% Uncomment the following line to tune the best value of alpha %%%%%%
    % alphaList = [1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3]*mean(DICpara.winstepsize);
    
    Err1List = zeros(length(alphaList),1); Err2List=Err1List; UList = cell(length(alphaList),1); FList=UList;
     
    % ------------------------------------------------ 
    for alphaInd = 1:length(alphaList)
        
        tic; alpha = alphaList(alphaInd); 
        % Solve displacement U with each alpha
        [U, normOfW, timeICGN] = funGlobalICGN(DICmesh,Df,fNormalized,gNormalized,U0,alpha,DICpara.tol,DICpara.maxIter);
        % Compute F deformation gradient with solved U
        DICpara.GaussPtOrder = 2; [F] = funGlobal_NodalStrainAvg(DICmesh,U,DICpara.GaussPtOrder);
        
        Err1List(alphaInd) = norm(U-U0,2);
        Err2List(alphaInd) = norm(F,2); 
         
    end
    
    % ====== Tune the coefficient of |grad u| regularizer ======
    ErrSumList = Err1List + 1*mean(DICpara.winstepsize)*Err2List;  % 10 is an empirical number
    [~,indexOfalpha] = min(ErrSumList); 
    
    try
        [fitobj] = fit(log10(alphaList(indexOfalpha-1:1:indexOfalpha+1))',ErrSumList(indexOfalpha-1:1:indexOfalpha+1),'poly2');
        p = coeffvalues(fitobj); alpha_best = 10^(-p(2)/2/p(1)); 
    catch
        alpha_best = alphaList(indexOfalpha); 
    end
    DICpara.alpha = alpha_best; % If you just have one item in the alphaList, this line doens't change anything.
    
    
    % ====== Re-run global DIC iterations with tuned alpha_best ======
    if abs(alpha_best - alpha) > abs(eps)
        [U, normOfW, timeICGN] = funGlobalICGN(DICmesh,Df,fNormalized,gNormalized,U0,alpha_best,DICpara.tol,DICpara.maxIter);
        DICpara.GaussPtOrder = 2; [F] = funGlobal_NodalStrainAvg(DICmesh,U,DICpara.GaussPtOrder);    
    end
    
    % ------- Smooth strain field --------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; 
    F = funSmoothStrain(F,DICmesh,DICpara);
     
    % ------- Save data -------
    ResultDisp{ImgSeqNum-1}.U = full(U);
    ResultDefGrad{ImgSeqNum-1}.F = full(F); % tempFoamAL;
    
    fprintf('------------ Section 4 Done ------------ \n \n')
    
    
end    


% ------ Plot ------
UWorld = U; UWorld(2:2:end) = -U(2:2:end);
FWorld = F; % Transform into the physical world coordinates
close all; Plotuv(UWorld,DICmesh.x0,DICmesh.y0World); 
Plotdisp_show(UWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
Plotstrain_show(F,coordinatesFEMWorld,elementsFEM);


%% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize),'_alpha',num2str(DICpara.alpha),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','normOfW','timeICGN');

   
%% Section 5
fprintf('------------ Section 5 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = funParaInput('StrainMethodOp'); 
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType');
% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
DICpara.OrigDICImgTransparency = 1; 
if DICpara.MethodToSaveFig == 1  
    DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');         
end

% ------ Start main part ------
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    if DICpara.ImgSeqIncUnit > 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit);
    elseif DICpara.ImgSeqIncUnit == 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)-1;
    end
    FEMeshInd = FEMeshIndLast + 1;
    
    if FEMeshInd == 1
        USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
        coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
        elementsFEM = ResultFEMesh{1}.elementsFEM;
        if (ImgSeqNum-1 == 1) || (DICpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    else
        USubpb2 = ResultDisp{ImgSeqNum-1}.U;
        if mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0
            coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
            elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
            coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
            UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
            xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
            UFEMesh = 0*USubpb2;
            UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
            UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
        end
        USubpb2 = USubpb2 + UFEMesh;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xList = min(coordinatesFEM(:,1)):DICpara.winstepsize:max(coordinatesFEM(:,1)); M = length(xList);
    yList = min(coordinatesFEM(:,2)):DICpara.winstepsize:max(coordinatesFEM(:,2)); N = length(yList);
    %[x0,y0] = ndgrid(xList,yList); 
    [x0,y0] = ndgrid(xList,yList); 
    x0 = x0-reshape(UFEMesh(1:2:end),size(x0,1),size(x0,2)); 
    y0 = y0-reshape(UFEMesh(2:2:end),size(y0,1),size(y0,2)); 
    y0World = (size(ImgNormalized{1},2)+1-y0); 
    coordinatesFEMWorld = [coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];

    % ------ Plotting and Compute Strain-------
    M = size(x0,1); N = size(x0,2);
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; FLocal = FSubpb2.FSubpb2; 
    else
        ULocal = USubpb2; FLocal = FSubpb2;
    end
    % ULocal(1:2:end)= -ULocal(1:2:end); 
    % FLocal(1:4:end)= -1*FLocal(1:4:end); FLocal(3:4:end)= -1*FLocal(3:4:end); % because u is flipped sign
    UWorld = ULocal; UWorld(2:2:end) = -UWorld(2:2:end); FWorld = FLocal; close all; Plotuv(UWorld,x0,y0World);
    % tic; D = funDerivativeOp(M,N,winstepsize); toc

    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDisp(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
    % ----- Compute strain field ------
    ComputeStrain; % Compute strain components
    % % ------- Add additional filter and plot strain field -------
    % Plotstrain_Fij; % Optional

    % ------ Plot disp and strain ------
    if DICpara.MethodToSaveFig > 1
        Plotdisp_show(UWorld,coordinatesFEMWorld,elementsFEM);
        Plotstrain0(FStraintemp,x0(1+Rad:M-Rad,1+Rad:N-Rad),y0(1+Rad:M-Rad,1+Rad:N-Rad),size(ImgNormalized{1}),...
        file_name{1,ImgSeqNum},DICpara.OrigDICImgTransparency);  
    else
        Plotdisp(UWorld,x0,y0,size(ImgNormalized{1}),file_name{1,ImgSeqNum},DICpara.OrigDICImgTransparency);
        Plotstrain(UWorld,Rad,FStraintemp,x0(1+Rad:M-Rad,1+Rad:N-Rad),y0(1+Rad:M-Rad,1+Rad:N-Rad),size(ImgNormalized{1}),...
        file_name{1,ImgSeqNum},DICpara.OrigDICImgTransparency); 
    end

    ResultStrain{ImgSeqNum-1}.Strain = FStraintemp;
    
    % ------ Save figures ------
    SaveFigFiles;
    
end
% ------ End of image sequence! ------
fprintf('------------ Section 5 Done ------------ \n \n')

  

%% Save data again including strain solve method
results_name = ['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize),'_alpha',num2str(DICpara.alpha),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultStrain','ResultFEMesh',...
     'normOfW','timeICGN');

 

