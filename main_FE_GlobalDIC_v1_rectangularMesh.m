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
addpath('./func','./src','./plotFiles/','./plotFiles/export_fig-d966721/');  
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
fNormalized = ImgNormalized{1}; % Load the first referece image
% ====== Initialize variable storage ======
ResultDisp = cell(length(ImgNormalized)-1,1); 
ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrainWorld = cell(length(ImgNormalized)-1,1);
ResultStressWorld = cell(length(ImgNormalized)-1,1);
ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/DICpara.ImgSeqIncUnit),1); % For incremental DIC mode
ResultFEMeshEachFrame = cell(length(ImgNormalized)-1,1); % Needs future improvment: to store FE-mesh for each frame
ResultAlpha = cell(length(ImgNormalized)-1,1);
ResultNormOfW = cell(length(ImgNormalized)-1,1);
ResultTimeICGN = cell(length(ImgNormalized)-1,1);
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    gNormalized = ImgNormalized{ImgSeqNum}; % Load current deformed image frame 
    
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of the unknown displacements.
    % The key idea is to either to use a new FFT-based cross correlation peak fitting,
    % or use the results from the last frame as the new initial guess for the next frame;
    % Particularly in the incremental mode DIC, the reference image can also be updated, e.g.,
    % " fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)}; "
    %
    % DICpara.NewFFTSearch = 0; % If you want to apply the FFT-based cross correlation to 
    % compute the initial guess for each frame, please make sure that "DICpara.NewFFTSearch = 0". 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1 % Apply FFT-based cross correlation to compute the initial guess 
        % ====== Integer Search ======
        [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara);
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
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % To update ref image in incremental mode
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1,  fNormalizedNewIndex = fNormalizedNewIndex-1; end
        fNormalized = ImgNormalized{fNormalizedNewIndex}; % Update reference
        [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
        U0 = zeros(2*size(DICmesh.coordinatesFEM,1),1); % PlotuvInit;
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % Use the solved results from the last frame as the new initial guess
        if ImgSeqNum < 7 % Import previous U for ImgSeqNum [2,6] 
            U0 = ResultDisp{ImgSeqNum-2}.U;
             
        else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
            nTime = 5; np = length(ResultDisp{ImgSeqNum-2}.U)/2; % "nTime" value 5 is an empirical value, can be changed.
            T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np); 
            for tempi = 1:nTime
                T_data_u(tempi,:) = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(1:2:np*2)';
                T_data_v(tempi,:) = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(2:2:np*2)';
            end
            nB = 3; t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]'; t_pre = [ImgSeqNum-1]';
            [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
            [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
            tempu = u_pred(1,:); tempv = v_pred(1,:);
            U0 = [tempu(:),tempv(:)]'; U0 = U0(:);
         
            % %%%%% After running the new ImgSeqNum, you can uncomment these 
            % %%%%% lines to compare how the initial guess has been improved.  
            % Plotdisp_show(U0-ResultDisp{ImgSeqNum-1}.U,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
            % Plotdisp_show(ResultDisp{ImgSeqNum-2}.U-ResultDisp{ImgSeqNum-1}.U,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
        end
    end
    
    % ====== Compute image gradients ======
    Df = funImgGradient(fNormalized,gNormalized); % % using finite difference;
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM);
    fprintf('------------ Section 3 Done ------------ \n \n')

    
    %% Section 4
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finite element based global DIC iterations
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DICpara.tol = 1e-3;     % TODO: Iteration stopping threshold 
    DICpara.maxIter = 100;  % TODO: Maximum IC-GN iterations in IC-GN iterations
    DICpara.alpha = 10;     % TODO：Regularization coefficient, alpha, the default value is 10
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
        DICpara.GaussPtOrder = 2; [F] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,U,DICpara.GaussPtOrder);
        
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
        DICpara.GaussPtOrder = 2; [F] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,U,DICpara.GaussPtOrder);    
    end
    
    % ------- Smooth strain field --------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; 
    F = funSmoothStrain(F,DICmesh,DICpara);
     
    % ------- Save data -------
    ResultDisp{ImgSeqNum-1}.U = full(U);
    ResultDefGrad{ImgSeqNum-1}.F = full(F); % tempFoamAL;
    ResultAlpha{ImgSeqNum-1}.alpha = alpha_best;
    ResultNormOfW{ImgSeqNum-1}.normOfW = full(normOfW);
    ResultTimeICGN{ImgSeqNum-1}.timeICGN = full(timeICGN);
    
    fprintf('------------ Section 4 Done ------------ \n \n')
    
    
end    


% ------ Plot ------
% Transform into the physical world coordinates
UWorld = U; UWorld(2:2:end) = -U(2:2:end);
FWorld = F; FWorld(2:4:end) = -F(2:4:end); FWorld(3:4:end) = -F(3:4:end); 
close all; Plotuv(UWorld,DICmesh.x0,DICmesh.y0World); 
Plotdisp_show(UWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
Plotstrain_show(FWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);


%% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize),'_alpha',num2str(DICpara.alpha),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','normOfW','timeICGN');

   
%% Section 5
fprintf('------------ Section 5 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain fields and plot disp and strain results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Convert units from pixels to the physical world ------
DICpara.um2px = funParaInput('ConvertUnit');
% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = funParaInput('StrainMethodOp');
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType');
% ------ Choose image to plot (first only, second and next images) ------
if length(ImgNormalized)==2, DICpara.Image2PlotResults = funParaInput('Image2PlotResults');
else DICpara.Image2PlotResults = 1; % Plot over current, deformed image by default
end
% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
% ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1;
if DICpara.MethodToSaveFig == 1
    DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');
end

% ------ Start main part ------
for ImgSeqNum = 2 : length(ImgNormalized)
    
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
    
    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    xList = min(coordinatesFEM(:,1)):DICpara.winstepsize:max(coordinatesFEM(:,1)); M = length(xList);
    yList = min(coordinatesFEM(:,2)):DICpara.winstepsize:max(coordinatesFEM(:,2)); N = length(yList);
    [x0,y0] = ndgrid(xList,yList);
    x0 = x0-reshape(UFEMesh(1:2:end),size(x0,1),size(x0,2));
    y0 = y0-reshape(UFEMesh(2:2:end),size(y0,1),size(y0,2));
    x0World = DICpara.um2px*x0;
    y0World = DICpara.um2px*y0; % Ignore this: (size(ImgNormalized{1},2)+1-y0);
    coordinatesFEMWorld = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------ Plotting and Compute Strain-------
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; FLocal = FSubpb2.FSubpb2;
    else
        ULocal = USubpb2; FLocal = FSubpb2;
    end
    UWorld = DICpara.um2px*ULocal; UWorld(2:2:end) = -UWorld(2:2:end);  % close all; Plotuv(UWorld,x0,y0World);
    
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt);
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDisp(ULocal,DICmesh,DICpara);
            %%DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
     
    % ----- Compute strain field ------
    ComputeStrain; % Compute strain: run ComputeStrain.m
    % %%%%% Add filter and plot strain field %%%%%
    % %%%%% Plotstrain_Fij; %%%%%
    
    % ------ Plot disp and strain ------
    close all; % Plotuv(ULocal,x0,y0); 
    
    if DICpara.OrigDICImgTransparency == 1
        Plotdisp_show(UWorld,coordinatesFEMWorld,elementsFEM,DICpara);
        [strainxCoord,strainyCoord,dispu,dispv,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy,strain_principal_max, ...
            strain_principal_min,strain_maxshear,strain_vonMises]  =  Plotstrain0( ...
            UWorld,FStrainWorld,Rad,x0World,y0World,size(ImgNormalized{1}),DICpara);
    
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image
            Plotdisp(UWorld,x0World,y0World,size(ImgNormalized{1}),file_name{1,1},DICpara);
            [strainxCoord,strainyCoord,dispu,dispv,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy,strain_principal_max, ...
                strain_principal_min,strain_maxshear,strain_vonMises]  =  Plotstrain( ...
                UWorld,FStrainWorld,Rad,x0World,y0World,size(ImgNormalized{1}),file_name{1,1},DICpara);
        
        else % Plot over second or next deformed images
            Plotdisp(UWorld,x0World,y0World,size(ImgNormalized{1}),file_name{1,ImgSeqNum},DICpara);
            [strainxCoord,strainyCoord,dispu,dispv,dudx,dvdx,dudy,dvdy,strain_exx,strain_exy,strain_eyy,strain_principal_max, ...
                strain_principal_min,strain_maxshear,strain_vonMises]  =  Plotstrain( ...
                UWorld,FStrainWorld,Rad,x0World,y0World,size(ImgNormalized{1}),file_name{1,ImgSeqNum},DICpara);
        end
    end
    
    % ----- Save strain results ------
    ResultStrainWorld{ImgSeqNum-1} = struct('strainxCoord',strainxCoord,'strainyCoord',strainyCoord, ...
        'dispu',dispu,'dispv',dispv,'dudx',dudx,'dvdx',dvdx,'dudy',dudy,'dvdy',dvdy, ...
        'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy, ...
        'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
        'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
    
    % ------ Save figures for tracked displacement and strain fields ------
    SaveFigFilesDispAndStrain;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 5 Done ------------ \n \n')


%% Save data again including stress solve method
results_name = ['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize(1)),'_alpha',num2str(DICpara.alpha),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultStrainWorld', ...
     'ResultFEMesh','ResultFEMeshEachFrame',...
     'ResultAlpha','ResultNormOfW','ResultTimeICGN');
 

%% Section 6: Compute stress
fprintf('------------ Section 6 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute stress fields and plot stress fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Choose material model ------
DICpara.MaterialModel = funParaInput('MaterialModel');
% ------ Define parameters in material models ------
if (DICpara.MaterialModel == 1) || (DICpara.MaterialModel == 2) % Linear elasticity
    fprintf('Define Linear elasticity parameters \n')
    fprintf("Young's modulus (unit: Pa): \n"); prompt = 'Input here (e.g., 69e9): '; 
    DICpara.MaterialModelPara.YoungsModulus = input(prompt); 
    fprintf("Poisson's ratio: \n"); prompt = 'Input here (e.g., 0.3): '; 
    DICpara.MaterialModelPara.PoissonsRatio = input(prompt);
    fprintf('------------------------------------- \n');
end

% ------ Start main part ------
for ImgSeqNum = 2 : length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); close all;
    
    % ------ Plot stress ------
    if DICpara.OrigDICImgTransparency == 1
        [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
            stress_principal_min_xyplane, stress_maxshear_xyplane, ...
            stress_maxshear_xyz3d, stress_vonMises]  =  Plotstress0( ...
            DICpara,ResultStrainWorld{ImgSeqNum-1},size(ImgNormalized{1}));
        
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = Plotstress( ...
                DICpara,ResultStrainWorld{ImgSeqNum-1},size(ImgNormalized{1}),file_name{1,1});
            
        else % Plot over second or next deformed images
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = Plotstress( ...
                DICpara,ResultStrainWorld{ImgSeqNum-1},size(ImgNormalized{1}),file_name{1,ImgSeqNum});
            
        end
    end
    
    
    % ------ Save figures for computed stress fields ------
    SaveFigFilesStress;
    
    % ----- Save strain results ------
    ResultStressWorld{ImgSeqNum-1} = struct('stressxCoord',ResultStrainWorld{ImgSeqNum-1}.strainxCoord,'stressyCoord',ResultStrainWorld{ImgSeqNum-1}.strainyCoord, ...
        'stress_sxx',stress_sxx,'stress_sxy',stress_sxy,'stress_syy',stress_syy, ...
        'stress_principal_max_xyplane',stress_principal_max_xyplane, 'stress_principal_min_xyplane',stress_principal_min_xyplane, ...
        'stress_maxshear_xyplane',stress_maxshear_xyplane,'stress_maxshear_xyz3d',stress_maxshear_xyz3d, ...
        'stress_vonMises',stress_vonMises);
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 6 Done ------------ \n \n')



%% Save data again including stress solve method
results_name = ['results_FE_globalDIC_',imgname,'_st',num2str(DICpara.winstepsize(1)),'_alpha',num2str(DICpara.alpha),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultStrainWorld','ResultStressWorld', ...
    'ResultFEMesh','ResultFEMeshEachFrame',...
     'ResultAlpha','ResultNormOfW','ResultTimeICGN');

 

