function paraInput = funParaInput(paraName)
% FUNCTION paraInput = funParaInput(paraName)
% To define DIC parameters from user's input
% ----------------------------------------------
%
%   INPUT:  paraName        Parameter name
%
%   OUTPUT: paraInput       Input DIC parameters
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
% ==============================================

%%
switch paraName
    
    case 'InitFFTSearchMethod' % FFT-solver to compute the initial guess
        fprintf('\n');
        fprintf('--- Method to compute an initial guess of displacements --- \n');
        fprintf('    0: Multigrid search based on an image pyramid  \n');
        fprintf('    1: Whole field search for all the subsets  \n');
        fprintf('    2: Search near manually clicked seeds and then interpolate for the full-field  \n');
        prompt = 'Input here: '; InitFFTSearchMethod = input(prompt); paraInput = InitFFTSearchMethod;
        fprintf('\n');
    
        
    case 'NewFFTSearch' % For an image sequence: re-run FFT for an initial guess or not?
        fprintf('\n'); 
        fprintf('Since we are dealing with an image sequence, for each new frame,   \n')
        fprintf('do we use last frame result as the initial guess or   \n')
        fprintf('Redo FFT initial guess for every new frame? \n')
        fprintf('    0: Use last frame (by default); \n');
        fprintf('    1: Redo initial guess.  \n');
        prompt = 'Input here: '; StillFFTSearch = input(prompt); paraInput = StillFFTSearch;
        fprintf('\n');
        
        
    case 'Subpb2FDOrFEM' % Solve ALDIC Subproblem 2 via a 'dinite difference' or 'fnite element' solver?
        fprintf('\n'); 
        fprintf('--- Method to solve ALDIC global step Subproblem 2 ---    \n')
        fprintf('    1: Finite difference (Recommend)   \n')
        fprintf('    2: Finite element method  \n')
        prompt = 'Input here: ';
        Subpb2FDOrFEM = input(prompt); paraInput = Subpb2FDOrFEM;
        
        
    case 'ClusterNo' % Parallel threads #
        fprintf('\n'); disp('--- Set up Parallel pool ---');
        fprintf('How many parallel pools to open? (Put in 1 if no parallel computing) \n');
        prompt = 'Input here: ';
        ClusterNo = input(prompt);
        % if ClusterNo > 1 % Backup codes
        %     delete(gcp); myCluster = parcluster('local'); delete(myCluster.Jobs);
        %     parpool(ClusterNo,'SpmdEnabled',false);
        % end
        paraInput = ClusterNo;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Section 8 %%%%% 
    case 'ConvertUnit' % Convert units from pixels to the physical world
        fprintf("Convert units from pixels to the physical world. \n");
        display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf("Results in Section 5 (ResultStrain), Section 6 \n(ResultStress), and all the plots will be converted \nto the physical world unit instead of the pixel unit.  \n");
        display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf("If you want to keep the pixel unit, enter '1'. \n")
        prompt = 'Input here (e.g., mm/px, um/px): ';
        um2px = input(prompt);
        paraInput = um2px;
        fprintf('------------------------------------- \n');
    
    case 'SmoothDispOrNot' % Smooth displacements or not
        fprintf('Do you want to smooth displacement? (0-yes; 1-no) \n');
        prompt = 'Input here: ';
        DoYouWantToSmoothOnceMore = input(prompt);
        paraInput = DoYouWantToSmoothOnceMore;
        fprintf('------------------------------------- \n');
        
        
    case 'RegularizationSmoothness' % Curvature regularization of displacements ?
        fprintf('Regularization smoothness \n');
        prompt = 'Input here (e.g., 0 or 1e-5~1e-3): ';
        DoYouWantToSmoothOnceMore = input(prompt);
        paraInput = DoYouWantToSmoothOnceMore;
        fprintf('------------------------------------- \n');
        
        
    case 'StrainMethodOp' % Choose strain computation method  
        fprintf('What method to use to compute strain? \n');
        fprintf('    0: Finite element (Recommend); \n');
        fprintf('    1: Finite difference ; \n');
        fprintf('    2: Plane fitting; \n');
        prompt = 'Input here: ';
        MethodToComputeStrain = input(prompt);
        while (MethodToComputeStrain ~= 0) && (MethodToComputeStrain ~= 1) && (MethodToComputeStrain ~= 2) && (MethodToComputeStrain ~= 3)
            disp('****** Wrong input! ******')
            fprintf('What method to use to compute strain? \n');
            fprintf('    0: Finite element (Recommend); \n');
            fprintf('    1: Finite difference ; \n');
            fprintf('    2: Plane fitting; \n');
            prompt = 'Input here: ';
            MethodToComputeStrain = input(prompt);
        end
        paraInput = [MethodToComputeStrain];
        fprintf('------------------------------------- \n');
        
        
    case 'StrainType' % Define strain type
        fprintf('Infinitesimal strain or finite strain? \n');
        fprintf('    0: Infinitesimal strain; \n');
        fprintf('    1: Eulerian strain; \n');
        fprintf('    2: Green-Lagrangian strain; \n');
        fprintf('    3: Others: code by yourself; \n');
        prompt = 'Input here: ';
        StrainType = input(prompt);
        while (StrainType ~= 0) && (StrainType ~= 1) && (StrainType ~= 2) && (StrainType ~= 3)
            disp('****** Wrong input! ******')
            fprintf('Infinitesimal strain or finite strain? \n');
            fprintf('    0: Infinitesimal strain; \n');
            fprintf('    1: Eluerian strain; \n');
            fprintf('    2: Green-Lagrangian strain; \n');
            fprintf('    3: Others: code by yourself; \n');
            prompt = 'Input here: ';
            StrainType = input(prompt);
        end
        paraInput = StrainType;
        fprintf('------------------------------------- \n');
        
        
    case 'Image2PlotResults' % Choose image to plot results (first only, second and next images) ------
        fprintf('Over which image(s) you want to plot the results? \n');
        fprintf('    0: First image only (if you are using only two images); \n');
        fprintf('    1: Second and next images; \n');
        prompt = 'Input here: '; Image2PlotResults = input(prompt);
		paraInput = Image2PlotResults;
        fprintf('------------------------------------- \n');
        
        
    case 'SaveFigFormat' % Finally saved figure format
        fprintf('Save figures into the format: \n');
        fprintf('    1: jpg(Choose transparency 0~1) \n');
        fprintf('    2: pdf(Choose transparency = 1) \n'); 
        fprintf('    3: Others: Edit codes in ./plotFiles/SaveFigFiles.m \n'); 
        prompt = 'Input here: '; MethodToSaveFig = input(prompt);
        paraInput = MethodToSaveFig;
        fprintf('------------------------------------- \n');
        
        
    case 'OrigDICImgTransparency' % DIC results plotting transparency
        fprintf('Define transparency for overlaying original images: \n')
        fprintf('Input a real number between 0(Only original images) \n')
        fprintf('and 0.99 (Non-transparent deformation results).\n')
        prompt = 'Input here(e.g., 0.5): '; OrigDICImgTransparency = input(prompt);
        paraInput = OrigDICImgTransparency;
        fprintf('------------------------------------- \n');
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Section 9 %%%%%
    case 'MaterialModel' % Material model: plane stress or plane strain?
        fprintf('Material model to compute Cauchy stress fields: \n');
        fprintf('    1: Linear elasticity -- Plane stress \n');
        fprintf('    2: Linear elasticity -- Plane strain \n');
        fprintf('    3: Others: User needs to edit codes in ./plotFiles/Plotstress.m   \n');
        prompt = 'Input here: '; MaterialModel = input(prompt);
        paraInput = MaterialModel;
        fprintf('------------------------------------- \n');
        
        
    otherwise
end



