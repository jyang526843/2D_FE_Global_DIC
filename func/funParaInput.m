function paraInput = funParaInput(paraName)

switch paraName
    
    case 'NewFFTSearch'
        fprintf('\n'); 
        fprintf('Since we are dealing with image sequences, for each new frame,   \n')
        fprintf('do we use last frame result as the initial guess or   \n')
        fprintf('Redo FFT initial guess for every new frame? \n    0: Use last frame (by default); \n    1: Redo initial guess.  \n')
        prompt = 'Input here: '; StillFFTSearch = input(prompt); paraInput = StillFFTSearch;
        fprintf('\n');
        
    case 'Subpb2FDOrFEM'
        fprintf('\n'); 
        fprintf('--- Method to solve ALDIC global step Subproblem 2 ---    \n')
        fprintf('    1: Finite difference(Recommended)   \n    2: Finite element method  \n')
        prompt = 'Input here: ';
        Subpb2FDOrFEM = input(prompt); paraInput = Subpb2FDOrFEM;
        
    case 'ClusterNo'
        fprintf('\n'); disp('--- Set up Parallel pool ---');
        fprintf('How many parallel pools to open? (Put in 1 if no parallel computing) \n');
        prompt = 'Input here: ';
        ClusterNo = input(prompt);
        % if ClusterNo > 1
        %     delete(gcp); myCluster = parcluster('local'); delete(myCluster.Jobs);
        %     parpool(ClusterNo,'SpmdEnabled',false);
        % end
        paraInput = ClusterNo;

    
    case 'SmoothDispOrNot' % Smooth displacements
        prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
        DoYouWantToSmoothOnceMore = input(prompt);
        paraInput = DoYouWantToSmoothOnceMore;
    
    case 'StrainMethodOp' % Choose strain computation method  
        fprintf('What method to use to compute strain? \n');
        fprintf('    0: Finite element (Recommended); \n');
        fprintf('    1: Finite difference; \n');
        fprintf('    2: Plane fitting; \n');
        prompt = 'Input here: ';
        MethodToComputeStrain = input(prompt);
        while (MethodToComputeStrain ~= 0) && (MethodToComputeStrain ~= 1) && (MethodToComputeStrain ~= 2) && (MethodToComputeStrain ~= 3)
            disp('****** Wrong input! ******')
            fprintf('What method to use to compute strain? \n');
            fprintf('    0: Finite element (Recommended); \n');
            fprintf('    1: Finite difference; \n');
            fprintf('    2: Plane fitting; \n');
            prompt = 'Input here: ';
            MethodToComputeStrain = input(prompt);
        end
        paraInput = [MethodToComputeStrain];
          
    case 'StrainType' % Choose strain computation method again
        fprintf('Infinitesimal stran or finite strain? \n');
        fprintf('    0: Infinitesimal stran; \n');
        fprintf('    1: Eulerian strain; \n');
        fprintf('    2: Green-Lagrangian strain; \n');
        fprintf('    3: Others: code by yourself; \n');
        prompt = 'Input here: ';
        StrainType = input(prompt);
        while (StrainType ~= 0) && (StrainType ~= 1) && (StrainType ~= 2) && (StrainType ~= 3)
            disp('****** Wrong input! ******')
            fprintf('Infinitesimal stran or finite strain? \n');
            fprintf('    0: Infinitesimal stran; \n');
            fprintf('    1: Eluerian strain; \n');
            fprintf('    2: Green-Lagrangian strain; \n');
            fprintf('    3: Others: code by yourself; \n');
            prompt = 'Input here: ';
            StrainType = input(prompt);
        end
        paraInput = [StrainType];
            
        
    case 'SaveFigFormat'
        fprintf('Save figures into different format: \n');
        fprintf('    1: jpeg (Choose transparency 0~1) \n');
        fprintf('    2: pdf (Choose transparency = 1) \n'); 
        fprintf('    3: pdf (Not to overlay original images) \n'); 
        fprintf('    4: Others: Edit codes in ./plotFiles/SaveFigFiles.m \n'); 
        prompt = 'Input here: '; MethodToSaveFig = input(prompt);
        paraInput = MethodToSaveFig;
        
    case 'OrigDICImgTransparency'
        fprintf('Define transparency for overlaying original images: \n')
        fprintf('Input a real number between 0 (Only original images) \n')
        fprintf('                        and 1 (Only deformation results).\n')
        prompt = 'Input here (e.g. 0.5): '; OrigDICImgTransparency = input(prompt);
        paraInput = OrigDICImgTransparency;
        
        
    otherwise
end
