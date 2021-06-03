%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update reference image for incremental mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ULast,DICpara,DICmesh)

ImgSeqIncROIUpdateOrNot = DICpara.ImgSeqIncROIUpdateOrNot;
gridxyROIRange = DICpara.gridxyROIRange;
winstepsize = DICpara.winstepsize;
coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;
M = DICmesh.M; N = DICmesh.N; 

%% ------ Use fixed mesh -or- Manually update mesh ------ 
if ImgSeqIncROIUpdateOrNot == 1
    gridxyROIRange = funReadImageRefUpdate(file_name{1,ImgSeqNum-1});
    gridxROIRangeNew = [ceil(gridxyROIRange.gridx(1)) : winstepsize : floor(gridxyROIRange.gridx(2))];
    gridyROIRangeNew = [ceil(gridxyROIRange.gridy(1)) : winstepsize : floor(gridxyROIRange.gridy(2))];
    [x0Newtemp,y0Newtemp] = ndgrid(gridxROIRangeNew,gridyROIRangeNew); x0Newtemp=x0Newtemp'; y0Newtemp=y0Newtemp';
    
    [DICmesh] = MeshSetUp(x0Newtemp,y0Newtemp,DICpara);
    DICpara.gridxyROIRange = gridxyROIRange;
    
%% ------ Update ROI automatically ------
elseif ImgSeqIncROIUpdateOrNot == 2
    % Try to automatically update FEM mesh when updating reference img,
    % but sometimes doesn't work very well and error could accumulate.
    % ---------------------------------------------------------------
    %ULast = ResultDisp{ImgSeqNum-2}.U;
    coordinatesFEMOld = coordinatesFEM; elementsFEMOld = elementsFEM;
    coordinatesFEM(:,1) = (coordinatesFEM(:,1) + ULast(1:2:end));
    coordinatesFEM(:,2) = (coordinatesFEM(:,2) + ULast(2:2:end));
    
    coordinatesFEMNewBottomInd = 1:1:M; coordinatesFEMNewTopInd = M*N-(M-1):1:M*N;
    coordinatesFEMNewLeftInd = 1:M:M*N; coordinatesFEMNewRightInd = M:M:M*N;
    gridxyROIRange.gridxROIRange = [max(coordinatesFEM(coordinatesFEMNewLeftInd,1)), min(coordinatesFEM(coordinatesFEMNewRightInd,1))];
    gridxyROIRange.gridyROIRange = [max(coordinatesFEM(coordinatesFEMNewBottomInd,2)), min(coordinatesFEM(coordinatesFEMNewTopInd,2))];
    
    DICpara.gridxyROIRange = gridxyROIRange;
    DICmesh.coordinatesFEM = coordinatesFEM; 
    
end
 