% ==============================================
% function funSmoothDispCrack
% ==============================================
function FLocal = funSmoothStrain(FLocal,DICmesh,DICpara)

coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;
winstepsize = DICpara.winstepsize;
StrainFilterSize = DICpara.StrainFilterSize;
StrainFilterStd = DICpara.StrainFilterStd;

FilterStd = StrainFilterStd; FilterSizeInput = StrainFilterSize; LevelNo = 1;
% switch nargin
%     case 5
%         FilterSizeInput = varargin{1};
%     case 6
%         FilterSizeInput = varargin{1}; FilterStd = varargin{2};
%     case 7
%         FilterSizeInput = varargin{1}; FilterStd = varargin{2}; LevelNo = varargin{3};
%     otherwise
%         disp('Wrong input in funSmoothStrain!');
% end

%%
% close all; Plotdisp_show(ULocal,elementsFEM,coordinatesFEM);

% prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
DoYouWantToSmoothOnceMore = 0; % DoYouWantToSmoothOnceMore = input(prompt);
if DoYouWantToSmoothOnceMore == 0  
    if isempty(FilterStd) == 1
        prompt = 'Choose filter standard deviation(0-default): ';
        FilterStd = input(prompt);
        if FilterStd == 0
            FilterStd = 0.5; 
        end
    else
        if FilterStd == 0
            FilterStd = 0.5;
        end
    end
    if isempty(FilterSizeInput) == 1
        prompt = 'Choose Gaussian filter size(0-default): ';
        FilterSizeInput = input(prompt);
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1; 
        end
    else
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1;
        end
    end
end

SmoothTimes = 1;
while (DoYouWantToSmoothOnceMore==0)
    Coordxnodes = [min(coordinatesFEM(:,1)):winstepsize/(2^(LevelNo-1)):max(coordinatesFEM(:,1))]'; 
    Coordynodes = [min(coordinatesFEM(:,2)):winstepsize/(2^(LevelNo-1)):max(coordinatesFEM(:,2))]';
    Iblur_Top11 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(1:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
    Iblur_Top11=Iblur_Top11';
    Iblur_Top22 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(4:4:end),Coordxnodes,Coordynodes,'regularizer','springs');  
    Iblur_Top22=Iblur_Top22';
    Iblur_Top21 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(2:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
    Iblur_Top21=Iblur_Top21';
    Iblur_Top12 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(3:4:end),Coordxnodes,Coordynodes,'regularizer','springs');
    Iblur_Top12=Iblur_Top12';
    % Iblur_Top = nan(size(ULocal,1),1); Iblur_Top(2*CoordCrackTop-1) = ULocal(2*CoordCrackTop-1); Iblur_Top(2*CoordCrackTop) = ULocal(2*CoordCrackTop); 
    % Iblur_Top10 = reshape(Iblur_Top(1:2:end),M,N); Iblur_Top20 = reshape(Iblur_Top(2:2:end),M,N);
    % -------------------------------------------------------
    % Iblur_Top1 = reshape(imgaussfilt(Iblur_Top10,FilterStd,'FilterSize',FilterSizeInput,'FilterDomain','spatial'), M*N, 1);
    % Iblur_Top2 = reshape(imgaussfilt(Iblur_Top20,FilterStd,'FilterSize',FilterSizeInput,'FilterDomain','spatial'), M*N, 1);
    % -------------------------------------------------------
    imageFilter=fspecial('gaussian',FilterSizeInput,FilterStd);
    % Iblur_Top1 = reshape(nanconv(Iblur_Top10,imageFilter,'edge','nanout'), M*N, 1);
    % Iblur_Top2 = reshape(nanconv(Iblur_Top20,imageFilter,'edge','nanout'), M*N, 1);
    % ULocal(2*CoordCrackTop-1) = Iblur_Top1(CoordCrackTop); ULocal(2*CoordCrackTop) = Iblur_Top2(CoordCrackTop);
    Iblur_Top1 = nanconv(Iblur_Top11,imageFilter,'edge','nanout');
    Iblur_Top4 = nanconv(Iblur_Top22,imageFilter,'edge','nanout');
    Iblur_Top2 = nanconv(Iblur_Top21,imageFilter,'edge','nanout');
    Iblur_Top3 = nanconv(Iblur_Top12,imageFilter,'edge','nanout');
    
    for tempi = 1:size(coordinatesFEM,1)
        [row1,col1] = find(Coordxnodes==coordinatesFEM((tempi),1));
        [row2,col2] = find(Coordynodes==coordinatesFEM((tempi),2));
        FLocal(4*(tempi)-3) = Iblur_Top1(row1,row2);
        FLocal(4*(tempi))   = Iblur_Top4(row1,row2);
        FLocal(4*(tempi)-2) = Iblur_Top2(row1,row2);
        FLocal(4*(tempi)-1) = Iblur_Top3(row1,row2);
    end
    
    % close all; Plotuv(ULocal,x,y); % Plotting u and v
    % close all; Plotdisp_show(ULocal,elementsFEM,coordinatesFEM);
     
    % prompt = 'Do you want to smooth displacement once more? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = SmoothTimes+1;
    if SmoothTimes > 2
        DoYouWantToSmoothOnceMore = 1;
    end
    
end