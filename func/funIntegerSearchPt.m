

%% Integer Search on a non-uniform mesh

function [uSeedPt,vSeedPt,PhiSeedPt,tempSizeOfSearchRegion] = funIntegerSearchPt(f,g,winsize,varargin)
  
    seedPtCoords = varargin{1}; uSeedPt = zeros(size(seedPtCoords,1),1); vSeedPt = uSeedPt; PhiSeedPt = uSeedPt;

    fprintf('--- The size of initial guess search zone (pixels)? ---  \n')
    fprintf('User should start to try a small integer value, and gradually increase the value of \n');
    fprintf('the search zone size until it is larger than the magnitudes of |disp u| and |disp v|. \n');
    fprintf('User could also input [size_x, size_y] to search in a rectangular zone. \n');
    prompt = 'Input here: ';
    tempSizeOfSearchRegion = input(prompt);
    % if length(tempSizeOfSearchRegion) == 1, tempSizeOfSearchRegion = tempSizeOfSearchRegion*[1,1]; end


    for tempi = 1:length(seedPtCoords)
        if seedPtCoords(tempi,1)>0 
        if ceil(seedPtCoords(tempi,1)-winsize/2)-tempSizeOfSearchRegion < 1 || ...
                floor(seedPtCoords(tempi,1)+winsize/2)+tempSizeOfSearchRegion > size(g,1) || ...
                ceil(seedPtCoords(tempi,2)-winsize/2)-tempSizeOfSearchRegion < 1 || ...
                floor(seedPtCoords(tempi,2)+winsize/2)+tempSizeOfSearchRegion > size(g,2)
            
            uSeedPt(tempi) = nan;
            vSeedPt(tempi) = nan;
            PhiSeedPt(tempi) = nan;
            continue;
        else

            C = f(ceil(seedPtCoords(tempi,1)-winsize/2):floor(seedPtCoords(tempi,1)+winsize/2), ...
                ceil(seedPtCoords(tempi,2)-winsize/2):floor(seedPtCoords(tempi,2)+winsize/2));
            D = g(ceil(seedPtCoords(tempi,1)-winsize/2)-tempSizeOfSearchRegion:floor(seedPtCoords(tempi,1)+winsize/2)+tempSizeOfSearchRegion, ...
                ceil(seedPtCoords(tempi,2)-winsize/2)-tempSizeOfSearchRegion:floor(seedPtCoords(tempi,2)+winsize/2)+tempSizeOfSearchRegion);

            XCORRF2OfCD0 = normxcorr2(C,D);

            [v1temp, u1temp, max_f] = findpeak(XCORRF2OfCD0(winsize:end-winsize+1,winsize:end-winsize+1),1);

            zero_disp = ceil(size(XCORRF2OfCD0(winsize:end-winsize+1,winsize:end-winsize+1))/2);

            uSeedPt(tempi) = u1temp-zero_disp(1);
            vSeedPt(tempi) = v1temp-zero_disp(2);
            PhiSeedPt(tempi) = max_f;

        end
        else
            uSeedPt(tempi) = nan;
            vSeedPt(tempi) = nan;
            PhiSeedPt(tempi) = nan;
        end

    end


end


%% ==============================================

%%
function qfactors = compute_qFactor(cc,qnum)

%get peak locations and cc_min maps (i.e. cc - cc(min))
[peak,cc_min] = cellfun(@(x) cc_max_find(double(x)),cc.A,'UniformOutput',0);

%compute two primary quality metrics, as given in "Xue Z, Particle Image
% Velocimetry Correlation Signal-to-noise Metrics, Particle Image
% Pattern Mutual Information and Measurement uncertainty Quantification.
% MS Thesis, Virginia Tech, 2014.

%peak to corr. energy ratio
pce = cellfun(@(x,y) (abs(y)^2)/(1/numel(x)*(sum(abs(x(:)).^2))),cc_min,peak,'UniformOutput',0);
%min value -> 1 (worst case)

%peak to entropy ratio
ppe = cellfun(@(x) q_entropy(double(x)),cc_min,'UniformOutput',0);%peak to cc (information) entropy
%min value -> 0 (worst case)

qfactors = cell2mat(...
    cellfun(@(x,y) [x(:);y(:)], pce,ppe,'UniformOutput',0))';

end

function [peak,cc_min] = cc_max_find(cc)
%find the peak and zero-adjusted cc map

cc_min = cc - min(cc(:));%zero-adjust
% cc_filt = imgaussfilt3(cc_min); %filter to remove noise from peak value

[peak,~] = max(cc_min(:)); %get the index of the peak


end

function [ppe] = q_entropy(cc_min)
%compute entropy q-factor for a given cc map

[cc_hist,~] = histcounts(cc_min,30); %get histogram values

entropy = 0;
p = cc_hist/sum(cc_hist); %compute probablities
for i = 1:length(p)%compute entropy
    if p(i) == 0
        entropy = entropy+p(i);
    else
        entropy = entropy+p(i)*log(1/p(i));
    end
end

ppe = 1/entropy; %peak to cc (information) entropy
%min value -> 0 (worst case)


end









