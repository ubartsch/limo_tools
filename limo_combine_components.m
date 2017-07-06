function combined_ica = limo_combine_components(varargin)

% compute a weighted average of components
%
% FORMAT combined_ica = limo_combine_components(data,weights,invweights,whichic)
%
% INPUT data a 3D matrix representing set of components 
%            dim compoments*time frames*epochs
%       weights are the columns EEG.weights
%       invweights are the columns EEG.icawinv
%       whichic are the com ponents to 'combine'
%       method 'sum' 'max' 'maxvar' (default)
%
% OUTPUT combined_ica is the weighted average of comnponents
%
% Cyril Pernet & Ramon Martinez-Cancino 23-10-2014 updates for components (ICA)
%% ----------------------------------
%  Copyright (C) LIMO Team 2010

%% check data orientation
STUDY        = varargin{1};
clsindx      = varargin{2};
subject_name = varargin{3};
data         = varargin{4};
weights      = varargin{5};
invweights   = varargin{6};
whichic      = varargin{7};


% Setting default method
if nargin == 7
    method = 'closercentroid';
elseif nargin == 7
    method = varargin{8};
end
clear varargin
nb_comp = size(data,1);
[n,p]=size(weights);

%% down to business
switch method
    
    %------
    case {'sum'} 
        % average power of ica weights on eletrodes
        W = sqrt(sum((invweights(:,whichic)).^2,1)); 
        combined_ica = mean(data(whichic,:,:).*repmat(W',[1,size(data,2),size(data,3)]),1);

    case ('max') 
        % take ica which has max variance
        icaact = weights*mean(data,3);
        varica = var(icaact(whichic,:), [], 2);      
        [tmp ind] = max(varica);
        combined_ica = data(ind,:,:);
        
    case ('maxvar') 
        % takes the ica which has the max mean variance accounted for
        icaact = weights*mean(data,3);
        for iComp = 1:length(whichic)
            comp = whichic(iComp);
            dataMinusIca = mean(data,3) - invweights(:,comp)*icaact(comp,:);
            varica(iComp) = mean(var(dataMinusIca, [], 2));
        end
        [tmp ind] = max(varica);
        combined_ica = data(ind,:,:);
        
    case ('closercentroid')
        [all_ickeepindx,ickeep_sub,datasetinfo_first] = limo_reducecluster(STUDY,clsindx);
        subjindxtmp = find(strcmp(ickeep_sub,subject_name));
        ickeepindx = all_ickeepindx(subjindxtmp);
        combined_ica = data(ickeepindx,:,:);
end


% figure;plot(mean(data,3)','--','LineWidth',1.5)
% hold on; plot(mean(combined_ica,3),'k','LineWidth',3); axis tight
