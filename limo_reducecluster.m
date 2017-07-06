function [ickeepindx,ickeep_sub,datasetinfo_first] = limo_reducecluster(STUDY,clsindx,varargin)
% Given the STUDY with clusters and the index of the cluster to analyze.
% This function will reduce the cluster to have only one IC per subject by
% picking up the IC of the subject that lies closer to the centroid of the
% cluster in the joint space(from the preclustering step).As a result the
% function will return:
% ickeepindx          - Index in the ICs to keep from cluster's IC. The same
%                       order of the matrix of measures and STUDY
% ickeep_sub          - Subject names associated to the index in ickeepindx
% datasetinfo_first   - Index in STUDY.datasetinfo of the first dataset in
%                       the sets of data for each component. Useful to retrieve information of
%                       this specific subject. Since the subject information is present in every
%                       set associated to each component, only one is provided as output
% See also: 
%
% Author: Ramon Martinez-Cancino, SCCN, 2017
%
% Copyright (C) 2017  Ramon Martinez-Cancino, 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('limo_reducecluster() error: calling convention {''key'', value, ... } error'); return;
end;
% try g.cluster; catch, g.cluster =  [];         end;

% Preparing input
if ~isfield(STUDY.etc, 'preclust') && ~isfield(STUDY.etc.preclust, 'preclustparams') && isempty(STUDY.etc.preclust.preclustparams)
    error('Not enough information provided in STUDY structure');
end
if isempty(clsindx)
    display('std_cca_reducecluster: Index of cluster is a required input');
    return
end

if ~isempty(STUDY.cluster(clsindx).child)
    display('std_cca_reducecluster: Invalid cluster. Could not harvest Parent clusters');
    return
end

datasetinfo_tmp = [];
for iset = 1:size(STUDY.cluster(clsindx).sets,2)
    if ~isfield('cell',STUDY.design)
        % For new EEGLAB versions with no 'cell' field under STUDY.design
        subjcluster(iset)  = unique( {STUDY.datasetinfo(STUDY.cluster(clsindx).sets(iset)).subject});
    else
        % For old EEGLAB versions (backcompatibility)
        tmpstruct  = std_setcomps2cell(STUDY, STUDY.cluster(clsindx).sets(:,iset), STUDY.cluster(clsindx).comps(:,iset));
        cellinds   = [ tmpstruct.setinds{:} ];
        cells      = STUDY.design(STUDY.currentdesign).cell(cellinds);
        
        subjcluster(iset) = unique( {STUDY.datasetinfo([cells.dataset]).subject});
        datasetinfo_tmp   = cat(2,datasetinfo_tmp,cells(1).dataset);

    end
end
datasetinfo_tmp = STUDY.cluster(clsindx).sets;
uniquesubj      = unique(subjcluster);

% Measuring euclidian distance to centroids in the joint space
preclustdata        = STUDY.cluster(clsindx).preclust.preclustdata;
jointspace.centroid = mean(preclustdata,1);

ickeepindx = []; ickeep_sub= {}; datasetinfo_first = [];
 for isub = 1:length(uniquesubj)
     subindxtmp  = find(strcmp(subjcluster,uniquesubj{isub}));
     ickeep_sub  = cat(2,ickeep_sub,uniquesubj(isub));
     if length(subindxtmp) == 1
         ickeepindx     = cat(2,ickeepindx,subindxtmp);
         datasetinfo_first = cat(2,datasetinfo_first,datasetinfo_tmp(subindxtmp));
     else
         Sjoint = squareform(pdist([jointspace.centroid' preclustdata(subindxtmp,:)']','euclidean'));
         [~,ictmp]         = min(Sjoint(2:end,1));
         ickeepindx        = cat(2,ickeepindx,subindxtmp(ictmp));
         datasetinfo_first = cat(2,datasetinfo_first,datasetinfo_tmp(subindxtmp(ictmp)));
     end
 end