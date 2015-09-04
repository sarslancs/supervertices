function [ dists ] = dist_func( XI, XJ, dg, R )
%DIST_FUNC Hybrid distance function that combines functional and spatial
%smilarity
%   [ DISTS ] = DIST_FUNC( XI, XJ, DG, R ) implements the distance function
%   D = sqrt(((dc / Nc) .^ 2) + ((dg / Ng) .^ 2)) where dc and DG 
%   correspond to the functional and spatial distance measures, 
%   respectively. Functional similarity is measured by the Pearson’s 
%   distance transformation. Spatial proximity is measured by the geodesic 
%   distance  along the cortical surface, approximated as the length of the 
%   shortest path between the nodes. Nc and Ng refer to the normalization 
%   factors. XI is the timeseries of the supervertex, wheras XJ is the set 
%   of timeseries of the cortical vertices in its search space.


% Internal parameters
Nc = 2; % Functional normalization factor
Ng = R; % Spatial normalization factor

if (size(dg,1) == 1)
    dg = dg';
end

XI = repmat(XI, size(XJ,1), 1);
dc = 1 - diag(tril(corr(XI', XJ'),0)); 
dists = sqrt(((dc / Nc) .^ 2) + ((dg / Ng) .^ 2)); 




