function [ geodesics ] = compute_geodesics( adjWeighted, corticalMask, R, verbose )
%COMPUTE_GEODESICS Compute geoedesic distances within the given range
%   [ GEODESICS ] = COMPUTE_GEODESICS( ADJWEIGHTED, CORTICALMASK, R, VERBOSE )
%   returns a cell structure of size 2 x 1, which stores the geodesic
%   distances of the cortical vertices to their neighbours across the 
%   cortical mesh. The distances are computed only for the cortical area
%   masked with CORTICALMASK and within the local neighbourhood of each
%   vertex, specified by R in mm. 
%
%   Set VERBOSE = 0 in order to avoid the progress information being  
%   printed on the console
%    
%   CAUTION: Due to memory issues, geodesic distance is not computed
%   between all pairs at once using graphallshortestpaths. Instead, 
%   graphshortestpath has been used within a for loop, which could take up
%   to a few hours to finish for the entire cortex, depending on your
%   computational power.

maskedIds = find(corticalMask == 1);
localIdx = cell(length(corticalMask), 1);
geodesicDists = cell(length(corticalMask), 1); 

nMsg = 0; % For console output only
for i = 1 : length(maskedIds)
    id = maskedIds(i);
    [dist, ~] = graphshortestpath(adjWeighted, id, maskedIds, 'Directed', false);  
    [sortedDist,sortedIds] = sort(dist);    
    sortedIds(sortedDist > R) = [];
    sortedDist(sortedDist > R) = [];
    localIds = maskedIds(sortedIds);
    localIdx{id} = localIds;
    geodesicDists{id} = sortedDist';
    
    % Print progress on the console
    if nargin == 3 || verbose == 1
        msg = [sprintf('\nProcessing node %d of %d\n', i, length(maskedIds)) ...
                      repmat('|', 1, ceil(i/100))];
        fprintf(repmat('\b',1, nMsg));
        fprintf(msg);
        nMsg = numel(msg);
    end

end

geodesics = cell(2,1);
geodesics{1} = localIdx;
geodesics{2} = geodesicDists;

