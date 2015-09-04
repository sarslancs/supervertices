function [ superLabels, superSeries, superVertices, superIdx ]  = supervertex_clustering( ...
          timeseries, verticalCoords, corticalMask, superVertices, superIdx,  adj, R, geodesics )
%SUPERVERTEX_CLUSTERING Cluster cortical vertices into supervertices
%   SUPERVERTEX_CLUSTERING groups the cortical vertices into a predefined 
%   numbers of supervertices. 
% 
%   --------------------------------INPUT----------------------------------
%   TIMESERIES: Timeseries data that will be clustered 
% 
%   VERTICALCOORDS: x, y, z, coordinates of the cortical vertices 
% 
%   CORTICALMASK: A binary mask that separates the cortical vertices from
%   the medial wall.
% 
%   SUPERVERTICES: x, y, z coordinates of the initial seeds. 
%
%   SUPERIDX: The cortical ids of the intial seeds 
% 
%   ADJ: A binary matrix of the adjacent vertices in the cortical mesh    
%
%   GEODESICS = A cell structure of size 2 x 1, which stores the geodesic
%   distance of their cortical vertices and their neighbours on the 
%   cortical mesh.
% 
%   R: Range of the local search space in mm
% 
%   --------------------------------OUTPUT--------------------------------
%   SUPERLABELS: [29K x 1] label vector, showing the cluster membership of 
%   the cortical vertices. Each label is within [1, nSupervs].
% 
%   SUPERSERIES: [nSupervs x d] timeseries matrix, in which each timeseries
%   is associated with a supervertex.
% 
%   SUPERVERTICES: [nSupervs x 3] matrix. x, y, z coordinates of the 
%   supervertices.
% 
%   SUPERIDX: [nSupervs x 1] vector. Cortical IDs of the supervertices.
% 
%   ------------------------------ALGORITHM------------------------------- 
%   A k-means alike algorithm is proposed that parcellates the rs-fMRI
%   data into sub-regions with respect to their functional correlations
%   and spatial proximity. It is different from the classical k-means in 
%   two aspects: 1) The search space is limited to reduce the number of  
%   distance calculations and 2) the distance function is capable of 
%   matching highly correlated vertices with each other while enforcing 
%   spatial contiguity within a cluster. 

% Internal parameters
endThr      = 10;  % If the number of vertices that have been assigned 
                    % to another cluster in the previous iteration is less 
                    % than ENDTHR, then STOP. In the paper, it is set to 1 
                    % for the finest possible clusters, but its impact on 
                    % the final parcellations is tiny. 
changed     = 1; % Has anything changed between iterations or not?
verbose     = 1; % Write output on console?

% Mappers
[ map32to29, map29to32 ] = cortical_mappers(corticalMask);

% The number of supervertices
nSupervs = length(superVertices);

superSeries = timeseries(map32to29(superIdx),:);
superDists = Inf(length(timeseries),1);
superLabels = zeros(length(timeseries),1);

% Limit search space to reduce computational cost
geodesicsLocal = get_vertices_within_search_space(superVertices, superIdx, geodesics, R);

idsLocal = geodesicsLocal{1};
distsLocal = geodesicsLocal{2};

newSuperVertices = zeros(size(superVertices));
newSuperSeries = zeros(size(superSeries));
oldSuperLabels = superLabels;

fprintf('\nSupervertices are being computed');
while (changed)  
    
    % Print progress on the console
    if verbose == 1
        fprintf('.');
    end
    
    % Real action starts from here
    for i = 1 : nSupervs
       indx = map32to29(idsLocal{i}');  
       neighTimeseries = timeseries(indx,:);  
       distances = dist_func(superSeries(i,:), neighTimeseries, distsLocal{i}, R); 
       di = superDists(indx);
       Li = superLabels(indx);
       compared = distances < di;
       di(compared == 1) = distances(compared == 1);
       Li(compared == 1) = i;
       superDists(indx) = di;
       superLabels(indx) = Li;    
    end
     
    % Update cluster centers and bunch time series        
    for i = 1 : nSupervs
        idx = find(superLabels == i);
        if isempty(idx)
            idx = map32to29(superIdx(i));
            superLabels(idx) = i;
        end
        bunch = timeseries(idx,:);
        newSuperSeries(i,:) = mean(bunch,1);
        mapped = map29to32(idx);        
        newSuperVertices(i, :) = mean(verticalCoords(mapped,:), 1);
    end
    
    % If the program steps into this, it is because one or more 
    % supervertices have been undertaken by the others. This happens 
    % when R is too small compared to nSupervs. My solution is to locate
    % these supervertices and keep them as singleton clusters.
    if sum(superLabels == 0) > 0
        zs = map29to32(superLabels == 0);       
        for i = 1 : length(zs)
            id = zs(i);
            ids = nonzeros(map32to29(adj(id,:)>0));
            [ uniqs, ~ ] = count_unique_elements(nonzeros(superLabels(ids)));
            superLabels(map32to29(id)) = uniqs(1);
        end        
    end
        
    % The recently computed supervertex centroids are being matched with the
    % closest vertices on the cortical surface.
    temp = verticalCoords;
    temp(corticalMask == 0,:) = Inf; % Wall vertices are assigned to Inf.
    superIdx = knnsearch(temp, newSuperVertices);
        
    newSuperVertices = verticalCoords(superIdx, :);
    
    % Getting ready for the next iteration
    geodesicsLocal = get_vertices_within_search_space(newSuperVertices, superIdx, geodesics, R);
    idsLocal = geodesicsLocal{1};
    distsLocal = geodesicsLocal{2};
    
    superSeries = newSuperSeries;
    superVertices = newSuperVertices;    

    % Calculate if any vertex has changed its clsuter 
    if sum(oldSuperLabels ~= superLabels) < endThr
        changed = 0;
    else
        oldSuperLabels = superLabels;
    end
        
end

fprintf('\nDone.\n');

