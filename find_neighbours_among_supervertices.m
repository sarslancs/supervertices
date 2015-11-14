function [ N ] = find_neighbours_among_supervertices( nSupervs, superLabels, corticalMask, faces )
%FIND_NEIGHBOURS_AMONG_SUPERVERTICES Generates the binary neighbourhood
%matrix for the supervertices
%	Returns N, the [nSupervs x nSupervs] neighbourhood matrix, for the
%	supervertices. This matrix is necessary to apply a spatially
%	constrained hierarchical clustering algorithm on the supervertices if
%	the merging is to be performed only between the adjacent supervertices.

% Neighbourhood matrix N
N = zeros(nSupervs,nSupervs);

% Mappers
[ map32to29, map29to32 ] = cortical_mappers(corticalMask);

for i = 1 : nSupervs    
    maps = map29to32(superLabels == i);
    idx = [cell2mat(arrayfun(@(x) find(faces(:,1) == x, 1,'first'), maps, 'UniformOutput', false));...
            cell2mat(arrayfun(@(x) find(faces(:,2) == x, 1,'first'), maps, 'UniformOutput', false));...
            cell2mat(arrayfun(@(x) find(faces(:,3) == x, 1,'first'), maps, 'UniformOutput', false))];
    
    nodes = unique(faces(idx,:));
    
    idx = [cell2mat(arrayfun(@(x) find(faces(:,1) == x, 1,'first'), nodes, 'UniformOutput', false));...
            cell2mat(arrayfun(@(x) find(faces(:,2) == x, 1,'first'), nodes, 'UniformOutput', false));...
            cell2mat(arrayfun(@(x) find(faces(:,3) == x, 1,'first'), nodes, 'UniformOutput', false))];        
    
    nodes = unique(faces(idx,:));
    nidx = unique(superLabels(nonzeros(map32to29(nodes))));
    nidx(nidx == i) = [];
    N(i,nidx) = 1;  
   
end

N = logical(N | N');



