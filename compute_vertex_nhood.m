function [ adj, adj_weighted ] = compute_vertex_nhood( verticalCoords, meshFaces )
%COMPUTE_VERTEX_NEIGHBOURHOOD Compute the adjacency matrices mapping the 
%neighbour vertices and their Eucledian distance 
%   [ ADJ, ADJ_WEIGHTED ] = COMPUTE_VERTEX_NEIGHBOURHOOD( VERTICES, FACES) 
%   returns two, n x n sparse matrices. ADJ is a binary matrix of the 
%   adjacent vertices in the cortical mesh, whereas ADJ_WEIGHTED is a 
%   weight matrix of the distances between the adjacent vertices. The
%   cortical mesh is given its cooridantes VERTICALCOORDS and and faces 
%   MESHFACES.


nVertices = length(verticalCoords);
dynamicDists = zeros(nVertices*8,3);
ptrStart = 1;
ptrEnd = 0;

for i = 1 : nVertices;
    idx = [find(meshFaces(:,1)==i); ... 
           find(meshFaces(:,2)==i); ...
           find(meshFaces(:,3)==i)];
    
     neighs = unique(meshFaces(idx,:));
     dists = pdist2(verticalCoords(i,:),verticalCoords(neighs,:));
     ptrEnd = ptrEnd + length(neighs);
     dynamicDists(ptrStart:ptrEnd,1) = i;
     dynamicDists(ptrStart:ptrEnd,2) = neighs;  
     dynamicDists(ptrStart:ptrEnd,3) = dists; 
     ptrStart = ptrStart + length(neighs);
end

dynamicDists(dynamicDists(:,1) == 0,:) = [];
adj = sparse(dynamicDists(:,1),dynamicDists(:,2),ones(length(dynamicDists),1),nVertices,nVertices);
adj_weighted = sparse(dynamicDists(:,1),dynamicDists(:,2),dynamicDists(:,3),nVertices,nVertices);
   

