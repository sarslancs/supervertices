function [ geodesicsLocal ] = get_vertices_within_search_space( superVertices, superIdx, geodesics, R )
%GET_VERTICES_WITHIN_SEARCH_SPACE Get the vertices within the search space
%of the supervertices
%   For each supervertex in SUPERVERTICES, computes the vertices within 
%   a search space R, centred on the supervertex. The ids of the vertices 
%   and their geodesic distance are returned in a cell structure 
%   GEODESICSLOCAL of size [2 x 1]. 

neighs = geodesics{1};
neigh_idx = neighs(superIdx);
dist_idx = geodesics{2};
dist_idx = dist_idx(superIdx);

n = size(superVertices, 1);

for j = 1 : n
    neighs = neigh_idx{j};
    dists = dist_idx{j};
    elim = dists > R;
    dists(elim == 1) = [];
    neighs(elim == 1) = [];
    neigh_idx{j} = neighs;
    dist_idx{j} = dists;
end

geodesicsLocal = cell(2,1);
geodesicsLocal{1} = neigh_idx;
geodesicsLocal{2} = dist_idx;



