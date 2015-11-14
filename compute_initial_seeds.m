function [ seeds, idx ] = compute_initial_seeds( vSphere, fSphere, vGray, corticalMask, ratio )
%COMPUTE_INITIAL_SEEDS Computes initial seeds for clustering
%   [ SEEDS, IDX ] = COMPUTE_INITIAL_SEEDS( VSPHERE, FSPHERE, VGRAY, CORTICALMASK, RATIO )
%   conputes the initial seeds using mesh resampling and knn search. The
%   spherical mesh represented by its vertices (VSPHERE) and faces
%   (FSPHERE) is uniformly subsampled into some percantage of the original
%   mesh, indicated by the decimation rate RATIO. The resulting vertices are 
%   then masked with CORTICALMASK to exlude the vertices with no timeseries 
%   data. Returns the coordinates of the seed vertices (SEEDS) and their
%   ids (IDX) on the cortical surface.
%
%   CAUTION: To run this function you have to add the iso2mesh library to 
%   the Path. Please find the library at: http://iso2mesh.sourceforge.net/
   
[v,f]=meshresample(vSphere,fSphere,ratio);
idx = knnsearch(vSphere,v);
elim = find(corticalMask == 0);
neighs = arrayfun(@(x) find(idx == x, 1,'first'), elim, 'UniformOutput', false);
A = cell2mat(neighs);
idx(A) = [];
seeds = vGray(idx,:);
    


