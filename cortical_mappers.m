function [ map32to29, map29to32 ] = cortical_mappers( corticalMask )
%CORTICAL_MAPPERS Maps vertices between 29K and 32K cortical structures
%   [ MAP32TO29, MAP29TO32 ] = CORTICAL_MAPPERS( CORTICALMASK ) returns the
%   mappers MAP32TO29, MAP29TO32 that work like a dictionary, telling you
%   the corresponding vertex in the other data structure. The HCP pipelines
%   represent the cortical surfaces as 32K standard meshes, in which the
%   medial wall vertices are not associated with any data source. For
%   example, the resting-state fMRI timeseries for the left hemisphere is
%   stored in a 29K x 1200 matrix, since no timeseries data is available
%   for the ~3K medial wall vertices. Moreover, there is no direct
%   correspondance between the 32K cortical vertices and their formation 
%   (or rows) in the 29K data matrix. For a given vertex ID, MAP32TO29
%   returns its associated row in the 29K data matrix, and vice versa.

rois = (1:length(corticalMask))';
drois = (1:length(corticalMask(corticalMask == 1)))';
rois(corticalMask == 1) = drois;
rois(corticalMask == 0) = 0;
map29to32 = arrayfun(@(x) find(rois == x, 1,'first'), ...
                    (1:length(corticalMask))', 'UniformOutput', false);
map29to32 = cell2mat(map29to32);

rois = find(corticalMask == 1);
map32to29 = corticalMask;
map32to29(map32to29 == 1) = 1:length(rois);



