function [ geodesics ] = compute_geodesics_wb_command( subjectID, inFolder, hem, corticalMask, R, verbose )
%COMPUTE_GEODESICS Compute geoedesic distances within the the given range
%   [ GEODESICS ] = COMPUTE_GEODESICS_WB_COMMAND
%                   ( SUBJECTID, INFOLDER, HEM, CORTICALMASK, R, VERBOSE ) 
%   returns a cell structure of size 2 x 1, which stores the geodesic 
%   distances of each cortical vertex and their closest neighbours on the 
%   cortical mesh. The distances are computed only for the cortical area 
%   masked with CORTICALMASK and within the local neighbourhood of each 
%   vertex, specified by R mm. SUBJECTID specifies the id of the subject to
%   be processed. INFOLDER is the location of the surface file, which will
%   be passed as a paramater to wb_command. HEM is the hemisphere.
%
%   Set VERBOSE = 0 in order to avoid the progress information being  
%   printed on the console
%
%   DIFFERENCE FROM COMPUTE_GEODESIC: This function directly makes use of
%   the native wb_command utility, distributed as part of the HCP Workbench,
%   available at http://www.humanconnectome.org/software/get-connectome-workbench.html
%   In order to use this function, the "wb_command" variable below should 
%   be set to the path the wb_command executable has been extracted to. 
%    
%   CAUTION: Similar to compute_geodesic, geodesic distances are not 
%   computed between all pairs at once. Instead, it computes the geodesic
%   distance from the specified vertex to all others within the given range 
%   R. The result is a single column metric file, which is then loaded into
%   the workspace and stored in a cell structure. Due to intensive I/O and
%   the dependency to an external program, the process could take up
%   to a few hours to finish for the entire cortex, depending on your
%   computational power.


ext = 'midthickness.32k_fs_LR.surf.gii';
wb_command = '../Workbench/exe_linux64/wb_command';
% Please set this to the path "wb_command" exe has been extracted to. 

out = [inFolder 'out/'];

if ~isdir(out)
    mkdir(out);
end

surface = [inFolder subjectID '.' hem '.' ext];

maskedIds = find(corticalMask == 1);
localIdx = cell(length(corticalMask), 1);
geodesicDists = cell(length(corticalMask), 1); 
   
nMsg = 0; % For console output only
for i = 1 : length(maskedIds)
    id = maskedIds(i);
    
    unix( [wb_command ' -surface-geodesic-distance ' surface ' '   ...
           num2str(id-1) ' ' out hem '_geodesics -limit ' num2str(R) ]);

    g = gifti([out hem '_geodesics']);
    dist = g.cdata;
    dist(corticalMask == 0) = -1;
    dist(dist == -1) = Inf;
    [sortedDist,sortedIds] = sort(dist);
    sortedIds(sortedDist > R) = [];
    sortedDist(sortedDist > R) = [];
    localIdx{id} = sortedIds;
    geodesicDists{id} = sortedDist;
    
    % Print progress on the console
    if nargin == 5 || verbose == 1
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

