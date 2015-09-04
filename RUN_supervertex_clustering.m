%% Set Path
% addpath('../lib/mesh/iso2mesh'); % http://iso2mesh.sourceforge.net/
% addpath('../lib/@gifti'); % http://www.artefact.tk/software/matlab/gifti/

%% Set parameters
subjectID       = '100307'; % Subject ID
hem             = 'L'; % Which hemisphere?
R               = 20; % Local search space in mm
nSupervs        = 1000; % Number of supervertices
keepRatio       = lookup_ratio_for_sampling(nSupervs, hem); % Keep ratio for sampling
saveOutput      = 0; % Save the generated output?  

%% Set input/output directories
readFrom = ['/data/HCP100/' subjectID];
% Please set the directory the files will be read from, if not specified,
% the output directory will be set automatically.

writeTo = [readFrom 'out/'];

if ~isdir(writeTo)
    mkdir(writeTo);
end
     
%% Load surface files
fileName = [subjectID '.' hem '.sphere.32k_fs_LR.surf.gii'];
g = gifti([readFrom fileName]);  
vSphere = g.vertices;
fSphere = g.faces;

fileName = [subjectID '.' hem '.midthickness.32k_fs_LR.surf.gii'];
g = gifti([readFrom fileName]);  
vGray = g.vertices;
fGray = g.faces;

fileName = [subjectID '.' hem '.atlasroi.32k_fs_LR.shape.gii'];
g = gifti([readFrom fileName]);  
corticalMask = g.cdata;

%% Compute the adjacency neighbourhood matrix
[ adj, adjWeighted ] = compute_vertex_nhood( vGray, fGray ); 

%% Load data
geodesics = compute_geodesics( adjWeighted, corticalMask, R );

fileName = [readFrom 'rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'];
dtseries = read_cifti(fileName);
timeseries = get_cortical_timeseries( dtseries, hem );

%% Compute initial seed vertices
[ seedCoors, seedIdx ] = compute_initial_seeds( vSphere, fSphere, vGray, corticalMask, keepRatio );
assert(nSupervs == length(seedCoors));

%% Compute supervertices
[ superLabels, superBunches, superVertices, superIdx ] = supervertex_clustering(timeseries, vGray, corticalMask, seedCoors, seedIdx, adj, R, geodesics);

%% Compute the adjacency matrix for supervertices
N = find_neighbours_among_supervertices(nSupervs, superLabels, corticalMask, fGray);            

%% Save output data structures for further analysis (e.g. the second clustering stage)
if saveOutput
    save([writeTo 'superBunches_n' num2str(nSupervs) '_' hem '.mat'],'superBunches');
    save([writeTo 'superVertices_n' num2str(nSupervs) '_' hem '.mat'],'superVertices');
    save([writeTo 'superLabels_n' num2str(nSupervs) '_' hem '.mat'],'superLabels');
    save([writeTo 'superIdx_n' num2str(nSupervs) '_' hem '.mat'],'superIdx');
    save([writeTo 'N_n' num2str(nSupervs) '_' hem '.mat'], 'N');
    disp('All generated files have been saved...');
end


           
                                                                                                                     
                          
