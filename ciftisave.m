function [ output_args ] = ciftisave(cifti,filename,caret7command)
%Save a CIFTI file as a GIFTI external binary and then convert it to CIFTI.
%Provided as part of the Human Connectome Project pipelines

tic
save(cifti,[filename '.gii'],'ExternalFileBinary')
toc

tic
unix([caret7command ' -cifti-convert -from-gifti-ext ' filename '.gii ' filename]);
toc

unix([' /bin/rm ' filename '.gii ' filename '.dat ']);

end

