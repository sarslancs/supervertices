function [ cdata ] = read_cifti( filename, wb_command )
%READ_CIFTI Import dscalar.nii files into Matlab
%   CDATA = READ_CIFTI(FILENAME, WB_COMMAND) reads the file given as
%   FILENAME using the native wb_command program, the path of which is
%   given either as a string in WB_COMMAND or specified below as default.

if nargin == 1
    wb_command = '../Workbench/exe_linux64/wb_command';
    % Please set this to the path "wb_command" exe has been extracted to. 
end

cii = ciftiopen(filename, wb_command);  
cdata = cii.cdata;

end

