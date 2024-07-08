function copySoRoSimExamples(path)
% copySoRoSimExamples Copies the example files to a specified folder
%
% copySoRoSimExamples( path ) copies the SoRoSim example files to the folder specified by path
%
%
% Example:
%
%  mkdir 'C:\Temp\SoRoSim_examples'
%  cd 'C:\Temp\SoRoSim_examples'
%  copySoRoSimExamples('.')
%

% arguments
%     path (1,1) string {mustBeFolder}
% end

if nargin==0
    path = uigetdir;
end

soRoSimTopFolder = fileparts(mfilename("fullpath"));
exampleFolder = fullfile(soRoSimTopFolder, "Examples");

copyfile(exampleFolder, path, 'f')
disp(['Examples copied to ',path])

end