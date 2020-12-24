function [ baseDir ] = environmentPaths
%ENVIRONMENTPATHS Adds all necessary paths for the codebase to MATLAB
%
% Syntax: 
%   baseDir = environmentPaths
%
% Inputs:
%   none
%
% Outputs:
%   baseDir - the base directory for the environment
%
% Examples:
%   environmentPaths
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

%% Environment paths
fprintf('Setting up user paths ... ')
baseDir = baseDirectory;
P = genpath(baseDir);
addpath(P);
fprintf('done!\n');

%% Clean up
% Remove paths starting with '.' or '~'
fprintf('Cleaning up paths ... ');
removePath = [ ];
[string, remString] = strtok(P, pathsep);

while ~isempty(string);
    if ~isempty(strfind(string, [filesep '.'])) ...
    || ~isempty(strfind(string, [filesep '~']))
        removePath = cat(2, removePath, pathsep, string);
    end
    [string, remString] = strtok(remString, pathsep); %#ok
end

if ~isempty(removePath)
    rmpath(removePath);
end

fprintf('done!\n');

end