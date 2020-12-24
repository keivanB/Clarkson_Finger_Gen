%STARTUP Sets up the MATLAB environment for fingerprint quality development
%   Checks for license dependencies, performs a selftest, and 
%   Add the directory of this file to MATLAB userpath and this function
%   will run automatically upon starting MATLAB.
%
% Syntax: 
%   startup
%
% Inputs:
%   none
%
% Outputs:
%   none
%
% Examples:
%   startup
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

close all;
clear all;
clc;

fprintf('Startup script for finger image quality codebase (version %s).\n\n', '1.0');
fprintf('See http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055\n')
fprintf('\tand http://nislab.no/software/fingerprintquality/\n')
fprintf('\tfor algorithm descriptions and evaluation.\n\n')
tStart = tic;

%% Add paths
rootDirectory = environmentPaths();

%% Set prefs
warning on all;
warning on backtrace;

%% Check if licenses for required toolboxes are available
toolboxDependencies = { ...
    'MATLAB', ...
    'statistics_toolbox', ...
    'image_toolbox', ...
    'signal_toolbox'};

hasLicense = zeros(length(toolboxDependencies), 1);
for ind = 1:length(toolboxDependencies)
    tb = toolboxDependencies{ind};
    hasLicense(ind) = license('checkout', tb);
end
clear ind tb;

if all(hasLicense)
    display('Licenses for all dependencies were found.');
else
    warning(['Not all licenses for required toolboxes were found. ', ...
        'Code may not run properly or at all.']);
    fprintf('Did not find following toolboxes:\n');
    fprintf('\t%s\n', (toolboxDependencies{~hasLicense}));
end
clear hasLicense toolboxDependencies;

%% Reset the random number generator
rng('default');

%% Go to codebase root
cd(rootDirectory);
tEnd = toc(tStart);
fprintf('Startup completed in %0.2f seconds.\n', tEnd);

clear tStart tEnd;
fprintf('Current working directory: %s\n', pwd);
