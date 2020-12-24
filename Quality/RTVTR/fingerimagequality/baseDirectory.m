function [ rootDirectory ] = baseDirectory
%BASEDIRECTORY Tells the location of this file
%
% Syntax: 
%   rootDirectory = baseDirectory
%
% Inputs:
%   none
%
% Outputs:
%   rootDirectory - A string containing the path to this file
%
% Examples:
%   baseDirectory
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

rootDirectory = fileparts(mfilename('fullpath'));
end