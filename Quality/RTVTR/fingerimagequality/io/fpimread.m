function gsim = fpimread(pathstr)
%FPIREAD Read and convert image it into grayscale 8bpp
%
% Syntax:
%   gsim = fpiread(pathstr)
%
% Inputs:
%   pathstr - string with the path to image file
%
% Outputs:
%   gsim    - matrix of grayscale image file
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Master Thesis, Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany

gsim = imread(pathstr);
if size(gsim,3) > 1
    gsim = rgb2gray(gsim);
end 