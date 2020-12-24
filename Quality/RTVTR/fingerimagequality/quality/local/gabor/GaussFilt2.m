function [ImageOut] = GaussFilt2(Image,Sigma),
% Filter image by a Gaussian window.
%
% Syntax:
%   [ImageOut] = GaussFilt2(Image,Sigma);
%
% Inputs:
%   Image   - input image
%   Sigma   - parameter Sigma of the Gaussian window
%
% Outputs:
%   GF      - output image
%
% Updates:
%   2011-08-15          Initial version
%   2012-04-26          Cleanup

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% Haiyun Xu (haiyun.xu@gmail.com)
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

    if nargin < 2, Sigma     = 3; end

    Nstd  = 3;
    Supp  = -round(Nstd*Sigma):round(Nstd*Sigma);
    [x,y] = meshgrid(Supp,Supp);
    Gauss = 1 ./ (2*pi*Sigma^2) .* exp(-(x.^2 + y.^2) / (2*Sigma^2));

    ImageOut = FFTFilt2(Gauss,Image);
end