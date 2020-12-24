function [GaborShenQuality, GF, GaborAngle, gaborBlkStd, gaborFore, gaborPoor] = gaborShenQuality(im, Tb, Tq, angleNum, freq, sigma, blksze)
% By applying Gabor filter in Shen's paper to calculate the fingerprint quality.
% Implements algorithm by 
%    Shen, L.; Kot, A. & Koo, W. 
%    Quality Measures of Fingerprint Images 
%    PROC. AVBPA, SPRINGER LNCS-2091, 2001, 266-271
%
% Syntax:
%   GaborStdQuality = gaborShen(im, Tb, Tq, angleNum);
%
% Inputs:
%   im       - fingerprint image
%   Tb       - segmentation parameter for the background threshold
%   Tq       - good/poor quality threshold
%   angleNum - number of angles in gabor filter bank - default is 8
%
% Outputs:
%   GaborShenQuality - standard deviation of the Gabor responses at each point of
%   fingerprint
%
% Updates:
%   2011-08-31          Initial version
%   2011-09-13  Added configurable number of angles
%   2015-04-29  Using the blockwise standard deviation as the block value

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% Haiyun Xu (haiyun.xu@gmail.com)
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

im = double(im)/255;

[Mfp,Nfp] = size(im);

GF = zeros(5*sigma+1, 5*sigma+1, angleNum);
GaborAngle = zeros(Mfp,Nfp,angleNum);

% Calculate response of gabor filters
i=1;
for Theta = pi*(0:angleNum-1)/angleNum,
    % (1) Get Gabor filter with certain parameters
    GF(:,:,i) = GaborFilterShen(Theta,freq,sigma);
    % (2) Apply Gabor filter
    GaborAngle(:,:,i) = abs(filter2(GF(:,:,i),im));
    i=i+1;
end

% Standard deviation per block of filter responses
fun = @(block_struct) std(block_struct.data(:));
gaborBlkStd = blockproc(GaborAngle, [blksze blksze], fun);
gaborFore = gaborBlkStd > Tb; % Foreground blocks
gaborPoor = gaborBlkStd > Tb & gaborBlkStd < Tq; % Foreground blocks with low q

% 0 quality if no foreground found
if (sum(gaborFore(:)) == 0)
    GaborShenQuality = 0;
else
    GaborShenQuality = 1-sum(gaborPoor(:))/sum(gaborFore(:));
end
