function [fdaIQM, blockRotated, blockCropped, t, dft, dftAmp] = fda(block, orientation, v1sz, scanResolution)
% FDA Computes image quality measure (IQM) for the Frequency Domain Analysis of ridges and valleys.
% Returns fdaIQM by performing ridge-valley periodical signature analysis within
% a block of fingerprint image given as a parameter. 
%
% Syntax:
%   fdaIQM = fda(block, orientation, scanResolution)
%
% Inputs:
%   block           - square block of image (orientation block + border to fully cover rotated img)
%   orientation     - angle of the orientation line perpendicular to the ridge direction
%                     within the block [rad]
%   scanResolution  - scanner resolution [ppi]
%
% Outputs:
%   fdaIQM          - local quality score (of the block)
%
% Examples:
%   fdaIQM = fda([36 36], ang_in_deg, 500);

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

% sanity check: scanner resolution
if (scanResolution ~= 500)
    warning('Wrong scanner resolution - FDA script made for 500 ppi'); 
end

% sanity check: check block size
cBlock = size(block,1)/2; %square block
if fix(cBlock) ~= cBlock;
    warning('MATLAB:bsize', 'Wrong block size! Consider block with size of even number');
end

% rotate image to get the ridges horizontal using nearest-neighbor interpolation
blockRotated = imrotate(block, rad2deg(orientation+pi/2), 'nearest', 'crop');

% set x and y
xoff = v1sz(1)/2;
yoff = v1sz(2)/2;

% extract slanted block by cropping the rotated image: To encure that rotated image does
% not contain any invalid regions.
blockCropped = blockRotated(cBlock-(xoff-1):cBlock+xoff,cBlock-(yoff-1):cBlock+yoff);

% Signature
t = zeros(size(blockCropped,1),1); 

% Calculate mean
for x=1:size(t),
    t(x) = mean(blockCropped(x,:));
end

dft = fft(t');

% Amplitude, cutting out DC
dftAmp = abs(dft(1, 2:end)); 

[fMax,fMaxIndex] = max(dftAmp);
iqmDenom = sum( dftAmp(1, 1: floor(length(dftAmp)/2) ) );

if(fMaxIndex==1) 
    fdaIQM = 1;  
    return
end

fdaIQM = ( fMax + 0.3*(dftAmp(fMaxIndex-1) + dftAmp(fMaxIndex+1)) ) / iqmDenom;