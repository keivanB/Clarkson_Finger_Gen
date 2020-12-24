function [ratios, blockRotated, blockCropped, v3, x, dt1, dt, ridval, change] = rvu(block, orientation, v1sz, scanResolution)
% RVU computes ridge/valley ratios for a given block.
%
% Syntax:           - ratios = rvu(block, orientation, v1sz, scanResolution)
%
% Inputs:
%   block           - square block of image (orientation block + border to fully cover rotated img)
%   orientation     - angle of the orientation line perpendicular to the ridge direction
%                     within the block [rad]
%   v1sz            - size of slanted square to extract from block [px] (recommended 32x16)
%   scanResolution  - scanner resolution [ppi]
%
% Outputs:
%   ratios          - local ratios (ridge/valley) of the ridge valley structure
%
% Examples:
%   ratios = rvu([36 36], ang_in_deg, [32 16], 500);
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% Vladimir Smida

% sanity check: scanner resolution
if (scanResolution ~= 500)
    warning('Wrong scanner resolution - RVU script made for 500 ppi'); 
end

% sanity check: check block size
cBlock = size(block,1)/2; %square block
if fix(cBlock) ~= cBlock;
    warning('MATLAB:bsize', 'Wrong block size! Consider block with size of even number');
end

% rotate image to get the ridges vertical
blockRotated = imrotate(block, rad2deg(orientation), 'nearest', 'crop');

% set x and y
xoff = v1sz(1)/2;
yoff = v1sz(2)/2;

% extract slanted block by cropping the rotated image: To ensure that rotated image does
% not contain any invalid regions.
blockCropped = blockRotated(cBlock-(yoff-1):cBlock+yoff,cBlock-(xoff-1):cBlock+xoff); % v2

% average profile of blockCropped: Compute average of each colum to get a projection of the grey
% values down the ridges.
v3 = mean(blockCropped);

%% Linear regression using least square
x = 1:length(v3);
dt1 = [ones(length(x),1) x'] \ v3';

%% Block segmentation into ridge and valley regions
dt = x*dt1(2) + dt1(1);

ridval = (v3 < dt)'; % ridges = 1, valleys = 0

%% Ridge-valley thickness
change = xor(ridval,circshift(ridval,1)); % find the bin change
change(1) = []; % there can't be change in 1. element 
changeIndex = find(change == 1);    % find indices with changes

ratios = [];

if ~isempty(changeIndex) % changes found = ridge-val structure
    % non complete ridges/valleys are removed from ridval and changeIndex
    ridvalComplete = ridval(changeIndex(1)+1:changeIndex(end));
    changeIndexComplete = changeIndex - changeIndex(1);
    changeIndexComplete(1) = []; % removing first value
    
    if isempty(ridvalComplete) 
        return
    end;
    
    begrid = ridvalComplete(1); % begining with ridge?
        
    % changeIndex now represents the change values...
    changeIndexComplete(end:-1:2) = changeIndexComplete(end:-1:2)-changeIndexComplete(end-1:-1:1);
    
    for m=1:length(changeIndexComplete)-1,
        ratios(m) = changeIndexComplete(m)/changeIndexComplete(m+1);
    end;
    
    ratios(begrid+1:2:end) = 1 ./ ratios(begrid+1:2:end);
else % NOT ridge/valley structure, skip computation
    return
end

