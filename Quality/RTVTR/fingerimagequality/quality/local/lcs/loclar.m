function [lcs, blkrot, v2, v3, x, dt, ridval, ridmat, dtridmat, ridbad, valmat, dtvalmat, valbad] = loclar(blk, orang, v1sz, sc)
% LOCLAR Computes local clarity score (LCS) of ridges and valleys.
% Returns lcs [ 1(best) - 0(worst) ] by performing ridge-valley structure 
% analysis within a block of FP image given as a parameter.
%
% Syntax:
%   lcs = loclar(blk, orang, v1sz, sc)
%
% Inputs:
%   blk         - square block of image (orientation block + border to fully cover v1)
%   orang       - angle of the orientation line perpendicular to the ridge direction
%               - within the block [rad]
%   v1sz        - size of slanted square to extract from block [px] (recommended 32x16)
%   sc          - scanner resolution [ppi]
%
% Outputs:
%   lcs     	- local clarity score (of the block) [ 0(worst) - 1(best) ]
%
% Examples:
%   lcs = loclar([36 36], ang_in_deg, [32 16], 500);
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Master Thesis, Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany


Wrmax125 = 5; % max ridge for 125 ppi scanner
Wvmax125 = 5; % max valley for 125 ppi scanner

Wrmin = 3; % [px]
Wrmax = 10;
Wvmin = 2;
Wvmax = 10;

%%
cblk = size(blk,1)/2; %square block
if fix(cblk) ~= cblk;
    warning('MATLAB:bsize', 'Wrong block size! Consider block with size of even number');
end

% rotate image to get the ridges vertical
blkrot = imrotate(blk, rad2deg(orang), 'nearest', 'crop');

xoff = v1sz(1)/2;
yoff = v1sz(2)/2;

% extract slanted block by cropping the rotated image: To ensure that rotated image does
% not contain any invalid regions.
v2 = blkrot(cblk-(yoff-1):cblk+yoff,cblk-(xoff-1):cblk+xoff);

% average profile of v2: Compute average of each colum to get a projection of the grey
% values down the ridges.
v3 = mean(v2);

%% Linear regression using least square
% output = input * coefficients
% operator / "divide" the output by the input to get the linear coefficients
x = 1:length(v3);
% Append a column of ones before dividing to include an intercept, dt1 = [intercept
% coeficient]
dt1 = [ones(length(x),1) x'] \ v3';

%% Block segmentation into ridge and valley regions
dt = x*dt1(2) + dt1(1);
ridval = (v3 < dt)'; % ridges = 1, valleys = 0

%% Ridge-valley thickness
begrid = ridval(1); % begining with ridge?
change = xor(ridval, circshift(ridval, 1)); % find the bin change
change(1) = []; % there can't be change in 1. element
change = find(change == 1);    % find indices
if ~isempty(change) % changes found = ridge-val structure
    change1r = circshift(change,1); change1r(1) = 0;
    Wrv = change - change1r; % ridge and valley thickness
    if begrid
        Wr = Wrv(1:2:end); % odd indeces
        Wv = Wrv(2:2:end); % even indeces
    else
        Wv = Wrv(1:2:end); % odd indeces
        Wr = Wrv(2:2:end); % even indeces       
    end
    NWr = Wr / ((sc/125)*Wrmax125);
    NWv = Wv / ((sc/125)*Wvmax125);
    
    % normalized max/min
    NWrmin = Wrmin / ((sc/125)*Wrmax125);
    NWrmax = Wrmax / ((sc/125)*Wrmax125);
    NWvmin = Wvmin / ((sc/125)*Wrmax125);
    NWvmax = Wvmax / ((sc/125)*Wrmax125);
else % NOT ridge/valley structure, skip computation
    lcs = NaN;
    return
end

%% Clarity test
% NOTE: can be different strategy how to deal with out of limit ridge-valley thickness:
% NOTE: first and last region can be UNCOMPLETE -should be somhow excludided from the test
% 1: all should fall in (except first/last)
% 2: majority
% 3: mean/median of all
muNWr = mean(NWr);
muNWv = mean(NWv);
if muNWr >= NWrmin && muNWr <= NWrmax && muNWv >= NWvmin && muNWv <= NWvmax
    % ridge region
    ridmat = v2(:,ridval==1); % matix of ridge pxs (according to v3 and dt1)
    dtridmat = ones(size(ridmat))*diag(dt(ridval==1)); % dt1 tresh for each coloum of ridmat
    ridbad = ridmat >= dtridmat;
    %valley region
    valmat = v2(:,ridval==0); % matix of valley pxs (according to v3 and dt1)
    dtvalmat = ones(size(valmat))*diag(dt(ridval==0)); % dt1 tresh for each coloum of ridmat
    valbad = valmat < dtvalmat;
    
    alpha = mean(ridbad(:));
    beta = mean(valbad(:));
    lcs = 1 - ((alpha + beta)/2);
else
    lcs = 0;
end
