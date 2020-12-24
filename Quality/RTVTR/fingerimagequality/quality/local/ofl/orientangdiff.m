function loq = orientangdiff(blk)
%ORIENTANGDIFF Compute absolute differencies in orientation angles
%
% Syntax:
%   loq = orientangdiff(blk)
%
% Inputs:
%   blk        - image block
%
% Outputs:
%   loq        - orientation angle in radians
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Master Thesis, Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany

if isa(blk, 'struct')
    blk = blk.data;
end

blk = blk(:);
bsz = size(blk,1);
bci = bsz/2 + .5;

if fix(bci) ~= bci
    error('Wrong block size')
end

bcent = blk(round(bci));
blk(bci) = []; % exclude center from block
loq = sum(abs(bcent - blk)) / (bsz - 1);
