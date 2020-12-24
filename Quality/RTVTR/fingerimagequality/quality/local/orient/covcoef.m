function [a, b, c] = covcoef(blk)
%COVCOEF Computes covariance coefficients of grey level gradients
%   Computes coefficients of covariance matrix [a c; c b] of grey level 
%   gradients inside a block specified as a parameter
%
% Syntax:
%   [a b c] = covcoef(blk)
%
% Inputs:
%   blk       - block of the image
%
% Outputs:
%   a         - param a of covariance matix [a c; c b]
%   b         - param b of covariance matix [a c; c b]
%   c         - param c of covariance matix [a c; c b]

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Master Thesis, Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany

    [fx, fy] = gradient(double(blk));
    a = mean(fx(:).^2);
    b = mean(fy(:).^2);
    c = fx.*fy; c = mean(c(:));
end