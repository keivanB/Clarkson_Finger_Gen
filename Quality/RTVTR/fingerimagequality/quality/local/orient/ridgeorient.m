function orientang = ridgeorient(a, b, c)
%RIDGEORIENT Computes orientation angles
%    Computes angles of orientation lines perpendiular to the ridge direction
%
% Syntax:
%   orientang = ridgeorient(a, b, c)
%
% Inputs:
%   a         - param a of covariance matix [a c; c b]
%   b         - param b of covariance matix [a c; c b]
%   c         - param c of covariance matix [a c; c b]
%
% Outputs:
%   orientang - orientation angle of line perpendicular to the ridge flow

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany

    denom = sqrt(c.^2 + (a - b).^2) + eps;
    sin2theta = c./denom;
    cos2theta = (a-b)./denom;

    orientang = atan2(sin2theta,cos2theta)/2; 
end