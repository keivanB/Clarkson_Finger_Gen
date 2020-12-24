function oclNOTISO = ocl(a, b ,c)
%OCL Computes orientation certainity level (OCL) score
% Computes oclISO score according to ISO/IEC TR 29794-4 [ 1(worst) - 0(best) ] 
% from the covariance matix [a c; c b]
%
% Syntax:
%   oclv = ocl(a, b ,c)
%
% Inputs:
%   a         - param a of covariance matix [a c; c b]
%   b         - param b of covariance matix [a c; c b]
%   c         - param c of covariance matix [a c; c b]
%
% Outputs:
%   oclNOTISO - orientation certainity value [ 0(worst) - 1(best) ]

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Master Thesis, Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany

    eigvmax = ( (a + b) + sqrt((a-b).^2 + 4.*c.^2) )/2;
    eigvmin = ( (a + b) - sqrt((a-b).^2 + 4.*c.^2) )/2;

    oclISO = eigvmin / eigvmax; % [ 1(worst) - 0(best) ]
    oclNOTISO = 1 - oclISO; % [ 0(worst) - 1(best) ]
end