function [GF] = GaborFilterCx(Theta,Freq,Sigma)
% Calculating Gabor filter.
%
% Syntax:
%   [GF] = GaborFilterCx(Theta,Freq,Sigma);
%
% Inputs:
%   Theta   - parameter angle of the Gabor filter
%   Freq    - parameter Frequency of the Gabor filter
%   Sigma   - parameter Sigma of the Gabor filter
%
% Outputs:
%   GF      - Gabor filter
%
% Updates:
%   2011-08-15          Initial version
%   2012-04-26          Cleaning of code

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% Haiyun Xu (haiyun.xu@gmail.com)
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

    if nargin < 3, Sigma = 6;      end
    if nargin < 2, Freq  = 0.1;    end
    if nargin < 1, Theta = pi/2;   end

    [x,y] = meshgrid(-2.5*Sigma:2.5*Sigma, -2.5*Sigma:2.5*Sigma);

    x1 = x*sin(Theta) + y*cos(Theta);
    y1 = x*cos(Theta) - y*sin(Theta);

    GF = exp(-1/2 * (x1.^2 / Sigma^2 + y1.^2 / Sigma^2)) .* exp(1i*2*pi*Freq*x1);
end