function [GF] = GaborFilterShen(Theta,Freq,Sigma)
% Calculating Gabor filter using Shen's paper definition.
%
% Syntax:
%   [GF] = GaborFilterShen(Theta,Freq,Sigma);
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
%   2011-08-31          Initial version

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% Haiyun Xu (haiyun.xu@gmail.com)
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

if nargin < 3, Sigma = 4;      end
if nargin < 2, Freq  = 0.12;    end
if nargin < 1, Theta = pi/8;   end

[x,y] = meshgrid(-2.5*Sigma:2.5*Sigma,-2.5*Sigma:2.5*Sigma);

x1 = x*cos(Theta) + y*sin(Theta);
y1 = -x*sin(Theta) + y*cos(Theta);

GF = exp(-1/2 * (x1.^2 / Sigma^2 + y1.^2 / Sigma^2)) .* exp(1i*2*pi*Freq*x1);