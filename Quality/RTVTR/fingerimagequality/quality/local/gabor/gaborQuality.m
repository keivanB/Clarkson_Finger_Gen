function [GaborStdQuality, imGauss, GF, FP_filt, GaborAngle, response, GaborStd] = gaborQuality(im, freq, sigma, angleNum)
% By applying Gabor filter to calculate the fingerprint quality.
% Implements algorithm by 
%   Olsen, M.; Xu, H. & Busch, C. 
%   Gabor filters as candidate quality measure for NFIQ 2.0 
%   Biometrics (ICB), 2012 5th IAPR International Conference on, 2012, 158-163
%
% Syntax:
%   GaborStdQuality = gaborQuality(im, freq, sigma, angleNum)
%
% Inputs:
%   im       - fingerprint image
%   freq     - parameter Frequency of the Gabor filter
%   sigma    - parameter Sigma of the Gabor filter
%   angleNum - number of angles in the filter bank
%
% Outputs:
%   GaborStdQuality - standard deviation of the Gabor responses at each point of
%   fingerprint
%
% Updates:
%   2011-08-15  Initial version
%   2011-09-12  Added configurable number of angles (martin.olsen@cased.de)%    

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% Haiyun Xu (haiyun.xu@gmail.com)
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

im = double(im)/255;

% remove edge effect of gabor filtering.
imGauss = im - GaussFilt2(im,1);

[Mfp,Nfp] = size(imGauss);

GaborAngle = zeros(Mfp, Nfp, angleNum);
GF = zeros(5*sigma+1, 5*sigma+1, angleNum);
FP_filt = zeros(Mfp, Nfp, angleNum);

if nargout > 2
    response = zeros(Mfp, Nfp, angleNum);
end

% Calculate Gabor response for each of the angles.
i=1;
for Theta = pi*(0:angleNum-1)/angleNum,
    % (1) Get Gabor filter with certain parameters
    GF(:,:,i) = GaborFilterCx(Theta,freq,sigma);
    
    % (2) Apply Gabor filter
    FP_filt(:,:,i) = filter2(GF(:,:,i), imGauss);
    
    % (3) Convolve Gabor response with Gaussian filter. 
    GaborAngle(:,:,i) = GaussFilt2(abs(FP_filt(:,:,i)),4);
    if nargout > 1
        response(:,:,i) = GaussFilt2(FP_filt(:,:,i), 4); 
    end %FP_filt
    i=i+1;
end

% Standard deviation of filter responses.
GaborStd = std(GaborAngle, 0, 3);

GaborStdQuality = sum(GaborStd(:))/Mfp/Nfp; 
