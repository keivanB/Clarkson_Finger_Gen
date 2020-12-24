function [quality, imagepadded, imagewindowed, FS, pow, powpol, FOI, ...
    FOIPower, radialpowerspectrum, radialpowerspectrumbinned, fsmin, ...
    fsmax] = rps(image, rad, theta, fmin, fmax)
% RADIALPOWERSPECTRUM computes the quality of a finger image based on
% powerspectrum
%
% Inputs
%  image  - finge image
%  rad    - maximum number of bins between fmin and fmax in the polar
%           representation to be considered. If there are less possible 
%           bins than rad due to image input size then only the possible 
%           number of bins are included.
%  theta  - Angular resolution of polar spectrum. Sampling angles are in
%           the range 0:(pi/theta-1):pi
%  fmin   - frequency of innermost boundary of annular band in cycles/pixel from DC
%  fmax   - frequency of outermost boundary of annular band in cycles/pixel
%           from DC
%

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

    % Defaults for fmin and fmax are based on the expect ridge/valley width
    if ~exist('fmin', 'var') || isempty(fmin) || fmin < 0
        fmin = 0.01;
    end

    if ~exist('fmax', 'var') || isempty(fmax) || fmax > 0.5
        fmax = 0.5;
    end

    % pad the image if it is not square
    [r, c] = size(image);
    if r ~= c
        d = max(r, c);
        imagepadded = ones(d, d, 'uint8')*127;
        cx = fix(d/2 - c/2)+1;
        ry = fix(d/2 - r/2)+1;
        imagepadded(ry:ry+r-1, cx:cx+c-1) = image;
    else
        imagepadded = image;
    end
    
    imdimension = max(size(imagepadded));

    % filter to reduce leakage using Blackman window function
    h = blackman(size(imagepadded, 1), 'symmetric');
    filt = h * h';
    imagewindowed = double(imagepadded) .* filt;

    % Do the FFT
    F = fft2(imagewindowed);

    % Center DC component
    FS = fftshift(F);

    % Get the power spectrum
    pow = log(1 + abs(FS));
    
    % Get the range of frequency band of interest in pixels
    fsmin = max(floor(imdimension * fmin), 1);
    fsmax = min(ceil(imdimension * fmax), size(pow,1));

    % Convert to polar coordinates
    powpol = cart2pol(pow, theta);
    
    % Get the frequency band of interest
    FOI = powpol(fsmin:fsmax, :);
    
    % Get power for all bands in FOI
    FOIPower = sum(FOI, 2); 
     
    % Perform the binning of frequencies. Note that only rad bins can be
    % constructed iff m is 0, otherwise rad-m bins will be made.
    m = mod(length(FOIPower), rad);
    if length(FOIPower) <= rad 
        radialpowerspectrum = FOIPower;
        radialpowerspectrumbinned = FOIPower;
    else
        radialpowerspectrum = FOIPower(1:size(FOIPower, 1) - m, :);
        radialpowerspectrumbinned = sum(reshape(radialpowerspectrum, [], rad), 1);
    end
    quality = max(radialpowerspectrumbinned);
end