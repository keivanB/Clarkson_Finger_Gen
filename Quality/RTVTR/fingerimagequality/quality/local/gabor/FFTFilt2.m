function [ImageOut] = FFTFilt2(Filt,ImageIn)
% Filter image by a defined filter.
%
% Syntax:
%   [ImageOut] = FFTFilt2(Filt,ImageIn),
%
% Inputs:
%   Filt        - defined filter
%   ImageIn     - input image
%
% Outputs:
%   ImageOut    - output image
%
% Updates:
%   2011-08-15          Initial version
%   2012-04-26          Clean up the code

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% Haiyun Xu (haiyun.xu@gmail.com)
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

    [Mi,Ni] = size(ImageIn);
    [Mf,Nf] = size(Filt);

    F = zeros(Mi,Ni);
    F(round(Mi/2 - Mf/2) + (1:Mf), round(Ni/2 - Nf/2) + (1:Nf)) = conj(Filt);

    if isequal(F, real(F)) && isequal(ImageIn, real(ImageIn)),
        ImageOut = fftshift(real(ifft2(fft2(ImageIn) .* conj(fft2(F, Mi, Ni)))));
    else
        ImageOut = fftshift(ifft2(fft2(ImageIn) .* conj(fft2(F, Mi, Ni))));
    end
end