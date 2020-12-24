function impol = cart2pol(imcart, thetan)
%IMPOL  Interpolation from cart to polar coordinates, function used by f2radpow
% Syntax:
%    impol = cart2pol(imcart, thetan)
% Inputs:
%       imcart      image in cartesian coordinates
%       thetan      theta coarness (# of deltaTheta)
% Outputs:
%       impol       image in polar coordinate system

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2012 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}
% 2011 Master Thesis, Vladimir Smida, vladimir.smida@[cased.de|gmail.com]
% FIT VUT, Czech Republic & CASED, Germany

    % square image: width = height
    imsize = size(imcart,1);
    imcen = imsize/2 + 0.5;

    % sampling angle step
    thetadelta = pi/(thetan-1);
    % matrices of all thetas, all rads
    [theta,rad] = meshgrid(0:thetadelta:pi, 0:imcen-1);

    % relation catesian-polar coordinates
    % reconstruct cartesian indicies from polar
    imx = rad .* cos(theta) + imcen;
    imy = rad .* sin(theta) + imcen;

    % interpolate form cart to reconstructed catesian -> polar coordinates
    [x,y] = meshgrid(1:imsize);
    impol = interp2(x, y, double(imcart), imx, imy);
end