function [ nfiqScore ] = computeNfiq(image)
%COMPUTENFIQ Proxy to call the nfiq executable
%   Computes the quality of IMAGE and returns it as score 1-5. 
%
% Syntax: 
%   [ nfiqScore ] = nfiq( image );
%
% Inputs:
%   image      - Path to image on disk or an image as a 2D uint8 matrix
%
% Outputs:
%   nfiqScore  - The score as determined by the NFIQ algorithm. Returns NaN
%   on failure.
%
% Examples:
%   nfiqScore = nfiq('fp_1.bmp');

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

nfiqScore = NaN;
tempImFilename = []; % Contains the path to input image for NFIQ.

% Get the location of nfiq executable
thisdir = fileparts(mfilename('fullpath'));
if ispc
    nfiqexe = 'nfiq.exe';
end
nfiqPath = strcat(thisdir, filesep, nfiqexe);

if exist(nfiqPath, 'file') ~= 2
    fprintf('Could not locate nfiq executable\n')
    return
end

%% Check if input is a path to an image on disk
if isa(image, 'char') ...
&& exist(image, 'file') == 2
    image = fpimread(image);
end

% The image must be stored to disk in 8-bit grayscale RAW for NFIQ
storedImage = false;
if isempty(tempImFilename) ...
&& isa(image, 'uint8')
    %Find a place to store our temporary file
    temporaryDir = pwd;
    fName = strcat('nfiqtmp', dec2hex(randi(2^32, 1)), '.raw');
    tempImFilename = strcat(temporaryDir, filesep, fName); 
    
    [imheight, imwidth] = size(image);
    imdepth = 8;
    % Transpose the image to account for NFIQ expecting row-major order in
    % the raw format
    image = image';
    % Write raw image file
    fid = fopen(tempImFilename, 'w');
    if fid < 0
        warning(['Can''t open file ''' tempImFilename ''' for writing.']);
        fclose(fid);
        return;
    end
    fwrite(fid, image(:), 'uint8'); 
    fclose(fid);
    storedImage = true;
else
    return
end

try
    %% Now call the NFIQ executable to get the score
    nfiqCmd = [nfiqPath, ...
        ' -raw ', ...
        num2str(imwidth),  ',', ...
        num2str(imheight), ',', ...
        num2str(imdepth), ...
        ' ', tempImFilename];
    [~, result] = system(nfiqCmd);
    nfiqScore = str2double(result);
catch ex
    disp(ex)
    if storedImage
        delete(tempImFilename);
    end
end

% Cleanup
if storedImage && (exist(tempImFilename, 'file') == 2)
    delete(tempImFilename);
end

end