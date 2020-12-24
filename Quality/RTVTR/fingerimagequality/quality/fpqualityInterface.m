function [out, qmap] = fpqualityInterface(metric, image, varargin)
% FPQUALITYINTERFACE is the common interface to compute quality for images.
%   The function accepts one or more fingerimage quality features in METRIC
%   and in IMAGE either 1) a full path to an image or 2) an image as a MxN
%   matrix. Quality values for the features listed in METRIC are returned
%   in the output vector OUT with the same length and order as METRIC. For
%   quality features which do compute local quality values, either per
%   pixel or block, the computed quality map is returned in QMAP.
%
% Syntax:
%   [q, map] = fpqualityInterface(metric, image);
%
% Example:
%   metrics  = {'orientationCertaintyLevel', 'gabor'};
%   root     = baseDirectory(); 
%   image    = fullfile(root, 'selftest', 'images', '0002.bmp')
%   [q, map] = fpqualityInterface(metric, image);
%
% Inputs (required):
%   metric  - a single string or a cell of strings with names of metrics to
%             be computed. Valid strings can be retrived by calling
%             fpqualityInterface with no parameters FPQUALITYINTERFACE()
%   image   - A path to a fingerprint image or a MxN matrix of type UINT8.
%
% Input (optional):
% The following parameters may be specified as key value pairs. Default
% values are used if the parameter is left unspecified.
%   segmentType       - kovesi    Use fingerprint segmentation algorithm by 
%                                 Peter Kovesi if none set in segmentMask 
%                                 parameter.
%   segmentMask       - NaN       User supplied segmentatation mask.
%   segmentBlocksize  - 32        Blocksize used for segmentation.
%   segmentThreshold  - 0.1       Threshold for segmentation.
%   qualityBlocksize  - 32        Blocksize for blockwise quality features.
%   OCL_blocksize     - 32        Blocksize for OCL.
%   OCL_sBlocksize    - [32 16]   Size of rotated blocks when aligning 
%                                 according to orientation for OCL.
%   LCS_blocksize     - 32        Blocksize for LCS.
%   LCS_sBlocksize    - [32 16]   Size of rotated blocks when aligning 
%                                 according to orientation for LCS.
%   OFL_blocksize     - 32        Blocksize for OFL.
%   OFL_sBlocksize    - [32 16]   Size of rotated blocks when aligning 
%                                 according to orientation for OFL.
%   OFL_border        - 1         Size of border in pixels for orientation 
%                                 flow quality feature.
%   OFL_angMin        - 0         Tolerance for angular change between 
%                                 neighboring blocks in orientation flow 
%                                 feature.
%   FDA_blocksize     - 32        Blocksize for FDA.
%   FDA_sBlocksize    - [32 16]   Size of rotated blocks when aligning 
%                                 according to orientation for FDA.
%   RVU_blocksize     - 32        Blocksize for RVU.
%   RVU_sBlocksize    - [32 16]   Size of rotated blocks when aligning 
%                                 according to orientation for RVU.
%   GAB_angles        - 8         Number of orientations in filter bank.
%   GAB_freq          - 0.1       GAB filter frequency in cycles per pixel.
%   GAB_sigma         - 6         GAB standard deviation of gaussian.
%   GSH_Tb            - 1         GSH threshold for bad blocks see Shen et al.
%   GSH_Tq            - 2         GSH threshold for quality see Shen et al.
%   GSH_angles        - 8         GSH number of orientations in filter bank.
%   GSH_freq          - 0.12      GSH filter frequency in cycles per pixel.
%   GSH_sigma         - 4         Standard deviation of gaussian.
%   GSH_blockSize     - 30        Blocksize for GSH
% 
% Outputs:
%   out     - A vector of quality values, one for each feature specified in
%             METRIC. The quality values appear in the same order as the
%             features appear in METRIC.
%   qmap    - A cell of matrices, one for each feature specified in METRIC.
%             Each matrix contains the quality map as determined by the
%             METRIC, if supported, otherwise empty.

% If you use this code in a publication please cite the following paper:
% Olsen, M. A.; Smida, V. & Busch, C. Finger image quality assessment features - definitions and evaluation IET Biometrics, Institution of Engineering and Technology, 2015
% The paper can be accessed for free via http://digital-library.theiet.org/content/journals/10.1049/iet-bmt.2014.0055
%
% 2015 Martin Aastrup Olsen, martin.olsen@{cased.de;hig.no}

    if nargin < 2
        out = supportedMetrics();
        warning('At least one valid metric and an image must be passed as input arguments')
        return
    end
    
    % Parse inputs
    opts = parseInput(metric, image, varargin);
    
    % Perform segmentation if no mask is set
    if isnan(opts.segmentMask)
        switch opts.segmentType
            case 'kovesi'
                [ ~, opts.segmentMask, ~] = ridgesegment(opts.image, ...
                                                         opts.segmentBlocksize, ...
                                                         opts.segmentThreshold);
            case 'none'
                opts.segmentMask = true(size(opts.image));
            otherwise
                opts.segmentMask = true(size(opts.image));
        end
    end
    
    out = NaN(numel(opts.metric), 1); % Preallocate output vector
    qmap = cell(numel(opts.metric), 1); % Preallocate output matrix
    for currentMetric = 1:numel(opts.metric)
        try
            switch opts.metric{currentMetric}
                case {'orientationCertaintyLevel', 'OCL'}
                    [out(currentMetric), qmap{currentMetric}] = ...
                        compOcl(opts.image, ...
                                opts.segmentMask, ...
                                opts.OCL_sBlocksize , ...
                                opts.OCL_blocksize);
                case {'localClarityScore', 'LCS'}
                    [out(currentMetric), qmap{currentMetric}] = ...
                        compLcs(opts.image, ...
                                opts.segmentMask, ...
                                opts.LCS_sBlocksize, ...
                                opts.LCS_blocksize);
                case {'orientationFlow', 'OFL'}
                    [out(currentMetric), qmap{currentMetric}] = ...
                        compOf( opts.image, ...
                                opts.segmentMask, ...
                                opts.OFL_sBlocksize, ...
                                opts.OFL_blocksize, ...
                                opts.OFL_angMin, ...
                                opts.OFL_border);
                case {'radialPowerSpectrum', 'RPS'}
                    out(currentMetric) = ...
                        compRPS(opts.image, ...
                                opts.RPS_rad, ...
                                opts.RPS_theta, ...
                                opts.RPS_fmin, ...
                                opts.RPS_fmax);
                case {'frequencyDomainAnalysis', 'FDA'}
                    [out(currentMetric), qmap{currentMetric}] = ...
                        compFda(opts.image, ...
                                opts.segmentMask, ...
                                opts.FDA_sBlocksize, ...
                                opts.FDA_blocksize);     
                case {'ridgeValleyUniformity', 'RVU'}
                    [out(currentMetric), qmap{currentMetric}] = ...
                        compRvu(opts.image, ...
                                opts.segmentMask, ...
                                opts.RVU_sBlocksize, ...
                                opts.RVU_blocksize);
                case {'gabor', 'GAB'}
                    [out(currentMetric), qmap{currentMetric}] = ...
                        compGabor(opts.image, ...
                                  opts.GAB_freq, ...
                                  opts.GAB_sigma, ...
                                  opts.GAB_angles);
                case {'gaborShen', 'GSH'}
                    out(currentMetric) = ...
                        compGaborShen(opts.image, ...
                                      opts.GSH_Tb, ...
                                      opts.GSH_Tq, ...
                                      opts.GSH_angles, ...
                                      opts.GSH_freq, ...
                                      opts.GSH_sigma, ...
                                      opts.GSH_blockSize);
                case {'sigma', 'SIG'}
                    [out(currentMetric)] = ...
                        compSigma(opts.image);
                case {'mu', 'MU'}
                    [out(currentMetric)] = ...
                        compMu( opts.image);
                case {'nfiq', 'NFIQ'}
                    out(currentMetric) = ...
                        computeNfiq(opts.image);
            end
        catch Ex
            out(currentMetric) = NaN;
            disp(Ex.identifier)
            disp(Ex.message)
        end     
    end

    function [Q_GSH] = compGaborShen(image, tb, tq, angles, freq, sigma, blocksize)
       [Q_GSH, GF, GaborAngle, gaborBlkStd, gaborFore, gaborPoor] = ...
           gaborShenQuality(image, tb, tq, angles, freq, sigma, blocksize);
    end
    
    function [Q_GAB, Q_GAB_local] = compGabor(image, freq, sigma, angles)
        [Q_GAB, imGauss, GF, FP_filt, GaborAngle, response, Q_GAB_local] = ...
            gaborQuality(image, freq, sigma, angles);
    end

    function [Q_PRS] = compRPS(image, rad, theta, fmin, fmax)
       [Q_PRS, imagepadded, imagewindowed, FS, pow, powpol, FOI, ...
        FOIPower, radialpowerspectrum, radialpowerspectrumbinned, fsmin, ...
        fsmax] = rps(image, rad, theta, fmin, fmax);
    end
    
    function [Q_MU] = compMu(image)
        Q_MU = mean2(double(image));
    end
 
    function [Q_SIG] = compSigma(image)
        Q_SIG = std2(double(image));
    end

    function [Q_RVU, Q_RVU_local] = compRvu(im, maskim, v1sz, blksz)
        [rows, cols] = size(im);
        eblksz = ceil(sqrt(sum(v1sz.^2)));
        blkoffset = ceil((eblksz - blksz)/2); % overlapping border
        mapsize = fix(([rows cols] - (eblksz - blksz))./blksz); 
        maskBseg = false(mapsize);
        % Set a resolution
        scres = 500; 
        blkorient = zeros(mapsize);
        rvuRatios = [];
        
        Q_RVU_local = zeros(mapsize);
        
        br = 1; bc = 1;
        for r = blkoffset+1:blksz:rows-(blksz+blkoffset-1)
            for c = blkoffset+1:blksz:cols-(blksz+blkoffset-1)
                blkim = im(r:r+blksz-1, c:c+blksz-1);
                maskB1 = maskim(r:r+blksz-1, c:c+blksz-1);
                maskBseg(br,bc) = all(maskB1(:));
                [cova, covb, covc] = covcoef(blkim);

                %% ridge ORIENT local
                blkorient(br,bc) = ridgeorient(cova, covb, covc);

                % overlapping windows (border = blkoffset)
                blkwim = im(r-blkoffset:r+blksz+blkoffset-1, c-blkoffset:c+blksz+blkoffset-1); 
                %% RVU local (only for foreground)
                if(maskBseg(br,bc)==1)
                    rvures = rvu(blkwim, blkorient(br,bc), v1sz, scres);
                    rvuRatios = [ rvuRatios rvures ];
                    Q_RVU_local(br, bc) = std(rvures);
                end;
                %%
                bc = bc+1;
            end
            br = br+1; bc = 1;
        end
        Q_RVU = std(rvuRatios(~isnan(rvuRatios)));
        if all(isnan(rvuRatios(:))) || isempty(rvuRatios)
            Q_RVU = 0;
        end
    end
    
    function [Q_FDA, Q_FDA_local] = compFda(im, maskim, v1sz, blksz)
        [rows, cols] = size(im);
        eblksz = ceil(sqrt(sum(v1sz.^2))); % block size for extraction of slanted block
        blkoffset = ceil((eblksz - blksz)/2); % overlapping border
        mapsize = fix(([rows cols] - (eblksz - blksz))./blksz);  
        maskBseg = false(mapsize);
        % Set a resolution
        scres = 500; 
        blkorient = zeros(mapsize);
        Q_FDA_local  = zeros(mapsize);
        br = 1; bc = 1;
        for r = blkoffset+1:blksz:rows-(blksz+blkoffset-1)
            for c = blkoffset+1:blksz:cols-(blksz+blkoffset-1)
                blkim = im(r:r+blksz-1, c:c+blksz-1);
                maskB1 = maskim(r:r+blksz-1, c:c+blksz-1);
                maskBseg(br,bc) = all(maskB1(:));
                [cova, covb, covc] = covcoef(blkim);

                %% ridge ORIENT local
                blkorient(br,bc) = ridgeorient(cova, covb, covc);

                % overlapping windows (border = blkoffset)
                blkwim = im(r-blkoffset:r+blksz+blkoffset-1, c-blkoffset:c+blksz+blkoffset-1); 

                %% FDA local
                Q_FDA_local(br,bc) = fda(blkwim, blkorient(br,bc), v1sz, scres);
 
                bc = bc+1;
            end
            br = br+1; bc = 1;
        end
        Q_FDA_local(~maskBseg) = NaN; % apply background mask
        Q_FDA =  mean(Q_FDA_local(~isnan(Q_FDA_local))) ;
        if all(isnan(Q_FDA_local(:))) || isempty(Q_FDA_local)
            Q_FDA = 0;
        end
    end
    
    function [Q_OFL, Q_OFL_local] = compOf(im, maskim, v1sz, blksz, angmin, border)
        [rows, cols] = size(im);
        eblksz = ceil(sqrt(sum(v1sz.^2)));
        blkoffset = ceil((eblksz - blksz)/2); 
        mapsize = fix(([rows cols] - (eblksz - blksz))./blksz); 
        maskBseg = false(mapsize);
        
        blkorient = zeros(mapsize);
        br = 1; bc = 1;
        for r = blkoffset+1:blksz:rows-(blksz+blkoffset-1)
            for c = blkoffset+1:blksz:cols-(blksz+blkoffset-1)
                blkim = im(r:r+blksz-1, c:c+blksz-1);
                maskB1 = maskim(r:r+blksz-1, c:c+blksz-1);
                maskBseg(br,bc) = all(maskB1(:));
                [cova, covb, covc] = covcoef(blkim);

                %% ridge ORIENT local
                blkorient(br,bc) = ridgeorient(cova, covb, covc);
                
                bc = bc+1;
            end
            br = br+1; bc = 1;
        end

        % get the diff of orient. angles from neighbouring blocks
        loqall = blockproc(blkorient, [1 1], ...
                           @orientangdiff, ...
                           'BorderSize', [border border], ...
                           'TrimBorder', false);

        angdiff     = deg2rad(90-angmin);
        angmin      = deg2rad(angmin);
        Q_OFL_local        = zeros(size(loqall));

        % overlapping window: if one of the surrouding blocks from which 
        % the anglediff was computed background, exclude window from comp.
        maskBloqseg = logical(blockproc(maskBseg, [1 1], ...
                                        @(x) all(x.data(:)), ...
                                        'BorderSize', [border border], ...
                                        'TrimBorder', false)); 
        
        % (angle) mask - only blocks with angle change > angmin
        maskBang = loqall > angmin; 
        % (local orientation quality) mask of FOREGROUND blocks > angmin deg
        maskBloq = maskBang & maskBloqseg;
        % map of local orientation quality scores
        Q_OFL_local(maskBloq) = (loqall(maskBloq) - angmin) ./ angdiff;

        % global orientation flow quality score (GOQS) only from FOREGROUND (all neighbouring)
        % blocks (background = 0) 
        Q_OFL = 1 - mean(Q_OFL_local(maskBloqseg));
        if ~any(maskBloqseg(:))
            Q_OFL = 0;
        end
    end
    
    function [Q_LCS, Q_LCS_local] = compLcs(im, maskim, v1sz, blksz)
        [rows, cols] = size(im);
        eblksz = ceil(sqrt(sum(v1sz.^2))); % block size for extraction of slanted block
        blkoffset = ceil((eblksz - blksz)/2); % overlapping border
        mapsize = fix(([rows cols] - (eblksz - blksz))./blksz);
        maskBseg = false(mapsize);
        
        % Set a resolution
        scres = 500; 
        blkorient = zeros(mapsize);
        Q_LCS_local = zeros(mapsize);
        br = 1; bc = 1; 
        for r = blkoffset+1:blksz:rows-(blksz+blkoffset-1)
            for c = blkoffset+1:blksz:cols-(blksz+blkoffset-1)
                blkim = im(r:r+blksz-1, c:c+blksz-1);
                maskB1 = maskim(r:r+blksz-1, c:c+blksz-1);
                maskBseg(br,bc) = all(maskB1(:));
                [cova, covb, covc] = covcoef(blkim);
                %% ridge ORIENT local
                blkorient(br,bc) = ridgeorient(cova, covb, covc);
                %% LCS local
                % overlapping windows (border = blkoffset)
                blkwim = im(r-blkoffset:r+blksz+blkoffset-1, c-blkoffset:c+blksz+blkoffset-1); 
                Q_LCS_local(br,bc) = loclar(blkwim, blkorient(br,bc), v1sz, scres);
                bc = bc+1;
            end
            br = br+1; bc = 1;
        end
        Q_LCS_local(~maskBseg) = NaN;
        Q_LCS = mean(Q_LCS_local(~isnan(Q_LCS_local)));
        if all(isnan(Q_LCS_local(:))) || isempty(Q_LCS_local)
            Q_LCS = 0;
        end
    end

function [Q_OCL, Q_OCL_local] = compOcl(im, maskim, v1sz, blksz)
    [rows, cols] = size(im);
    eblksz = ceil(sqrt(sum(v1sz.^2))); 
    blkoffset = ceil((eblksz - blksz)/2); 
    mapsize = fix(([rows cols] - (eblksz - blksz))./blksz);
    maskBseg = false(mapsize);

    Q_OCL_local = zeros(mapsize);
    br = 1; bc = 1; 
    for r = blkoffset+1:blksz:rows-(blksz+blkoffset-1)
        for c = blkoffset+1:blksz:cols-(blksz+blkoffset-1)
            blkim = im(r:r+blksz-1, c:c+blksz-1);
            maskB1 = maskim(r:r+blksz-1, c:c+blksz-1);
            maskBseg(br,bc) = all(maskB1(:));
            [cova, covb, covc] = covcoef(blkim);
            Q_OCL_local(br,bc) = ocl(cova, covb, covc);
            bc = bc+1;
        end
        br = br+1; bc = 1;
    end
    Q_OCL_local(~maskBseg) = NaN; % mask the background blocks
    Q_OCL = mean(Q_OCL_local(~isnan(Q_OCL_local))); 
    if all(isnan(Q_OCL_local(:))) || isempty(Q_OCL_local)
        Q_OCL = 0;
    end
end
    
%% INPUT HANDLING
    % parse the inputs
    function [opts] = parseInput(metric, image, args)
        if ischar(metric)
            metric = {metric};
        end
        if ischar(image)
            image = fpimread(image);
        end
        p = inputParser;
        
        p.addRequired('metric', @(x) areMetricsValid(x));
        p.addRequired('image', @(x) isa(x, 'uint8') && ismatrix(x));
        p.addParamValue('segmentType', 'kovesi');
        p.addParamValue('segmentMask', NaN, @(x) islogical(x) && ismatrix(x));
        p.addParamValue('segmentBlocksize', 32);
        p.addParamValue('segmentThreshold', 0.1);
        p.addParamValue('OCL_blocksize', 32);
        p.addParamValue('OCL_sBlocksize', [32 16]);
        p.addParamValue('LCS_blocksize', 32);
        p.addParamValue('LCS_sBlocksize', [32 16]);
        p.addParamValue('OFL_blocksize', 32);        
        p.addParamValue('OFL_sBlocksize', [32 16]);
        p.addParamValue('OFL_border', 1);
        p.addParamValue('OFL_angMin', 0);
        p.addParamValue('FDA_blocksize', 32);        
        p.addParamValue('FDA_sBlocksize', [32 16]);
        p.addParamValue('RVU_blocksize', 32);        
        p.addParamValue('RVU_sBlocksize', [32 16]);
        p.addParamValue('RPS_rad', 10);
        p.addParamValue('RPS_theta', 100);
        p.addParamValue('RPS_fmin', 0.06);
        p.addParamValue('RPS_fmax', 0.18)
        p.addParamValue('GAB_angles', 8); % GAB constants according to 
        p.addParamValue('GAB_freq', 0.1); % paper by M. A. Olsen et al.
        p.addParamValue('GAB_sigma', 6);  % see gaborQuality.m
        p.addParamValue('GSH_Tb', 1);  % GSH constants according to paper 
        p.addParamValue('GSH_Tq', 2);  % by L. Shen et al. 
        p.addParamValue('GSH_angles', 8); % see gaborShenQuality.m
        p.addParamValue('GSH_freq', 0.12);
        p.addParamValue('GSH_sigma', 4);
        p.addParamValue('GSH_blockSize', 30);
        p.parse(metric, image, args{:});
        opts = p.Results();
    end
    
    function [validMetrics] = supportedMetrics()
        validMetrics = { ...
            'orientationCertaintyLevel', 'OCL', ...
            'localClarityScore'        , 'LCS', ...
            'orientationFlow'          , 'OFL', ...
            'radialPowerSpectrum'      , 'RPS', ...
            'frequencyDomainAnalysis'  , 'FDA', ...
            'ridgeValleyUniformity'    , 'RVU', ...
            'gabor'                    , 'GAB', ...
            'gaborShen'                , 'GSH', ...
            'sigma'                    , 'SIG', ...
            'mu'                       , 'MU', ...
            'nfiq'                     , 'NFIQ', ...
            };
    end

    % Check if metrics is part of the valid metrics list
    function [metricsValid] = areMetricsValid(metrics)
        validMetrics = supportedMetrics;
        for ind = 1:numel(metrics)
            metricsValid = true;
            if ~any(strcmpi(metrics(ind), validMetrics))
                metricsValid = false;
                warning('metric not valid!');
                break;
            end
        end
    end
end