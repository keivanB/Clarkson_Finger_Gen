function [Q_RVU, Q_Ratios, Q_RVU_local, RTVTR_Ratios, Q_profile, Q_ridval, Q_change, Q_Blocks, Q_R_Blocks, Q_WLC, Q_RC] = RTVTR_Interface(metric, image, varargin)
% Edited by Keivan Bahmani June 2020
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
   
    
   [Q_RVU, Q_Ratios, Q_RVU_local, RTVTR_Ratios, Q_profile, Q_ridval, Q_change, Q_Blocks, Q_R_Blocks, Q_WLC, Q_RC] = compRTVTR(opts.image, opts.segmentMask, opts.RVU_sBlocksize, opts.RVU_blocksize);
    
    function [Q_RVU, Q_Ratios, Q_RVU_local, rvuRatios, Q_profile, Q_ridval, Q_change, Q_Blocks, Q_R_Blocks, Q_WLC, Q_RC] = compRTVTR(im, maskim, v1sz, blksz)
        [rows, cols] = size(im);
        eblksz = ceil(sqrt(sum(v1sz.^2)));
        blkoffset = ceil((eblksz - blksz)/2); % overlapping border
        mapsize = fix(([rows cols] - (eblksz - blksz))./blksz); 
        maskBseg = false(mapsize);
        % Set a resolution
        scres = 500; 
        blkorient = zeros(mapsize);
        rvuRatios = [];
        Q_Ratios  = {};  % RTVR for each block
        Q_profile = {}; % ridge valley profile
        Q_ridval = {};  % ridge valley profile after regression
        Q_Blocks = {};  % fingerprint block used for calculaition
        Q_R_Blocks = {};
        Q_change = {}; % change in the ridge valley 
        Q_WLC = {};    %  white lines count
        Q_RC = {};     %  ridge count
        Q_RVU_local = zeros(mapsize);
        %BLOCKS = cel(length(blkoffset+1:blksz:rows-(blksz+blkoffset-1)),lenght(blkoffset+1:blksz:cols-(blksz+blkoffset-1))); 
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
                    [rvures, blockRotated, blockCropped, v3, ~, ~, ~, ridval, change] = rvu(blkwim, blkorient(br,bc), v1sz, scres);
                    rvuRatios = [rvuRatios rvures];
                    Q_Ratios{br,bc} = rvures;
                    Q_profile{br,bc} = v3;
                    Q_ridval{br,bc} = ridval;  %ridges 1 valleys 0
                    Q_change{br,bc} = change;
                    Q_Blocks{br,bc} = blockCropped;
                    Q_R_Blocks{br,bc} = blockRotated;
                    Q_RVU_local(br, bc) = std(rvures);   %std of the RTVTR locals, will be used to select 15 highest quality blocks
                    %Q_WLC{br, bc}= length(ridval(ridval == 0));
                    n = 0;    % number of ridges
                    nwl = 0;  % number of white lines (valleys)
                    for i=1:length(ridval)
                        if (i==1)
                           if(ridval(i)==1)
                                %n = n+1; I am not counting the first pixel
                                %as a ridge.
                                b = 1;  % current column for the next loop
                           else
                                b = 0;
                           end
                        else
                           if(ridval(i)==1) % it is a ridge coloumn
                              if(b == 0)    % it is the start of a ridge patern
                                 n = n+1; 
                                 b = 1; 
                              end
                           else
                              if(b == 1)  % it is start of a valley coloumn
                                nwl = nwl+1;
                                b = 0;       
                              else
                                b = 0;  
                              end
                           end
                        end
                    end      
                    Q_RC{br, bc}= n;
                    Q_WLC{br, bc}= nwl;
                end
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
            'ridgevalleythiknessratio' , 'RTVTR', ...
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