function crf = camresponse(filenames, varargin)
%CAMRESPONSE Estimate camera response function curve
%   CRF = CAMRESPONSE(FILES) estimates the camera response function (CRF)
%   curve from a set of spatially registered low dynamic range images
%   listed in FILES. FILES is an array of strings or a cell array of
%   character vectors. These files must contain EXIF exposure metadata.
%
%   CRF = CAMRESPONSE(IMDS) estimates the CRF curve from the set of
%   spatially registered low dynamic range images in imageDatastore IMDS.
%
%   CRF = CAMRESPONSE(___, 'ExposureTimes', EXPTIMES) uses the values in
%   the EXPTIMES vector as the exposure times for the low dynamic range
%   images during computation. This parameter overrides EXIF exposure
%   metadata. The number of elements in EXPTIMES must be the same as the
%   number of low dynamic range images.
%
%   Class Support 
%   ------------- 
%   The input images can be color or grayscale. They can have any bit
%   depth. Preferred bit depth for low dynamic range image is 8 or 16.
%   The output CRF is a double-precision.
%
%   Notes
%   -----
%   [1] This function requires a minimum of two images with different
%   exposure times. A larger number of images yields a better estimate of
%   CRF at the expense of more processing time
%
%   References
%   ----------
%   [1] Paul E. Debevec and Jitendra Malik. "Recovering high dynamic range
%   radiance maps from photographs." ACM SIGGRAPH 2008 classes. ACM, 2008.
%
%   Example 1:
%   ----------
%   %  Estimate CRF from a series of six low dynamic range images that 
%   %  share the same F-Stop number and have different exposure times.
%
%   files = ["office_1.jpg", "office_2.jpg", "office_3.jpg", ...
%               "office_4.jpg", "office_5.jpg", "office_6.jpg"];
%
%   crf = camresponse(files);
%
%   figure; 
%   hold on
%   plot(crf(:,1), '-r', 'LineWidth', 2);
%   plot(crf(:,2), '-g', 'LineWidth', 2);
%   plot(crf(:,3), '-b', 'LineWidth', 2);
%   xlabel('Image intensity'), ylabel('log exposure, irradiance');
%   title('Camera Response Function');
%   axis square; grid on;
%   hold off
%
%   Example 2:
%   ----------
%   %  Estimate CRF curve from a series of six low dynamic range images in 
%   %  imageDatastore with given exposure time parameters.
%
%   setDir = fullfile(toolboxdir('images'), 'imdata','office_*');
%   imds = imageDatastore(setDir);
%   expTimes = [0.0333, 0.1000, 0.3333, 0.6250, 1.3000, 4.0000];
%   crf = camresponse(imds,'ExposureTimes',expTimes);
%
%   See also MAKEHDR, HDRREAD, HDRWRITE.

%   Copyright 2018 The MathWorks, Inc.

% Parse and check inputs.
narginchk(1,3)
validateattributes(filenames,{'cell','string','matlab.io.datastore.ImageDatastore'},...
    {'nonempty'},mfilename,'files',1);
if isa(filenames,'matlab.io.datastore.ImageDatastore')
    validateattributes(filenames.Files, {'cell'}, ...
        {'nonempty'}, mfilename, 'files',1);
end
if iscell(filenames)
    numImage = numel(filenames);
    fname = filenames{1};
elseif isstring(filenames)
    filenames = cellstr(filenames);
    numImage = numel(filenames);
    fname = filenames{1};
elseif isa(filenames,'matlab.io.datastore.ImageDatastore')
    filenames = filenames.copy();
    filenames.reset();
    numImage = numel(filenames.Files);
    fname = filenames.Files{1};
end

% Get the metadata from the first image file
meta = getMetaData(fname);

% Number of channels
if isfield(meta, 'BitsPerSample')
    numChannels = numel(meta.BitsPerSample);
elseif strcmp(meta.ColorType, 'truecolor')
    numChannels = 3;
else
    numChannels = 1;
end

varargin = matlab.images.internal.stringToChar(varargin);
options = parseArgs(varargin{:});
validateOptions(meta, options, numImage);

if ~isempty(options.ExposureTimes)
    exposureTime = options.ExposureTimes;
else
    exposureTime = zeros(1,numImage);
    for p = 1:numImage
        if isa(filenames,'matlab.io.datastore.ImageDatastore')
            fname = filenames.Files{p};
        else
            fname = filenames{p};
        end
        exposureTime(p) = getExposure(fname);
    end
end

% Loading logarithm exposure times
logExpo = log(exposureTime);

% Weighting Function
bitsPerSample = meta.BitDepth / numChannels;
maxVal = 2^bitsPerSample-1;
weight = 0:maxVal;
weight(ceil(maxVal/2)+1:maxVal+1) = maxVal - weight(ceil(maxVal/2)+1:maxVal+1);

% Loading the images
image = zeros(meta.Height,meta.Width,numChannels,numImage);
imageStack = loadImage(filenames, meta, numImage, numChannels);
for itr = 1:numImage
    image(:,:,:,itr) = imageStack{itr};
end

% regularization for CRF estimation
lambda = 40;

% Total number of pixels per image per channel
nPixChannel = meta.Height*meta.Width;
% Storing each channel pixels separately
nAllPix = reshape(image, [nPixChannel, numChannels, numImage]);

% Selecting limited pixels for CRF estimation
samplePix = floor(nPixChannel*0.002); % sample rate
pix = 1:samplePix:nPixChannel;

nSelPix = numel(pix); % number of pixels selected

nlevel = 2^bitsPerSample;

% system Ax = b
if bitsPerSample == 8
    A = zeros((nSelPix*numImage) + nlevel + 1, nlevel + nSelPix);    
else
    A = sparse((nSelPix*numImage) + nlevel + 1, nlevel + nSelPix);
end
b = zeros(size(A,1),1);
crf = zeros(nlevel,numChannels);

for itrChannel = 1:numChannels
    imgPix = reshape(nAllPix(pix, itrChannel, :), [nSelPix, numImage]);
    crf(:,itrChannel) = crfPerChannel(imgPix,A,b,nlevel,weight,logExpo,lambda);
end


%--------------------------------------------------------------------------
function crfCurve = crfPerChannel(pix,A,b,nlevel,weight,logExpo,lambda)
% Data-fitting equations [w(Zij)*(g(Zij) - ln(Ei)) = w(Zij)*ln(deltatj)]
count = 1;
for itr=1:size(pix,1)
    for n=1:size(pix,2)
        index = pix(itr,n)+1;
        A(count,index) = weight(index);
        A(count,nlevel+itr) = -weight(index);
        b(count) = weight(index) * logExpo(n);
        count = count+1;
    end
end

% Constraint [g(Zmid) = 0]
A(count,nlevel/2+1) = 1;
count=count+1;

% Smoothness equations [g(z-1) - 2g(z) + g(z+1) = 0]
for itr=1:nlevel-2
    A(count,itr) = lambda*weight(itr+1);
    A(count,itr+1) = -2*lambda*weight(itr+1);
    A(count,itr+2) = lambda*weight(itr+1);
    count = count+1;
end

% Solve the system Ax = b
% Disable the warnings about conditioning for singular and
% nearly singular matrices
warnState = warning('off', 'MATLAB:rankDeficientMatrix');
x = A\b;
% Restore the warning states to their original settings
warning(warnState) 

crfCurve = x(1:nlevel);


%--------------------------------------------------------------------------
function meta = getMetaData(filename)

try
    meta = imfinfo(filename);
catch ME
    if (isequal(ME.identifier, 'MATLAB:imagesci:imfinfo:fileOpen'))
        error(message('images:makehdr:fileNotFound', filename));
    else
        % Unexpected error
        rethrow(ME)
    end
end

% If there are several images in the file,
% use the meta data of the first one
if ~isscalar(meta)
    meta = meta(1);
end


%--------------------------------------------------------------------------
function exif = getExposureDataFromFile(filename)
% Extract exposure metadata from a file containing EXIF.

meta = getMetaData(filename);
if isfield(meta, 'DigitalCamera')
    exif = meta.DigitalCamera;
else
    error(message('images:camresponse:exifFormat', filename));
end

if (isempty(exif) || ...
        ~isstruct(exif) || ...
        ~isfield(exif, 'FNumber') || ...
        ~isfield(exif, 'ExposureTime'))
    
    error(message('images:camresponse:noExposureMetadata', filename))
end


%--------------------------------------------------------------------------
function [baseTime, baseFStop] = getExposure(filename)
% Extract the exposure values from a file containing EXIF metadata.

exif = getExposureDataFromFile(filename);
baseFStop = exif.FNumber;
baseTime = exif.ExposureTime;


%--------------------------------------------------------------------------
function options = parseArgs(varargin)
% Parse the parameter-value pairs, getting default values.

parser = inputParser();
parser.FunctionName = mfilename;

% NameValue 'ExposureTimes'
defaultExposureTimes = [];
validateExposureTimes = @(x) validateattributes(x, ...
    {'single','double'}, ...
    {'nonempty', 'vector', 'real','nonnan', 'finite','positive'}, ...
    mfilename,'ExposureTimes');
parser.addParameter('ExposureTimes', ...
    defaultExposureTimes, ...
    validateExposureTimes);

parser.parse(varargin{:});
options = parser.Results;


%--------------------------------------------------------------------------
function validateOptions(meta, options, numImage)
% make sure images are 'truecolor' or 'grayscale' only
if ~(strcmp(meta.ColorType, 'truecolor') ...
        || strcmp(meta.ColorType, 'grayscale'))
    error(message('images:validate:invalidImageFormat',meta.Filename))
end

% Make sure that at least two images are provide for CRF estimation
if (numImage < 2)
    error(message('images:camresponse:tooLessFiles')) 
end

% Make sure that the correct number of exposure-related values are given.
if (~isempty(options.ExposureTimes) ...
        && (numel(options.ExposureTimes) ~= numImage))
    error(message('images:camresponse:wrongExposureCount')) 
end


%--------------------------------------------------------------------------
function images = loadImage(files, meta, numImage, numPlanes)
% Read all LDR images
images = cell(1,numImage);

% Load all images
if isa(files,'matlab.io.datastore.ImageDatastore')
    for imgCount = 1:numImage
        images{imgCount} = readimage(files,imgCount);
    end
else
    for imgCount = 1:numImage
        images{imgCount} = imread(files{imgCount});
    end
end

if isa(files,'matlab.io.datastore.ImageDatastore')
    files = files.Files'; % To validate dimensions and bitdepth
end

% Make sure all images have same dimensions
for imgCount = 1:numImage
    if ~isequal(size(images{imgCount},1),meta.Height) ...
            || ~isequal(size(images{imgCount},2),meta.Width) ...
            || ~isequal(size(images{imgCount},3),numPlanes)
        error(message('images:makehdr:imageDimensions', files{imgCount}));
    end
end

% Make sure all images have same bitdepth
for imgCount = 1:numImage
    metatemp = getMetaData(files{imgCount});
    if ~isequal(metatemp.BitDepth,meta.BitDepth)
        error(message('images:makehdr:imageBitDepth', files{imgCount}));
    end
end

