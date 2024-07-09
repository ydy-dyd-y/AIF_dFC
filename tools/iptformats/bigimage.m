classdef bigimage < handle
    %bigimage Out-of-core processing of very large images
    %   The bigimage class is for processing large images that are too big
    %   to fit in memory or when processing the images might use up all the
    %   free memory. These large images are represented as blocks of data
    %   that can be independently loaded and processed. They also often
    %   contain more than one level, which show the same image at different
    %   magnification values. This class leverages these storage formats to
    %   make the processing more efficient.
    %
    %   IMG = bigimage(FILENAME) creates a bigimage object that wraps the
    %   large image contained in FILENAME. IMG is a 1-by-1 object with
    %   various properties describing the image and methods to manipulate
    %   the data it contains. IMG can correspond to either a
    %   multi-resolution image or an image that holds only one resolution
    %   level. The supported formats of FILENAME are TIFF and BigTIFF.
    %
    %   IMG = bigimage(DIRNAME) creates a bigimage object that wraps the
    %   bigimage directory created by using the OutputDirectory parameter in
    %   the apply method, or by the write method.
    %
    %   IMG = bigimage(VARNAME) creates a bigimage object from a variable
    %   in the workspace. VARNAME is either MxN or a MxNxChannels numeric
    %   matrix.
    %
    %   IMG = bigimage(SPATIALREF, NUMCHANNELS, DATATYPE) creates a
    %   writable bigimage object with NUMCHANNELS channels and pixels of
    %   type DATATYPE. The number of resolution levels and their sizes are
    %   based on the N element imref2d array SPATIALREF. The image data for
    %   IMG is uninitialized at first but can be set with the setBlock
    %   function.
    %
    %   IMG = bigimage(LEVELSIZES, NUMCHANNELS, DATATYPE) creates a
    %   writable bigimage object with NUMCHANNELS channels and pixels of
    %   type DATATYPE. The number and size of the resolution levels are
    %   given by the N-by-2 array of integers LEVELSIZES, where N is the
    %   number of resolution levels. The first column of LEVELSIZES is the
    %   number of rows for each level, and the second column is the number
    %   of columns. The image data for IMG is uninitialized at first but
    %   can be set with the setBlock function.
    %
    %   IMG = bigimage(___, Name, Value) can be used to specify additional
    %   parameters. 
    %   Parameters include:
    %
    %    BlockSize          - The block size. A 1-by-2 or a L-by-2 array of
    %                         blocksizes. Where L is equal to number of
    %                         levels. The block size is the smallest unit
    %                         of data that can be read or written. When
    %                         specified as a constructor parameter, it
    %                         should be 1-by-2. Default value is
    %                         automatically chosen based on data source.
    %
    %    UnloadedValue      - Default pixel value for unloaded blocks. A
    %                         scalar or 1-by-1-by-Channels value of
    %                         ClassUnderlying type that is used to fill a
    %                         block when it does not exists in the
    %                         Datasource. This can happen when a bigimage
    %                         is created by one of the write constructors
    %                         or for blocks of a bigimage that did not
    %                         complete processing due to an aborted apply
    %                         invocation. A scalar value will be extended
    %                         to 1-by-1-by-Channels if Channels>1. The
    %                         default value is 0.
    %
    %    SpatialReferencing - imref2d object specifying pixel location and
    %                         extents. An array of imref2d objects
    %                         describing the location and pixel spacing of
    %                         each resolution level.
    %
    %
    %   bigimage properties:
    %       SpatialReferencing    - imref2d object specifying pixel location and extents.
    %       UnloadedValue         - Default pixel value for unloaded blocks
    %
    %   bigimage properties: (ReadOnly)
    %       BlockSize               - Block size
    %       Channels                - Number of channels
    %       ClassUnderlying         - Pixel datatype
    %       CoarsestResolutionLevel - Lowest resolution image level
    %       DataSource              - Location for data backing the image
    %       SourceDetails           - Source metadata
    %       FinestResolutionLevel   - Highest resolution image level
    %       LevelSizes              - Image size at each level
    %
    %   bigimage methods:
    %       apply        - Process image
    %       isequal      - Compare two bigimage objects for equality
    %       getRegion    - Get arbitrary portion of image
    %       getBlock     - Get specific block of image
    %       getFullLevel - Get full image at specified level
    %       setBlock     - Set data into specific block of image
    %       write        - Write image contents to new file
    %
    %   Example 1
    %   ---------
    %   % Create and visualize a bigimage
    %     bim = bigimage('tumor_091R.tif');
    %     disp(bim);
    %     bh = bigimageshow(bim);
    %     % Zoom in to see finer details, check the resolution being used
    %     % at any time by inspecting this property:
    %     disp(bh.ResolutionLevel);  
    %
    %   Example 2
    %   ---------
    %   % Create a bigimage from a variable
    %     im = imread('concordaerial.png');
    %     bim = bigimage(im);
    %     disp(bim);
    %     % Default spatial referencing sets up the image origin at (0.5,
    %     % 0.5), with pixels having unit length and width.
    %     disp(bim.SpatialReferencing)
    %
    %   Example 3
    %   ---------
    %   % Create an empty, writable bigimage.
    %     bim = bigimage([1000 1500], 3, 'uint8',...
    %                      'UnloadedValue', uint8([0 255 0]),...
    %                      'BlockSize',[100 100]);
    %     % Write a 'red' block 
    %     rblock = zeros(100,100,3,'uint8');    
    %     rblock(:,:,1) = 255;
    %     % Note that the default coordinate system starts at [0.5, 0.5]
    %     setBlock(bim, 1, [400.5,400.5], rblock);
    %    
    %     bigimageshow(bim);     
    %
    %   Example 4
    %   ---------
    %   % Setting custom spatial referencing information
    %     bim = bigimage('tumor_091R.tif');        
    %     % Distance between pixel centers at finest level. Obtained from
    %     % raw data available at
    %     % https://camelyon17.grand-challenge.org/Data/. 
    %     pext = 0.000226316; % (in milli-meters)
    %
    %     % Assume the top left edge of the first pixel starts a (0,0),
    %     % compute the bottom right edge of the last pixel and set the
    %     % spatial referencing:   
    %     imageSize = bim.SpatialReferencing(1).ImageSize;
    %     bim.SpatialReferencing(1).XWorldLimits = [0 pext*imageSize(2)];
    %     bim.SpatialReferencing(1).YWorldLimits = [0 pext*imageSize(1)];
    %
    %     % Since the other levels cover the same world extents, set their
    %     % spatial extents to match    
    %     bim.SpatialReferencing(2).XWorldLimits = [0 pext*imageSize(2)];
    %     bim.SpatialReferencing(2).YWorldLimits = [0 pext*imageSize(1)];
    %     bim.SpatialReferencing(3).XWorldLimits = [0 pext*imageSize(2)];
    %     bim.SpatialReferencing(3).YWorldLimits = [0 pext*imageSize(1)];
    %
    %     % View image with updated coordinates
    %     h = bigimageshow(bim);
    %
    %   See also bigimageDatastore, bigimageshow.
    
    % Copyright 2018-2019 The MathWorks, Inc.
    
    properties (Access = public)
        
        %SpatialReferencing - imref2d object specifying pixel location and extents
        %  An array of imref2d objects describing the location and pixel
        %  spacing of each resolution level.
        SpatialReferencing
        
        %UnloadedValue - Default pixel value for unloaded blocks
        %  A scalar or 1-by-1-by-Channels value of ClassUnderlying type
        %  that is used to fill a block when it does not exists in the
        %  Datasource. This can happen when a bigimage is created by one of
        %  the write constructors or for blocks of a bigimage that did not
        %  complete processing due to an aborted apply invocation. A scalar
        %  value will be extended to 1-by-1-by-Channels if Channels>1. The
        %  default value is 0.
        UnloadedValue
    end
    
    
    properties(Dependent, SetAccess=private)
        %LevelSizes - Image size at each level
        % A L-by-2 array of integers, where L is the number of resolution
        % levels. The first column of LevelSizes is the number of rows for
        % each level, and the second column is the number of columns.
        LevelSizes
        
        %SourceDetails - Source metadata
        %  Struct containing metadata from the source.
        SourceDetails
        
        %Channels - Number of channels
        %  The number of channels or hyperspectral bands in the
        %  image.
        Channels
               
        %CoarsestResolutionLevel - Lowest resolution image level
        %  The level of the image where each pixel covers the largest area
        %  in world coordinates. That is, the level with the lowest or
        %  coarsest resolution data, which corresponds to the smallest size
        %  image.
        CoarsestResolutionLevel
        
        %FinestResolutionLevel - Highest resolution image level
        %  The level of the image where each pixel covers the smallest area
        %  in world coordinates. That is, the level with the highest or
        %  finest resolution data, which corresponds to the largest size
        %  image.
        FinestResolutionLevel
        
        %DataSource - Location for data backing the image
        %  String containing the location from where data is read. If the
        %  bigimage object was not constructed using a filename, this
        %  property will be empty.
        DataSource
    end
    
    properties(Dependent, Hidden)
        NumLevels
    end
    
    
    properties (Hidden)
        %Adapter - Adapter object used for all file I/O
        %  A subclass of images.internal.adapters.BigImageAdapter
        %
        Adapter
    end
    
    properties (Access = private)
        ClassNames
        PixelLabelIDs
        
        WaitBarVars
        
        CleanUpAdapterDir = false
    end
    
    properties (Transient, Access = private, NonCopyable)
        Futures
        DoneDataQueue
    end
    
    properties (GetAccess = public, SetAccess = private)
        % BlockSize - The block size
        %  A 1-by-2 or a L-by-2 array of blocksizes. Where L is equal to
        %  number of levels. The block size is the smallest unit of data
        %  that can be read or written. When specified as a constructor
        %  parameter, it should be 1-by-2.
        %  Default value is automatically chosen based on data source.
        BlockSize
        
        %ClassUnderlying - Pixel datatype
        %  Datatype of the pixels in the image. This is a MATLAB datatype,
        %  such as "uint8" or "single".
        ClassUnderlying       
    end
    
    methods (Access = public)
        function obj = bigimage(firstArg, varargin)
            %   IMG = bigimage(FILENAME) creates a bigimage object that
            %   wraps the large image contained in FILENAME. IMG is a
            %   1-by-1 object with various properties describing the image
            %   and methods to manipulate the data it contains. IMG can
            %   correspond to either a multi-resolution image or an image
            %   that holds only one resolution level. The supported formats
            %   of FILENAME are TIFF and BigTIFF
            %
            %   IMG = bigimage(DIRNAME) creates a bigimage object that
            %   wraps the big image directory created by using the
            %   OutputDirectory parameter in the apply method, or by the
            %   write method.
            %
            %   IMG = bigimage(VARNAME) creates a bigimage object from a
            %   variable in the workspace. VARNAME is either MxN or a
            %   MxNxChannels numeric matrix.
            %
            %   IMG = bigimage(SPATIALREF, NUMCHANNELS, DATATYPE) creates a
            %   writable bigimage object with NUMCHANNELS channels and
            %   pixels of type DATATYPE. The number of resolution levels
            %   and their sizes are based on the N element imref2d array
            %   SPATIALREF. The image data for IMG is uninitialized at
            %   first but can be set with the setBlock function.
            %
            %   IMG = bigimage(LEVELSIZES, NUMCHANNELS, DATATYPE) creates a
            %   writable bigimage object with NUMCHANNELS channels and
            %   pixels of type DATATYPE. The number and size of the
            %   resolution levels are given by the N-by-2 array of integers
            %   LEVELSIZES, where N is the number of resolution levels. The
            %   first column of LEVELSIZES is the number of rows for each
            %   level, and the second column is the number of columns.
            
            if (nargin == 1) && isstruct(firstArg)
                % Called by loadobj() --> empty initialization.
                return
            end
            
            % Datatype is a string, parser will think its a PV, so
            % convert it as one.
            if nargin>2 && isnumeric(varargin{1})
                % Three or more inputs, and second is numeric
                varargin{end+1} = 'Datatype';
                varargin{end+1} = varargin{2};
                varargin(2) = [];
            end
            
            parser = inputParser;
            parser.FunctionName = mfilename;
            parser.CaseSensitive = false;
            parser.PartialMatching = true;
            parser.KeepUnmatched = false;
            
            parser.addOptional('NumChannels', []);
            parser.addParameter('DataType', 'double');
            parser.addParameter('UnloadedValue', 0);
            parser.addParameter('SpatialReferencing', []);
            %parser.addParameter('ClassNames', string([]), @validateClassNames)
            %parser.addParameter('PixelLabelIDs', {},...
            %    @(p)validateattributes(p, {'numeric', 'cell'},...
            %    {'nonempty'}, mfilename, 'PixelLabelIDs'));
            parser.addParameter('BlockSize', [1024 1024]);
            
            
            parser.parse(varargin{:});
            results = parser.Results;
            
            % Validate
            if ~contains('NumChannels', parser.UsingDefaults)
                validateattributes(results.NumChannels, {'numeric'}, ...
                    {'scalar', 'positive', 'integer'},mfilename, 'NUMCHANNELS');
            end
            if ~contains('DataType', parser.UsingDefaults)
                results.DataType = validatestring(results.DataType,...
                    ["uint8", "int8", "uint16", "int16", "uint32",...
                    "int32", "double", "single", "logical"],...
                    mfilename, 'DATATYPE');
            end
            validateattributes(results.BlockSize, {'numeric'}, ...
                {'positive', 'numel', 2, 'row','integer'}, ...
                mfilename, 'BlockSize');
            
            % If NumChannels is specified, datatype MUST be specified
            isWritable = false;
            if ~contains('NumChannels', parser.UsingDefaults) % channels is specified
                if contains('DataType', parser.UsingDefaults) % datatype is not
                    error(message('images:bigimage:missingDatatype'));
                else
                    % Both channels and datatype is specified
                    isWritable = true;
                end
            end
            
                
            if isWritable
                obj.BlockSize = double(results.BlockSize);
                
                % firstArg is either spatial ref or levelsizes
                if isa(firstArg, 'imref2d')
                    spatialRef = firstArg;
                    validateattributes(spatialRef, {'imref2d'}, ...
                        {'nonempty', 'vector'}, 'bigimage', 'SpatialReferencing')                    
                else
                    validateattributes(firstArg, "numeric", ...
                        {"nonempty", "2d", "integer", "positive", "ncols", 2},...
                        mfilename, 'LEVELSIZES');
                    resolutionLevelSizes = double(firstArg);
                    spatialRef = createSpatialReferencing(resolutionLevelSizes);
                end
                
                % Create writable binary adapter
                obj.CleanUpAdapterDir = true;
                binAdapter = images.internal.adapters.BinAdapter(tempname,'w');               
                % Use this to set internal metadata, so getRegion can
                % return unloaded value
                for ind = 1:numel(spatialRef)
                    binAdapter.appendMetadata(spatialRef(ind).ImageSize,...
                        obj.getBlockSize(ind), results.NumChannels, results.DataType, []);
                end
                obj.Adapter = binAdapter;
                obj.ClassUnderlying = obj.Adapter.PixelDatatype;
                obj.SpatialReferencing = spatialRef;
                
            else % Readable adapters
                if isa(firstArg, 'images.internal.adapters.BigImageAdapter')
                    % Undocumented Adapter input syntax
                    obj.Adapter = firstArg;
                else
                    obj.Adapter = pickAdapter(firstArg, 'r');
                end
                obj.ClassUnderlying = obj.Adapter.PixelDatatype;

                usingDefaultBlockSize = contains('BlockSize',parser.UsingDefaults);
                if usingDefaultBlockSize
                    % Default block size is either [1024 1024] with no adapter
                    % or the closest multiple of underlying file to get to 1024
                    % or if underlying block is larger, pick that
                    optBlockSize = [1024 1024]; % for good performance
                    factor = floor(optBlockSize./obj.Adapter.IOBlockSize);
                    factor(factor<1)=1;
                    obj.BlockSize = obj.Adapter.IOBlockSize.*factor;
                else
                    obj.BlockSize = results.BlockSize;
                end
            end
            
            if contains('UnloadedValue', parser.UsingDefaults)
                % Cast it
                if (obj.ClassUnderlying == "categorical")
                    % <undefined>
                    results.UnloadedValue = categorical(0, 1);
                else
                    results.UnloadedValue = cast(results.UnloadedValue, obj.ClassUnderlying);
                end
            end
            obj.UnloadedValue = results.UnloadedValue;
            
            if ~contains('SpatialReferencing', parser.UsingDefaults)
                % Name-Value trumps default
                obj.SpatialReferencing = results.SpatialReferencing;
            end                        
        end
        
        function delete(obj)
            if isa(obj.Adapter,'images.internal.adapters.BinAdapter')...
                    && obj.CleanUpAdapterDir
                % If we created the adapter, delete its content.
                rmdir(obj.Adapter.DataSource,'s');
            end
        end
        
        function varargout = apply(obj, level, fcn, varargin)
            %apply   Process image
            %   NEWIMG = apply(IMG, LEVEL, FCN) calls the function with the
            %   handle FCN on all of the blocks at the specified resolution
            %   LEVEL. Processing happens one block at a time. The result is
            %   a new bigimage object NEWIMG. This new image only exists in
            %   a temporary form as a set of binary files until the object
            %   is explicitly written. The value of NEWIMG.ClassUnderlying
            %   is the type of the result of FCN, which is a function that
            %   takes a block and returns a block of the same size. NEWIMG
            %   has one resolution level.
            %
            %   [IMG1, IMG2, ...] = apply(___) yields new bigimage objects
            %   IMG1, IMG2, and so on for each output of FCN that has the
            %   same size as the input block.
            %
            %   [IMG, IMG2, ..., IMGN, OTHER] = apply(___) yields N
            %   bigimage objects and a table or cell array of OTHER outputs
            %   from FCN. The outputs of FCN that are the same size as the
            %   input block will become bigimage objects. If the other
            %   output from FCN is a struct, OTHER will become a table of
            %   results for each block. If the other outputs of FCN do not
            %   have the same size as the input block and is not a struct,
            %   then OTHER will be a cell array of the extra outputs. The
            %   last column of OTHER contains the world coordinates of the
            %   center of the top left pixel of each block of the output
            %   (Note: the order of processing the blocks is not
            %   guaranteed).
            %
            %   OTHER = apply(___) produces either a table or cell array of
            %   extra values if FCN does not produce any image outputs.
            %
            %   [___] = apply(IMG, LEVEL, FCN, EXTRAIMAGES) passes
            %   spatially related blocks from IMG and EXTRAIMAGES at the
            %   given LEVEL to FCN. EXTRAIMAGES is a vector of bigimage
            %   objects and each image in EXTRAIMAGES must have the same
            %   spatial extents as IMG. It is possible to mix images with
            %   different numbers of channels and ClassUnderlying types.
            %   Extra image blocks are provided to FCN in the same order as
            %   EXTRAIMAGES. While the blocks correspond to the same
            %   spatial area as the block from IMG, they need not be of the
            %   same size.
            %
            %   [___] = apply(___, Name, Value) specifies additional
            %   parameters. Supported parameters include:
            %
            %   'UseParallel'         A logical scalar specifying if a new
            %                         or existing parallel pool should be
            %                         used. If no parallel pool is active,
            %                         a new pool is opened based on the
            %                         default parallel settings. This
            %                         syntax requires Parallel Computing
            %                         Toolbox. The DataSource property of
            %                         all input bigimages should be valid
            %                         paths on each of the parallel
            %                         workers. Specify OutputFolder if
            %                         workers do not share the same
            %                         filesystem as the client process. If
            %                         relative paths are used, ensure
            %                         workers and client process are on the
            %                         same logical working directory.
            %                         Default value is false.
            %
            %   'BlockSize'           A 2-element vector, [ROWS COLUMNS]),
            %                         specifying block size passed to FCN.
            %                         All channels are always passed.
            %                         Default value is the BlockSize of IMG
            %                         at the specified LEVEL.
            %
            %   'BorderSize'          A 2-element vector [M N]. Adds M rows
            %                         and N columns from neighboring blocks
            %                         to each block passed to FCN. For
            %                         blocks that lie on the edge of an
            %                         image, data is padded according to
            %                         the value of the 'PadMethod'
            %                         parameter. Default value is [0, 0].
            %
            %   'DisplayWaitbar'      A logical scalar. When set to true,
            %                         a waitbar is displayed for long
            %                         running operations.If the waitbar is
            %                         cancelled, a partial output is
            %                         returned if available. The default
            %                         value is true.
            %
            %   'OutputFolder'        A folder/directory location to save
            %                         the output bigimages. This directory
            %                         will contain bigimage subdirectories,
            %                         one for each bigimage output. Default
            %                         value is obtained by invoking the
            %                         TEMPNAME function, and its contents
            %                         are cleared when the output bigimage
            %                         is cleared. However, contents of this
            %                         directory are not cleared
            %                         automatically if OutputFolder is
            %                         explicitly specified.
            %
            %   'PadMethod'           String specifying the kind of data
            %                         that should be added to blocks that
            %                         are on the edge of the image when
            %                         'BorderSize' contains a nonzero
            %                         value. This parameter has no effect
            %                         on interior blocks or on the sides of
            %                         edge blocks that point toward the
            %                         interior of the image. Accepted
            %                         values for kind are a numeric scalar,
            %                         'replicate', and 'symmetric'. The
            %                         default is 0.
            %
            %   'Mask'                A bigimage object, BIGLOGICAL, limits
            %                         processing to parts of IMG that
            %                         overlap the nonzero regions of
            %                         BIGLOGICAL, which is a bigimage
            %                         object with a ClassUnderlying value of
            %                         logical.
            %
            %   'InclusionThreshold'  A scalar ratio, RATIO, determines
            %                         whether an input block should be
            %                         considered for processing. Mask has
            %                         to be specified. For each input
            %                         block, the spatially corresponding
            %                         area of the Mask is tested against
            %                         the RATIO. RATIO is a double value in
            %                         the range [0, 1] and specifies the
            %                         minimum percentage for inclusion in
            %                         processing. A value of 0 implies that
            %                         only one mask pixel must be
            %                         true/nonzero. A value of 1 requires
            %                         that all mask pixels in the mask
            %                         region be true/nonzero. Default value
            %                         is 0.5.
            %
            %   'ExtraImageLevels'    A vector with the same number of
            %                         elements as the number of extra
            %                         bigimage inputs. These values are the
            %                         levels that will be read. The spatial
            %                         extents at these ExtraImageLevels
            %                         must match those of IMG at LEVEL. By
            %                         default, these values are the same as
            %                         LEVEL. Note: If the underlying image
            %                         sizes do not match, then the blocks
            %                         supplied to FCN will not have the
            %                         same size.
            %
            %   Notes:
            %   ------
            %   When apply is called, FCN is immediately evaluated.
            %   Multiple calls to apply will result in multiple traversals
            %   through the data. Data is passed to apply one block at a
            %   time in an undefined order that is consistent with the most
            %   efficient way to traverse the data in the file. Each block
            %   is guaranteed to be passed to fcn only once.
            %
            %   To optimize processing time, pack as much work as
            %   practicable into the function referenced by FCN. This
            %   minimizes disk I/O and ensures data locality. Processing
            %   time at a particular level can be further reduced by using
            %   masks created at lower resolution levels to exclude large
            %   regions of the image that don't require processing. To
            %   further optimize runtime, use BlockSizes as large as system
            %   memory and intermediate memory need of FCN allow. Its best
            %   to keep this a multiple of the underlying Datasource's
            %   physical block size.
            %
            %   If any of the image dimension at the processing level is
            %   not an integral multiple of the BlockSize, partial blocks
            %   will be returned along that dimension.
            %
            %   When BorderSize is specified and an output bigimage is
            %   desired, FCN must return a block with the border trimmed
            %   out, or a block with the same size as the input. apply
            %   will trim the border if the output block is the same size
            %   as the input.
            %
            %   The mask is not available within FCN. If needed,
            %   additionally pass it in as an ExtraImage.
            %
            %   Example 1
            %   ---------
            %   % Enhance image details to better visualize region boundaries
            %     bim = bigimage('tumor_091R.tif');
            %     benh= apply(bim,1, @(block)imguidedfilter(block,block,'DegreeOfSmoothing', 2000));
            %     ha1 = subplot(1,2,1);
            %     bigimageshow(bim,'ResolutionLevel',1);
            %     ha2 = subplot(1,2,2);
            %     bigimageshow(benh);
            %     linkaxes([ha1, ha2]);
            %     xlim([2100, 2600])
            %     ylim([1800 2300])
            %
            %   Example 2
            %   ---------
            %   % Use a mask to limit region used for processing to improve
            %   % processing efficiency            
            %     bim = bigimage('tumor_091R.tif');
            %
            %     % Create a mask at the coarsest level
            %     bmask = apply(bim, bim.CoarsestResolutionLevel,@(im)rgb2gray(im)<80);  
            %     figure
            %     bigimageshow(bmask)
            %
            %     % Use the mask to limit regions used in apply call
            %     % Note - use bigimageshow and showmask functions to
            %     % interactively determine a useful InclusionThreshold
            %     benh= apply(bim,1, @(block)imguidedfilter(block,block,'DegreeOfSmoothing', 2000),...
            %                 'Mask', bmask,'InclusionThreshold', 0.005);
            %     ha1 = subplot(1,2,1);
            %     bigimageshow(bim,'ResolutionLevel',1);
            %     ha2 = subplot(1,2,2);
            %     bigimageshow(benh);
            %     linkaxes([ha1, ha2]);
            
            %            
            %  See also bigimage, bigimageshow
                        
            validateattributes(level, {'numeric'}, ...
                {'integer', 'positive', 'scalar', '<=', obj.NumLevels}, ...
                'bigimage', 'level', 2)
            
            nargoutRequested = max(nargout,1);  % Support "ans"
            
            if nargin > 3 && isa(varargin{1}, 'bigimage')
                extraImages = varargin{1};
                varargin(1) = [];
            elseif nargin > 3 && isempty(varargin{1})
                error(message('images:bigimage:applyNotImageOrParameter'))
            else
                extraImages = [];
            end
            
            % Parse remainder of input arguments.
            parser = inputParser;
            parser.FunctionName = 'apply';
            parser.CaseSensitive = false;
            parser.PartialMatching = true;
            parser.KeepUnmatched = false;
            parser.addParameter('OutputFolder',tempname, @valdiateOutputFolder)
            parser.addParameter('BlockSize', [], @isTwoElementPositiveInteger)
            parser.addParameter('BorderSize', [0 0], @isTwoElementPositiveInteger)
            parser.addParameter('PadMethod', 0, @validatePadMethod)
            parser.addParameter('Mask', [], @validateMasks)
            parser.addParameter('InclusionThreshold', 0.5, @isZeroToOneScalar)
            parser.addParameter('ShowProgress', false, @isLogicalScalar)
            parser.addParameter('UseParallel', false, @validateUseParallel)
            parser.addParameter('DisplayWaitbar', true, @isLogicalScalar)
            parser.addParameter('ExtraImageLevels', [], @(extraLevels) obj.validateExtraImageLevels(extraLevels, extraImages, level))
            parser.parse(varargin{:});
            
            % Parse masks
            masks = parser.Results.Mask;
            obj.checkMaskCompatibility(masks, level)
            obj.validateExtraImages(extraImages, parser.Results.ExtraImageLevels, level)
            
            inclusionThreshold = parser.Results.InclusionThreshold;
            
            
            if any(strcmp(parser.UsingDefaults, 'BlockSize'))
                args.BlockSize = obj.getBlockSize(level);
            else
                args.BlockSize = parser.Results.BlockSize;
            end
            
            args.BorderSize = parser.Results.BorderSize;
            args.PadMethod = parser.Results.PadMethod;
            args.InclusionThreshold = inclusionThreshold;
            args.NargoutRequested = nargoutRequested;
            
            if any(strcmp(parser.UsingDefaults, 'ExtraImageLevels'))
                eLevels = repmat(level, numel(extraImages));
            else
                eLevels = parser.Results.ExtraImageLevels;
            end
            
            args.Levels = [level eLevels];
            args.OutputFolder = parser.Results.OutputFolder;                        
            args.TimeStamp = datestr(now,'ddd_HH_MM_SS');
            
            % Create output folder            
            [flag, fmsg] = mkdir(args.OutputFolder);
            if ~flag
                error(message('images:bigimage:coultNotCreateOutputDir', args.OutputFolder, fmsg));
            end            
            
            % Wait bar variables ------------------------------------------
            obj.WaitBarVars.waitBar = []; % handle to the hg waitbar object
            
            obj.WaitBarVars.NumCompleted = 0;
            obj.WaitBarVars.Total = 0;
            for lvlInd = 1:numel(level)
                % Compute total number of blocks across all requested
                % levels
                numBlocks = prod(ceil(obj.SpatialReferencing(level(lvlInd)).ImageSize ./ args.BlockSize));
                obj.WaitBarVars.Total = obj.WaitBarVars.Total+numBlocks;
            end
            
            % Update wait bar for first few blocks and then per percentage increment
            obj.WaitBarVars.updateIncrements = unique([1:20 round((0.01:0.01:1) .*obj.WaitBarVars.Total)]);
            obj.WaitBarVars.updateCounter = 1;
            obj.WaitBarVars.startTic = tic;
            
            obj.WaitBarVars.Enabled = parser.Results.DisplayWaitbar;
            obj.WaitBarVars.Aborted = false;
            
            % Callback for when each block is done processing
            updateWaitBarFcn = @()obj.updateWaitbar();
            % Callback to close wait bar on completion or on CTRL+C
            cleanUpFcn = onCleanup(@()obj.deleteWaitBar());
            %--------------------------------------------------------------
            
            
            if parser.Results.UseParallel
                % Set up parallel pool.
                p = gcp();
                if isempty(p)
                    error(message('images:bigimage:couldNotOpenPool'))
                end
                numWorkers = p.NumWorkers;
                
                locations = {obj.DataSource};
                if ~isempty(extraImages)
                    locations = [locations(:)', {extraImages.DataSource}];
                end
                if ~isempty(masks)
                     locations = [locations(:)', {masks.DataSource}];
                end
                validatePathsOnWorkers(p,locations, args.OutputFolder);
                
                % A data queue for the workers to single completion of one
                % block, used to trigger an update to the waitbar
                obj.DoneDataQueue = parallel.pool.DataQueue;
                afterEach(obj.DoneDataQueue, @(~)updateWaitBarFcn());
                
                % Put images and extra data on the workers.
                CMainImg = parallel.pool.Constant(obj);
                CExtraImg = parallel.pool.Constant(extraImages);
                CMasks = parallel.pool.Constant(masks);
                
                % Run fcn on each of the parallel workers.
                obj.Futures = parallel.FevalFuture.empty();
                for idx = numWorkers:-1:1
                    obj.Futures(idx) = parfeval(p, @applyWorkerParallel, ...
                        nargoutRequested, idx, numWorkers, fcn, CMainImg,...
                        CExtraImg, CMasks, args, obj.DoneDataQueue);
                end
                
                % In case user cancels before any data is processed
                out = cell(1, nargoutRequested);
                
                % Fetch data and merge the output values.
                for idx = 1:numWorkers
                    results = cell(1, nargoutRequested);
                    
                    try
                        % To fetch completed results one by one
                        [~, results{:}] = fetchNext(obj.Futures);
                    catch ERR
                        if ~isempty(ERR.cause) && strcmp(ERR.cause{1}.identifier, 'parallel:fevalqueue:ExecutionCancelled')
                            % Did the user cancel? if so, return
                            % partial results
                            break;
                        else
                            % Unexpected error, clean up and inform user
                            cancel(obj.Futures);
                            % Wait till all running onces are completed
                            % (else they may bring up a wait bar which wont
                            % get deleted)
                            wait(obj.Futures);
                            delete(obj.DoneDataQueue);
                            rethrow(ERR)
                        end
                    end
                    if idx == 1
                        out = results;
                    else
                        out = merge(out, results);
                    end
                end
                
                varargout = out;
            else
                % Call fcn for all data at once on the host.
                varargout = cell(1, nargoutRequested);
                isWaitBarAbortedFcn = @()obj.WaitBarVars.Aborted;
                [varargout{:}] = applyWorkerSequential(fcn, obj, extraImages, masks, args, updateWaitBarFcn, isWaitBarAbortedFcn);
            end
            
            
            for ind = 1:numel(varargout)
                if isa(varargout{ind}, 'bigimage')
                    % Workers create adapters first to ensure the local
                    % objects don't clean up the local parts. But the
                    % output from this function should clean itself up.
                    varargout{ind}.CleanUpAdapterDir = true;
                    if ~contains('OutputFolder', parser.UsingDefaults)
                        % Unless OutputFolder is specified explicitly
                        varargout{ind}.CleanUpAdapterDir = false;
                    end
                    % Overwrite the default blocksize to explicitly be the
                    % block size that was used for its creation.
                    varargout{ind}.BlockSize = args.BlockSize;
                end
            end
            
        end
        
        function tf = isequal(obj, otherObj)
            %isequal   Compare equality of two bigimage objects
            %  TF = isequal(IMG1, IMG2) compares the sizes and contents of
            %  two bigimage object, IMG1 and IMG2. The images are equal if
            %  their pixel data is the same and the following properties
            %  are equal: SpatialReferencing, ClassUnderlying, and
            %  BlockSize.
            %
            %  Pixel values are only compared if the metadata is equal, and
            %  comparison ends early when unequal pixel values are found.
            
            tf = isequal(obj.SpatialReferencing, otherObj.SpatialReferencing) && ...
                isequal(obj.Channels, otherObj.Channels) && ...
                isequal(obj.ClassUnderlying, otherObj.ClassUnderlying) && ...
                isequal(obj.BlockSize, otherObj.BlockSize);
            
            if tf
                % Metadata is the same. Interrogate each block.
                if obj.CoarsestResolutionLevel == obj.FinestResolutionLevel
                    allLevels = 1;
                else
                    % Always compare coarsest first and finest last.
                    allLevels = 1:obj.NumLevels;
                    allLevels(allLevels == obj.CoarsestResolutionLevel) = [];
                    allLevels(allLevels == obj.FinestResolutionLevel) = [];
                    allLevels = [obj.CoarsestResolutionLevel, allLevels, obj.FinestResolutionLevel];
                end
                
                for oneLevel = allLevels
                    ds1 = bigimageDatastore(obj, oneLevel);
                    ds2 = bigimageDatastore(otherObj, oneLevel);
                    while tf && hasdata(ds1) && hasdata(ds2)
                        tf = isequal(read(ds1), read(ds2));
                    end
                end
            end
        end
        
        function data = getRegion(obj, level, regionStartWorld, regionEndWorld)
            %getRegion   Get arbitrary portion of image
            %   PIXELS = getRegion(IMG, LEVEL, REGIONSTART, REGIONEND)
            %   returns all PIXELS at the specified LEVEL of IMG which have
            %   their center location between the points REGIONSTART and
            %   REGIONEND, inclusive. These points are (x,y) pairs of world
            %   coordinates.
            %
            %   Example
            %   -------
            %   % Read a sub-region of a bigimage
            %     bim = bigimage('tumor_091R.tif');
            %     xlimits = [2100 2600];
            %     ylimits = [1800, 2300];
            %     imL1 = getRegion(bim, 1,[xlimits(1), ylimits(1)], [xlimits(2), ylimits(2)]);
            %     imL2 = getRegion(bim, 2,[xlimits(1), ylimits(1)], [xlimits(2), ylimits(2)]);
            %     imL3 = getRegion(bim, 3,[xlimits(1), ylimits(1)], [xlimits(2), ylimits(2)]);
            %     montage({imL1, imL2, imL3},'Size', [1 3])
            %
            %   See also bigimage
            
            
            validateattributes(level, {'numeric'}, {'positive', 'scalar', 'integer', '<=', obj.NumLevels}, 'bigimage', 'level')
            validateattributes(regionStartWorld, {'numeric'}, {'vector', 'numel', 2}, 'bigimage', 'regionStart')
            validateattributes(regionEndWorld, {'numeric'}, {'vector', 'numel', 2}, 'bigimage', 'regionEnd')
            if ~all(regionEndWorld>=regionStartWorld)
                error(message('images:bigimage:endShouldBeGreaterThanOrEqualToStart'));
            end
            
            % Convert to r,c format
            regionStartWorld = [regionStartWorld(2), regionStartWorld(1)];
            regionEndWorld = [regionEndWorld(2), regionEndWorld(1)];
            
            ref = obj.SpatialReferencing(level);
            if any(regionStartWorld<[ref.YWorldLimits(1), ref.XWorldLimits(1)])...
                    ||any(regionEndWorld>[ref.YWorldLimits(2), ref.XWorldLimits(2)])
                error(message('images:bigimage:invalidLocation'));
            end
            
            [startIntrinsicR, startIntrinsicC ]= obj.SpatialReferencing(level).worldToSubscriptAlgo(regionStartWorld(2), regionStartWorld(1));
            [endIntrinsicR, endIntrinsicC ]= obj.SpatialReferencing(level).worldToSubscriptAlgo(regionEndWorld(2), regionEndWorld(1));
            
            data = obj.Adapter.readRegion(level, [startIntrinsicR, startIntrinsicC], [endIntrinsicR, endIntrinsicC]);
            if strcmp(obj.ClassUnderlying, 'categorical')
                data = categorical(data, obj.PixelLabelIDs, cellstr(obj.ClassNames));
            end
        end
        
        function data = getBlock(obj, level, locationWorld)
            %getBlock   Read specific block of image
            %   PIXELS = getBlock(IMG, LEVEL, LOCATION) reads the PIXELS of
            %   the block at the specified LEVEL of IMG that contains
            %   LOCATION, which is a (x,y) point specified in world
            %   coordinates.
            %
            %   Example
            %   -------
            %   % Read the block that belongs to the chosen location
            %     bim = bigimage('tumor_091R.tif');
            %     subplot(1,2,1)
            %     bigimageshow(bim, 'GridVisible','on','GridLevel', 1);
            %     title('Move point to select block');
            %     hp = drawpoint('Position', [1,1]);
            %     ha = subplot(1,2,2);
            %     addlistener(hp, 'ROIMoved', @(~,~)...
            %                 imshow(bim.getBlock(1, hp.Position), 'Parent', ha));
            %                   
            %   See also bigimage
            
            blockStartIntrinsic = obj.worldStartToBlockOriginInstrinsic(level, locationWorld);
            
            % Compute block end, clamp
            ref = obj.SpatialReferencing(level);
            bsize = obj.getBlockSize(level);
            blockEndIntrinsic = blockStartIntrinsic + bsize-1;
            blockEndIntrinsic = min(blockEndIntrinsic, ref.ImageSize);
            
            data = obj.Adapter.readRegion(level, blockStartIntrinsic, blockEndIntrinsic);
            
            if strcmp(obj.ClassUnderlying, 'categorical')
                data = categorical(data, obj.PixelLabelIDs, cellstr(obj.ClassNames));
            end
        end
        
        function image = getFullLevel(obj, level)
            %getFullLevel  Get full image at specified location
            %  IMAGE = getFullLevel(IMG) reads the full image at the
            %  coarsest resolution level of IMG. Use this method to extract
            %  a low resolution image from IMG for use with Apps.
            %
            %  IMAGE = getFullLevel(IMG, LEVEL) reads the full image at
            %  specified level.
            %
            %  Note: Check the LevelSizes property to ensure the level
            %  being read from is small enough to fit in system memory.
            %
            %   Example
            %   -------
            %   % Read a full level into memory and create a bigimage out of
            %   % it
            %     bim = bigimage('tumor_091R.tif');
            %
            %     % Create a mask at the coarsest level
            %     clevel = bim.CoarsestResolutionLevel;
            %     % Work on the full image since imbinarize relies on
            %     % global values
            %     imcoarse = getFullLevel(bim, clevel);
            %     stainMask = ~imbinarize(rgb2gray(imcoarse));
            %     % Retain the original spatial referencing information
            %     bmask = bigimage(stainMask,...
            %                'SpatialReferencing', bim.SpatialReferencing(clevel));
            %     figure
            %     bigimageshow(bmask);
            
            if nargin==1
                level = obj.CoarsestResolutionLevel;
            end
            ref = obj.SpatialReferencing(level);
            image = obj.getRegion(level,...
                [ref.XWorldLimits(1), ref.YWorldLimits(1)],...
                [ref.XWorldLimits(2), ref.YWorldLimits(2)]);
        end
        
        function setBlock(obj, level, locationWorld, data)
            % setBlock Set block content
            %  setBlock(IMG, LEVEL, LOCATION, DATA) sets the block content
            %  at the specified LEVEL of IMG that contains LOCATION.
            %  LOCATION is a (x,y) point specified in world coordinates.
            %  DATA is a M-by-N-by-Channels numeric data whose first two
            %  dimensions match the corresponding levels BlockSize. DATA
            %  must have the same type as ClassUnderlying.
            %
            %  Notes:
            %  setBlock can only be called on bigimages that were
            %  created with one of the writable constructor syntax.
            %  If size of DATA is less than BlockSize, it is padded with
            %  UnloadedValue.
            %  setBlock will internally trim data for partial edge blocks.
            %
            %   Example
            %   -------
            %   % Create a mask from an ROI
            %     bim = bigimage('tumor_091R.tif');
            %
            %     % Create an ROI
            %     h = bigimageshow(bim);
            %     hROI = drawcircle(gca, 'Radius', 470, 'Position', [1477 2284]);
            %
            %     % Choose level at which to create the mask
            %     maskLevel = 3;
            %
            %     ref = bim.SpatialReferencing(maskLevel);
            %     pixelExtent = [ref.PixelExtentInWorldX, ref.PixelExtentInWorldY];
            %
            %     % Create a writable bigimage
            %     bmask = bigimage(ref, 1, 'logical');
            %
            %     for cStart = 1:bmask.BlockSize(2):ref.ImageSize(2)
            %       for rStart = 1:bmask.BlockSize(1):ref.ImageSize(1)
            %           % Center of top left pixel of this block in world units
            %           xyStart = [cStart, rStart].* pixelExtent;
            %
            %           % Center of bottom right pixel of this block in world units
            %           bsize = fliplr(bmask.BlockSize); % r,c->x,y
            %           xyEnd = ([cStart, rStart] + (bsize-1)).* pixelExtent;
            %
            %           % Mesgrid in world units for all pixels in this block.
            %           [xgrid, ygrid] = meshgrid(xyStart(1):ref.PixelExtentInWorldX:xyEnd(1),...
            %               xyStart(2):ref.PixelExtentInWorldY:xyEnd(2));
            %
            %           % Create in/out mask for this block
            %           roiPositions = hROI.Vertices;
            %           tileMask = inpolygon(xgrid, ygrid,...
            %               roiPositions(:,1), roiPositions(:,2));
            %
            %           % Write out the block
            %           setBlock(bmask, 1, xyStart, tileMask);
            %       end
            %     end
            %
            %     figure
            %     bigimageshow(bmask)
            %
            % See also bigimage
            validateattributes(data, obj.ClassUnderlying,...
                {'nonempty'},mfilename, 'DATA');
            if size(data,3) ~= obj.Channels
                error(message('images:bigimage:dataChannelsIncorrect', obj.Channels));
            end
            
            if obj.Adapter.Mode ~= 'w'
                error(message('images:bigimage:readOnlyAdapter'));
            end
            
            blockStartIntrinsic = obj.worldStartToBlockOriginInstrinsic(level, locationWorld);
            
            bsize = obj.getBlockSize(level);
            if size(data,1) > bsize(1) || ...
                    size(data,2) > bsize(2) || ...
                    size(data,3) ~= obj.Channels
                error(message('images:bigimage:amountOfBlockData'))
            end
            
            obj.setBlockIntrinsic(level, blockStartIntrinsic, data);
        end
        
        function write(obj, filename)
            %write   Write image contents to new file
            %   write(IMG, TIFFFILE) writes the contents of IMG out as a
            %   TIFF file with LZW (lossless) compression. TIFFFILE is
            %   expected to have a .tif or .tiff extension.
            %
            %   write(IMG, DIRNAME) writes a binary block version of IMG to
            %   the directory DIRNAME. DIRNAME is created and should not
            %   already exist. A new bigimage can be created from this
            %   directory using the bigimage(DIRNAME) syntax.
            %
            %   Notes:
            %   To preserve custom SpatialReferencing, use the
            %   DIRNAME syntax.
            %   Tiff files are written using the LZW loseless compression
            %   scheme.
            %
            %   Example
            %   -------
            %   % Write out a bigimage to disk
            %     bim = bigimage('tumor_091R.tif');
            %     mask = apply(bim, 3, @(im)rgb2gray(im)<100);
            %    
            %     % Write to a TIFF file
            %     write(mask, 'maskOut.tif')
            %
            %     % Write to a bigimage directory
            %     write(mask, 'maskDir/');
            %
            %     % Load it back
            %     bm1 = bigimage('maskOut.tif');
            %     bm2 = bigimage('maskDir');            
            %
            %     % TIFF files do not retain the spatial referencing
            %     % information, hence switch to the default
            %     disp(bm1.SpatialReferencing)
            %
            %     % Restore from the one saved in the bigimage directory,
            %     % required to ensure isequal call works on the same
            %     % region.
            %     bm1.SpatialReferencing = bm2.SpatialReferencing
            %
            %     disp(isequal(bm1, bm2))
            %
            %  See also bigimage
            
            if isstring(filename) && isscalar(filename)
                filename = char(filename);
            end
            validateattributes(filename,{'char', 'string'}, {'nonempty'},...
                mfilename,'DIRNAME/TIFFFILE');
            
            % no-clobber
            if ~isempty(dir(filename))
                error(message('images:bigimage:destinationexists'));
            end
            
            % Get adapter ready.
            if contains(filename,'.tif')
                % Assume TIFFFILE
                writeAdapter = images.internal.adapters.TIFFAdapter(filename, 'w');
            else
                % Assume DIR
                writeAdapter = images.internal.adapters.BinAdapter(filename, 'w');
                writeAdapter.UnloadedValue = obj.UnloadedValue;
                writeAdapter.SpatialReferencing = obj.SpatialReferencing;
            end
            
            if isa(obj.Adapter,'images.internal.adapters.BinAdapter') ...
                    && isa(writeAdapter,'images.internal.adapters.BinAdapter')
                % Copy the source directory
                copyfile(obj.Adapter.DataSource, writeAdapter.DataSource);
                writeAdapter.finalizeWrite();
                return
            end
            
            try
                for lInd = 1:obj.NumLevels
                    bSize = obj.getBlockSize(lInd);
                    
                    if isa(writeAdapter,'images.internal.adapters.TIFFAdapter')
                        % Required block size to be a multiple of 16
                        bSize = ceil(bSize/16)*16;
                        if ~isequal(bSize, obj.getBlockSize(lInd))
                            warning(message('images:bigimage:changedBlockSize'));
                        end
                    end
                    
                    % For each level
                    imageSize = obj.SpatialReferencing(lInd).ImageSize;
                    writeAdapter.appendMetadata(imageSize, ...
                        bSize, obj.Channels, obj.ClassUnderlying, []);
                    
                    numblocksRC = ceil(imageSize ./ bSize);
                    
                    for r = 1:numblocksRC(1)
                        for c = 1:numblocksRC(2)
                            % Each raw block
                            blockStart = ([r c]-1).*bSize+[1 1];
                            blockEnd = blockStart + bSize -1;
                            % clamp
                            blockEnd = min(blockEnd, obj.LevelSizes(lInd,:));
                            data = obj.getRegionIntrinsic(lInd, blockStart, blockEnd);       
                            writeAdapter.writeBlock(lInd, blockStart, data)
                        end
                    end
                end                
                writeAdapter.finalizeWrite();
                
            catch ALL
                % Write failed.
                % Clean up partially written stuff
                clear writeAdapter
                if exist(filename,'file') 
                    delete(filename);
                end
                if exist(filename, 'dir')                    
                    rmdir(filename, 's');
                end
                rethrow(ALL);
            end
        end
        
        function s = saveobj(obj)
            s.BlockSize = obj.BlockSize;
            s.UnloadedValue = obj.UnloadedValue;
            s.ClassNames = obj.ClassNames;
            s.PixelLabelIDs = obj.PixelLabelIDs;
            s.Adapter = obj.Adapter;
            s.ClassUnderlying = obj.ClassUnderlying;
            s.SpatialReferencing = obj.SpatialReferencing;
        end
    end
    
    % Property access methods
    methods
        function set.Adapter(obj, value)
            validateattributes(value,...
                {'images.internal.adapters.BigImageAdapter'},...
                {'nonempty', 'scalar'}, mfilename, 'Adapter');
            obj.Adapter = value;           
            obj.ClassUnderlying = obj.Adapter.PixelDatatype;
            % Reload props (if) saved with the adapter
            if isprop(obj.Adapter, 'UnloadedValue') && ~isempty(obj.Adapter.UnloadedValue)
                obj.UnloadedValue = obj.Adapter.UnloadedValue; %#ok<MCSUP>
            end
            if isprop(obj.Adapter, 'SpatialReferencing') && ~isempty(obj.Adapter.SpatialReferencing)
                obj.SpatialReferencing = obj.Adapter.SpatialReferencing; %#ok<MCSUP>
            else
                % No, create defaults from image extents
                resolutionLevelSizes = [obj.Adapter.Height(:), obj.Adapter.Width(:)];
                obj.SpatialReferencing = createSpatialReferencing(resolutionLevelSizes);%#ok<MCSUP>
            end
        end
        
        function set.SpatialReferencing(obj, value)
            validateattributes(value, {'imref2d'}, ...
                {'nonempty', 'vector'}, 'bigimage', 'SpatialReferencing')
            
            adapter = obj.Adapter;%#ok<MCSUP>
            if ~isempty(adapter)
                % Validate against actual data
                validateattributes(value, {'imref2d'}, ...
                    {'nonempty', 'numel', numel(adapter.Height)}, mfilename, 'SpatialReferencing');
                for lInd = 1:numel(adapter.Height)
                    if ~isequal(value(lInd).ImageSize, [adapter.Height(lInd), adapter.Width(lInd)])
                        error(message('images:bigimage:imageSizeDoNotMatch'));
                    end
                end
            end           
            
            obj.SpatialReferencing = value;
        end
        
        function set.UnloadedValue(obj, value)
            
            dataType = char(obj.ClassUnderlying); %#ok<MCSUP>
            channels = obj.Channels; %#ok<MCSUP>
            
            % Either scalar or numChannels in numel
            attribute = {'scalar'};
            if ~isscalar(value)
                attribute = {'numel', channels};
            end
            validateattributes(value, {dataType}, attribute, 'bigimage');
            
            if isscalar(value) && channels>1
                % Convert to 1x1xchannels
                obj.UnloadedValue = repmat(value, [1, 1, channels]);
            else % already numChannels, make it right shape
                obj.UnloadedValue = reshape(value, [1, 1, channels]);
            end
            
            if isa(obj.Adapter, 'images.internal.adapters.BinAdapter') %#ok<MCSUP>
                % need to return 'virtual' blocks
                obj.Adapter.UnloadedValue = obj.UnloadedValue;%#ok<MCSUP>
            end
        end
    end
    
    % Get methods
    methods
        function ds = get.DataSource(obj)
            ds = obj.Adapter.DataSource;
        end
        
        function l = get.NumLevels(obj)
            l = numel(obj.SpatialReferencing);
        end
        
        function LevelSizes = get.LevelSizes(obj)
            LevelSizes = reshape([obj.SpatialReferencing.ImageSize],[2, obj.NumLevels])';
        end
        
        function numChannels = get.Channels(obj)
            numChannels = obj.Adapter.Channels;
        end
        
        function details = get.SourceDetails(obj)
            details = obj.Adapter.SourceMetadata;
        end
                       
        function l = get.CoarsestResolutionLevel(obj)
            px = [obj.SpatialReferencing.PixelExtentInWorldX];
            py = [obj.SpatialReferencing.PixelExtentInWorldY];
            parea = px.*py;
            [~, l] = max(parea);
        end
        
        function l = get.FinestResolutionLevel(obj)
            px = [obj.SpatialReferencing.PixelExtentInWorldX];
            py = [obj.SpatialReferencing.PixelExtentInWorldY];
            parea = px.*py;
            [~, l] = min(parea);
        end
    end
    
    % Hidden helper functions
    % These are undocumented, hidden syntaxes for internal use only.
    % Their implementations may change in future releases. Do not use
    % them or rely on their behavior.
    methods (Access = public, Hidden = true)
        
        function blockSize = getBlockSize(obj, level)
            % Blocksize is either levelx2 of 1x2 (same block size for all
            % levels)
            level = min(size(obj.BlockSize,1), level);
            blockSize = obj.BlockSize(level,:);
        end
        
        function tf = isLocationInImage(obj, level, location)
            refObj = obj.SpatialReferencing(level);
            
            tf = contains(refObj,location(2),location(1));
            
        end
        
        function addSubLevel(obj,  fromLevel, atImageSize)
            % addSubLevel Add a new level resampled from an existing level
            %   addSubLevel(fromLevel, atImageSize) creates a new level
            %   image level from an existing level at the specified image
            %   size. atImageSize is of the form [Rows Cols]. Either one
            %   Rows or Cols must be NaN. addSubLevel will attempt to create a
            %   new level with image size as close as possible to the
            %   specified non-nan dimension while retaining the reference
            %   levels aspect ratio.
            %
            %   Note: Nearest neighbor interpolation is used to perform the
            %   downsampling.
            %
            
            assert(isa(obj.Adapter,'images.internal.adapters.BinAdapter'));
            
            validateattributes(fromLevel, {'numeric'}, ...
                {'scalar', 'positive', '<=', numel(obj.SpatialReferencing)},...
                mfilename, 'fromLevel');
            
            validateattributes(atImageSize, {'numeric'},...
                {'numel', 2, 'positive'}, mfilename, 'atImageSize');
            if ~(sum(isnan(atImageSize))==1)
                error(message('images:bigimage:oneShouldbeNaN'));
            end
            
            refRef = obj.SpatialReferencing(fromLevel);
            
            % Resolve nan, keep aspect ratio as close as possible
            nonNanInd = ~isnan(atImageSize);
            factor = ceil(refRef.ImageSize(nonNanInd)/atImageSize(nonNanInd));
            scale = 1/factor;
            if factor<=1
                error(message('images:bigimage:onlyDownSize'));
            end
            
            
            % Refine image size
            atImageSize = round(refRef.ImageSize/factor);
            
            % Create new spatial ref for this level
            newRef = refRef;
            newRef.ImageSize = atImageSize;
            % Enforce the same world extents as the reference level
            newRef.XWorldLimits = refRef.XWorldLimits;
            newRef.YWorldLimits = refRef.YWorldLimits;
            
            % Update block size
            if size(obj.BlockSize,1)~=1
                % Pick block size of previous level
                obj.BlockSize(end+1,:) = obj.BlockSize(end,:);
            end
            newLevel = numel(obj.SpatialReferencing)+1;
            
            numblocksRC = obj.getNumberOfblocks(fromLevel);
            
            numDone = 0; %
            waitBar = iptui.cancellableWaitbar(getString(message('images:bigimage:addingNewLevel')),...
                getString(message('images:bigimage:processingNBlocks')),prod(numblocksRC),numDone);
            deleteWaitBar = onCleanup(@()waitBar.destroy);
            
            obj.Adapter.appendMetadata(newRef.ImageSize, obj.BlockSize(end,:), obj.Channels, obj.ClassUnderlying, []);
            
            % Set up an iteration scheme that goes through all the source
            % blocks required for one destination block, in steps of resize
            % factor.
            rDst = 1;
            for r = 1:factor:numblocksRC(1)
                cDst=1;
                for c = 1:factor:numblocksRC(2)
                    outblockInCells = cell.empty();
                    % Loop through source blocks required for one
                    % destination block
                    for rst = r: min(numblocksRC(1), r+factor-1)
                        for cst = c: min(numblocksRC(2), c+factor-1)
                            if waitBar.isCancelled()
                                break;
                            end
                            
                            startIntrinsic = ([rst cst]-1).*obj.getBlockSize(fromLevel)+[1 1];
                            endIntrinsic = startIntrinsic+ obj.getBlockSize(fromLevel)-1;
                            % clamp
                            endIntrinsic = min(endIntrinsic, refRef.ImageSize);
                            
                            loadedblock = obj.Adapter.readRegion(fromLevel, startIntrinsic, endIntrinsic);
                            
                            if obj.ClassUnderlying == "categorical"
                                loadedblock = uint8(loadedblock);
                            end
                            outblockInCells{rst-r+1, cst-c+1} = imresize(loadedblock, scale, 'nearest');
                            numDone = numDone+1;
                        end
                        waitBar.update(numDone);
                    end
                    
                    % Assembly a single destination block from all the resized
                    % source blocks
                    outblock = cell2mat(outblockInCells);
                    
                    % Pad any dimension that is smaller (due to partial
                    % source blocks)
                    outblockSize = [size(outblock,1), size(outblock,2)];
                    if any(outblockSize < obj.BlockSize(end,:))
                        reqPadding = obj.BlockSize(end,:) - outblockSize;
                        reqPadding(reqPadding<0) = 0;
                        paddedblockSize = max([size(outblock,1), size(outblock,2)], obj.BlockSize(end,:));
                        paddedblockSize(3) = size(outblock,3);
                        paddedblock = zeros(paddedblockSize, 'like', outblock);
                        for cInd = 1:size(outblock,3)
                            paddedblock(:,:,cInd) = ...
                                padarray(outblock(:,:,cInd), reqPadding, obj.UnloadedValue(1,1,cInd), 'post');
                        end
                        outblock = paddedblock;
                    end
                    
                    % Shrink any dimension that is larger (this should only
                    % be 1 pixel)
                    outblockSize = [size(outblock,1), size(outblock,2)];
                    if any(outblockSize > obj.BlockSize(end,:))
                        % Can happen if our factor (when ceiled) results in
                        % subblocks that are slightly larger the block size
                        % when assembled together
                        outblock = imresize(outblock, obj.BlockSize(end,:), 'nearest');
                    end
                    
                    blockOrigin = ([rDst cDst]-1).*obj.getBlockSize(newLevel)+[1 1];
                    obj.setBlockIntrinsic(newLevel, blockOrigin, outblock);
                    cDst = cDst+1;
                end
                rDst = rDst+1;
            end
            
            % Add new spatial ref entry (after adding the new level)
            obj.SpatialReferencing(end+1) = newRef;
        end        
        
        function blockStartIntrinsic = worldStartToBlockOriginInstrinsic(obj, level, locationWorld)
            validateattributes(level, {'numeric'}, {'positive', 'scalar', 'integer', '<=', obj.NumLevels}, 'bigimage', 'level')
            validateattributes(locationWorld, {'numeric'}, {'vector', 'numel', 2}, 'bigimage', 'location')
            
            % Convert to r,c format
            locationWorld = [locationWorld(2), locationWorld(1)];
            
            if (level < 1) || (level > obj.NumLevels)
                error(message('images:bigimage:invalidResolutionLevel', obj.NumLevels))
            end
            
            if ~obj.isLocationInImage(level, locationWorld)
                error(message('images:bigimage:invalidLocation'))
            end
            
            ref = obj.SpatialReferencing(level);
            bsize = obj.getBlockSize(level);
            [r,c] = ref.worldToSubscriptAlgo(locationWorld(2), locationWorld(1));
            
            % Find block start of this point
            inWholeBlockUnits = floor([r,c]./bsize);
            isExactMultiple = rem([r,c], bsize)==[0 0];
            if any(isExactMultiple)
                % Point falls on exact edge of block boundary, which belong
                % to the top-left of this block
                inWholeBlockUnits(isExactMultiple) = inWholeBlockUnits(isExactMultiple)-1;
            end
            
            wholeBlockEndBeforeThisPoint = inWholeBlockUnits.*bsize;
            blockStartIntrinsic = wholeBlockEndBeforeThisPoint+1;
        end
        
        function pctNNZ = computeWorldRegionNNZ(obj, level, maskStartWorld, maskEndWorld)
            % Given world extents of a region, compute the pct NNZ in that
            % region.            
            %  - fully out of bounds regions are treats as all false
            %  - partially out of bounds regions are clamped before
            %  computing nnz. (this will inflate the pct for these
            %  regions).
                                    
            % API Strictly for masks only
            assert(obj.Channels==1);
            
            maskRef = obj.SpatialReferencing(level);
                       
            % Check for intersection
            outOfBounds = ...
                maskStartWorld(1)>maskRef.XWorldLimits(2)...
                |maskRef.XWorldLimits(1)>maskEndWorld(1)...
                |maskStartWorld(2)>maskRef.YWorldLimits(2)...
                |maskRef.YWorldLimits(1)>maskEndWorld(2);
            
            % Default for out of bounds
            pctNNZ = 0;
            
            if ~outOfBounds
                % Convert intrinsic first (subs will result in NaNs)
                [x,y] = maskRef.worldToIntrinsicAlgo(maskStartWorld(1),maskStartWorld(2));
                maskStartIntrinsic =[x,y];
                [x,y] = maskRef.worldToIntrinsicAlgo(maskEndWorld(1), maskEndWorld(2));
                maskEndIntrinsic =[x,y];
                
                % Convert to subs (r,c) (integers, flipped order of (x,y))
                maskStartSub = round([maskStartIntrinsic(2), maskStartIntrinsic(1)]);
                maskEndSub= round([maskEndIntrinsic(2), maskEndIntrinsic(1)]);
                % Clamp to in-bounds
                maskStartSub = min(max(maskStartSub, [1 1]), maskRef.ImageSize);
                maskEndSub = min(max(maskEndSub, [1 1]), maskRef.ImageSize);
                
                pctNNZ = obj.Adapter.computeRegionNNZ(level,...
                    maskStartSub, maskEndSub);
            end
        end
        
        function data = getRegionPadded(obj, level, regionStartWorld, regionEndWorld, padMethod)
            %getRegion   Get arbitrary portion of the padded image
            %   PIXELS = getRegionPadded(IMG, LEVEL, REGIONSTART, REGIONEND, PADMETHOD)
            %   returns all PIXELS at the specified LEVEL of IMG which have
            %   their center location between the points REGIONSTART and
            %   REGIONEND, inclusive. These points are (x,y) pairs of world
            %   coordinates. Regions that fall outside the image are padded
            %   using padMethod.
            %
            %   See also bigimage
                     
            ref = obj.SpatialReferencing(level);
            imageWorldXLim = ref.XWorldLimits;
            imageWorldYLim = ref.YWorldLimits;
            imageStartWorld = [imageWorldXLim(1), imageWorldYLim(1)];
            imageEndWorld = [imageWorldXLim(2), imageWorldYLim(2)];
            
            % Compute required prepadding (in world coordinates)
            prepad = [0 0];
            if any(regionStartWorld<imageStartWorld)
                prepad = imageStartWorld-regionStartWorld;
                prepad(prepad<0)=0;
            end
            
            % Compute required postpadding
            postpad = [0 0];
            if any(regionEndWorld>imageEndWorld)
                postpad = regionEndWorld-imageEndWorld;
                postpad(postpad<0)=0;
            end
            
            % Clamp and read valid region
            clampedRegionStart = max(regionStartWorld, imageStartWorld);
            clampedRegionEnd = min(imageEndWorld, regionEndWorld);
            data = obj.getRegion(level, clampedRegionStart, clampedRegionEnd);
            
            % Convert pre/post pad region in world coordinates into pixel
            % count
            pixelDims = [ref.PixelExtentInWorldX, ref.PixelExtentInWorldY];
            % Choose ceil so we at least get one pixel. Extra images cannot
            % span _exact_ spatial extents as reference since we dont have
            % a way to support partial pixels.
            prepad = ceil(prepad./pixelDims);
            postpad = ceil(postpad./pixelDims);
            
            % Convert to row/col for use with padarray
            prepad = [prepad(2), prepad(1)];
            postpad = [postpad(2), postpad(1)];
            
            if any(prepad)
                data = padarray(data,prepad, padMethod, 'pre');
            end
            if any(postpad)
                data = padarray(data,postpad, padMethod, 'post');
            end
        end
    end
    
    % API for bigimageDatastore
    methods (Access = public, Hidden = true)
        function data = getRegionIntrinsic(obj, level, regionStartIntrinsic, regionEndIntrinsic)
            data = obj.Adapter.readRegion(level, regionStartIntrinsic, regionEndIntrinsic);
            if strcmp(obj.ClassUnderlying, 'categorical')
                data = categorical(data, obj.PixelLabelIDs, cellstr(obj.ClassNames));
            end
        end
        
    end
    
    methods (Static)
        function obj = loadobj(s)
            obj = bigimage(s);
            obj.ClassUnderlying = s.ClassUnderlying;
            obj.Adapter = s.Adapter;            
            obj.BlockSize = s.BlockSize;
            obj.UnloadedValue = s.UnloadedValue;            
            obj.ClassNames = s.ClassNames;
            obj.PixelLabelIDs = s.PixelLabelIDs;
            obj.SpatialReferencing = s.SpatialReferencing;
        end
    end
    
    % Internal, validation-free implementation helpers
    methods (Access = private)
        
        function setBlockIntrinsic(obj, level, blockOriginIntrinsic, data)
            obj.Adapter.writeBlock(level, blockOriginIntrinsic, data);
        end
        
        function updateWaitbar(obj)
            % Gets called after each block is processed, updates a waitbar
            % at a slower rate
            
            if ~obj.WaitBarVars.Enabled
                return;
            end
            
            % Increment the cached count of completed blocks
            obj.WaitBarVars.NumCompleted = obj.WaitBarVars.NumCompleted + 1;
            
            % Updates are expensive, do so intermittently
            if obj.WaitBarVars.NumCompleted >= obj.WaitBarVars.updateIncrements(obj.WaitBarVars.updateCounter)
                obj.WaitBarVars.updateCounter = obj.WaitBarVars.updateCounter + 1;
                
                if isempty(obj.WaitBarVars.waitBar)
                    % Wait bar not yet shown
                    elapsedTime = toc(obj.WaitBarVars.startTic);
                    % Decide if we need a wait bar or not,
                    remainingTime = elapsedTime / obj.WaitBarVars.NumCompleted * (obj.WaitBarVars.Total - obj.WaitBarVars.NumCompleted);
                    if remainingTime > 15 % seconds
                        if images.internal.isFigureAvailable()
                            obj.WaitBarVars.waitBar = iptui.cancellableWaitbar(getString(message('images:bigimage:waitbarTitleGUI')),...
                                getString(message('images:bigimage:waitbarTitle')),obj.WaitBarVars.Total,obj.WaitBarVars.NumCompleted);
                        else
                            obj.WaitBarVars.waitBar = iptui.textWaitUpdater(getString(message('images:bigimage:waitbarTitle')),...
                                getString(message('images:bigimage:waitbarCompletedTxt')),obj.WaitBarVars.Total,obj.WaitBarVars.NumCompleted);
                        end
                    end
                    
                elseif obj.WaitBarVars.waitBar.isCancelled()
                    % User clicked on 'Cancel'
                    obj.WaitBarVars.Aborted = true;
                    if ~isempty(obj.Futures)
                        % UseParallel == true mode.
                        delete(obj.DoneDataQueue);
                        cancel(obj.Futures);
                    end
                else
                    % Show progress on existing wait bar
                    obj.WaitBarVars.waitBar.update(obj.WaitBarVars.NumCompleted);
                    drawnow;
                end
            end
        end
        
        function deleteWaitBar(obj)
            if ~isempty(obj.WaitBarVars.waitBar)
                % Delete wait bar if one was created
                destroy(obj.WaitBarVars.waitBar)
            end
        end
    end
    
    methods (Access = private)
        
        function numblocksRC = getNumberOfblocks(obj, level)
            refObj = obj.SpatialReferencing(level);
            numblocksRC = ceil(refObj.ImageSize ./ obj.getBlockSize(level));
        end
        
        function checkMaskCompatibility(obj, masks, levels)
            for idx = 1:numel(masks)
                thisMask = masks(idx);
                if thisMask.Channels > 1
                    error(message('images:bigimage:maskChannels'))
                end
                
                thisLevel = levels(idx);
                refObjImage = obj.SpatialReferencing(thisLevel);
                refObjMask = thisMask.SpatialReferencing(thisMask.FinestResolutionLevel);
                if ~isequal(refObjImage.XWorldLimits, refObjMask.XWorldLimits) || ...
                        ~isequal(refObjImage.YWorldLimits, refObjMask.YWorldLimits)
                    error(message('images:bigimage:incompatibleMaskSizes'))
                end
            end
        end
        
        function validateExtraImageLevels(obj, extraLevels, extraImages, level)
            validateattributes(extraLevels, {'numeric'}, ...
                {'integer', 'positive', 'numel', numel(extraImages)}, ...
                'bigimage', 'ExtraImageLevels')
            obj.validateExtraImagesAndLevels(level, extraImages, extraLevels);
        end
        
        function validateExtraImages(obj, extraImages, extraLevels, level)
            if isempty(extraImages)
                return
            end
            
            if isempty(extraLevels)
                extraLevels = repmat(level, [numel(extraImages), 1]);
            end
            obj.validateExtraImagesAndLevels(level, extraImages, extraLevels);
        end
        
        function validateExtraImagesAndLevels(obj, level, extraImages, extraLevels)
            refObj = obj.SpatialReferencing(level);
            for i = 1:numel(extraImages)
                thisExtraLevel = extraLevels(i);
                thisExtraImage = extraImages(i);
                
                if  thisExtraLevel > thisExtraImage.NumLevels
                    error(message('images:bigimage:extraImageLevelsValue'))
                end
                
                exRef = thisExtraImage.SpatialReferencing(thisExtraLevel);
                if ~isequal(refObj.XWorldLimits, exRef.XWorldLimits)...
                        ||~isequal(refObj.YWorldLimits, exRef.YWorldLimits)
                    error(message('images:bigimage:incompatibleExtents'))
                end
            end
        end
    end
end


function tf = validateUseParallel(useParallel)

isLogicalScalar(useParallel);
if useParallel && ~matlab.internal.parallel.isPCTInstalled()
    error(message('images:bigimage:couldNotOpenPool'))
else
    tf = true;
end
end

function tf = valdiateOutputFolder(folder)
if ischar(folder)
    folder = string(folder);
end
validateattributes(folder, {'string', 'char'},...
    {'nonempty','scalar'}, mfilename, 'OutputFolder');
tf = true;
end

function tf = validateClassNames(classNames)

validateattributes(classNames, {'string', 'char', 'cell'}, {'nonempty'}, 'bigimage', 'ClassNames')
if iscell(classNames)
    tf = iscellstr(classNames);
else
    tf = true;
end

end


function tf = validatePadMethod(padMethod)

if isstring(padMethod) || ischar(padMethod)
    validatestring(padMethod, {'replicate', 'symmetric'}, 'bigimage');
else
    validateattributes(padMethod, "numeric", {'scalar', 'real', 'nonsparse'}, 'bigimage');
end

tf = true;

end


function validateClassesAndLabels(classNames, pixelLabelIDs)

if isempty(classNames) || isempty(pixelLabelIDs)
    error(message('images:bigimage:missingClassNamesOrLabelIDs'))
end
end


function tf = isZeroToOneScalar(input)

validateattributes(input, {'numeric'}, {'scalar', 'real', '>=', 0, '<=', 1})
tf = true;

end


function tf = isLogicalScalar(input)

validateattributes(input, {'logical'}, {'scalar'})
tf = true;

end


function tf = isTwoElementPositiveInteger(input)

validateattributes(input, "numeric", {'size', [1 2], 'integer', 'positive'})
tf = true;

end


function adapter = pickAdapter(source, mode)

if isnumeric(source)||islogical(source)||iscategorical(source)
    validateattributes(source, {'numeric', 'logical', 'categorical'},...
        {'nonempty'}, mfilename, 'VAR');
    if ndims(source)>3
        validateattributes(source, {'numeric', 'logical', 'categorical'},...
            {'nonempty', 'ndims', 3}, mfilename, 'VAR');
    end
    adapter = images.internal.adapters.MatrixAdapter('variable', source);
    return
end

if ischar(source)
    source = string(source);
end
validateattributes(source, {'string', 'char'},...
    {'nonempty','scalar'}, mfilename, 'SOURCENAME');

if isfolder(source)
    adapter = images.internal.adapters.BinAdapter(source, mode);
    return
end

source = char(source);
[~,~,ext] = fileparts(source);
ext(ext == '.') = '';

switch lower(ext)
    case {'tif', 'tiff'}
        adapter = images.internal.adapters.TIFFAdapter(source, mode);
    otherwise
        error(message('images:bigimage:unknownImageFormat', source))
        
end

end


function tf = willBeBigImage(x, sizeOfOneblock, borderSize)
% Output will be a bigimage if the output block size is either:
% - Same size as input (which may have borders included if bordersize was
% specified)
% - Or smaller than input with the borders trimmed where the user fcn
% trimmed off the additional border passed in
inputBlockSize = sizeOfOneblock(1:2);
inputBlockSizeWithoutBorders = sizeOfOneblock(1:2) - 2*borderSize;
tf = isequal(size(x,[1,2]), inputBlockSize)...
    || isequal(size(x,[1,2]), inputBlockSizeWithoutBorders);
end


function tf = compatibleSizes(inBlock, outBlock, borderSize)
% Output block size either matches the input block size or, if border size
% was specified, the trimmed intput block size
inBlockSize = size(inBlock,[1,2]);
inBlockSizeWithoutBorder = size(inBlock,[1, 2]) - 2*borderSize;
outBlockSize = [size(outBlock,1), size(outBlock,2)];

tf = isequal(inBlockSize, outBlockSize) ...
    || isequal(inBlockSizeWithoutBorder, outBlockSize);
end


function outblock = trimPadding(outputBlock, inputBlockSize, borderSize)

outblock = outputBlock;
outBlockSize = size(outputBlock,[1,2]);

if ~isequal(borderSize, [0 0]) && isequal(outBlockSize, inputBlockSize)
    % Border size was specified
    % AND the border was NOT trimmed by user callback, so trim it now.    
startIdx = borderSize + 1;
    endIdx = [size(outputBlock,1), size(outputBlock,2)] - borderSize;
    outblock = outputBlock(startIdx(1):endIdx(1), startIdx(2):endIdx(2), :);
end
end


function tf = validateMasks(mask)

validateattributes(mask, {'bigimage'}, {'scalar','nonempty'})
tf = true;

end


function varargout = applyWorkerParallel(workerIdx, numWorkers, fcn, CMainImg, CExtraImg, CMasks, args, pDataQ)

varargout = cell(1, args.NargoutRequested);

levels = args.Levels;

imgLocal = CMainImg.Value;

% bigimage object(s) did not de-serialize correctly, almost always due to
% loadobj erroring out.
if isempty(imgLocal.BlockSize)
    error(message('images:bigimage:parallelLoadFailed'));
end
extraImages = CExtraImg.Value;
for eInd = 1:numel(extraImages)
    if isempty(extraImages(eInd).BlockSize)
        error(message('images:bigimage:parallelLoadFailed'));
    end
end

numMasks = numel(CMasks.Value);
masks = [CMasks.Value];
if isempty(masks)
    masksLocal = [];
end
for i = 1:numMasks
    masksLocal(i) = masks(i); %#ok<AGROW>
end

% Create datastores for the images and split them.
ds = makeDatastoreForApply(imgLocal, levels(1), masksLocal, args);
numPartitions = numWorkers;

dsLocal = ds.partition(numPartitions, workerIdx);
if ~hasdata(dsLocal)
    return
end

% Execute the function, creating new images and results.
updateWaitBarFcn = @()send(pDataQ,1);
% Abort is handled in the client
isWaitBarAbortedFcn = @()false;


[varargout{:}] = applyWorkerCommon(fcn, dsLocal, extraImages, args, updateWaitBarFcn, isWaitBarAbortedFcn);

end


function varargout = applyWorkerSequential(fcn, img, extraImages, masks, args, updateWaitBarFcn, isWaitBarAbortedFcn)

levels = args.Levels;

% Create datastores for each image block contributor.
ds = makeDatastoreForApply(img, levels(1), masks, args);

if ~hasdata(ds)
    error(message('images:bigimage:nothingToApply'))
end

varargout = cell(1, args.NargoutRequested);
[varargout{:}] = applyWorkerCommon(fcn, ds, extraImages, args, updateWaitBarFcn, isWaitBarAbortedFcn);

end


function varargout = applyWorkerCommon(fcn, ds, extraImages, args, updateWaitBarFcn, isWaitBarAbortedFcn)

img = ds.Images;
levels = args.Levels;


% Process function once to inspect the outputs
[requiredblock, blockInfo, blockLocationsIntrinsic] = read(ds);
requiredblock = requiredblock{1};

extrablocks = cell(size(extraImages));
for i = 1:numel(extraImages)
    elevel = levels(i+1);
    extrablocks{i} = extraImages(i).getRegionPadded(elevel, ...
        blockInfo.DataStartWorld, blockInfo.DataEndWorld, args.PadMethod);
end
results = cell(1, args.NargoutRequested);
[results{:}] = fcn(requiredblock, extrablocks{:});
isBigOutput = cellfun(@(x) willBeBigImage(x, size(requiredblock), args.BorderSize), results);


% Initialize the outputs
varargout = cell(1, args.NargoutRequested);
for i = 1:numel(isBigOutput)
    if isBigOutput(i)
        ref = img.SpatialReferencing(levels(1));
        
        % Prepare the adapter
        outFolder = [args.OutputFolder, filesep, args.TimeStamp, '_outputArg_',num2str(i)];
        adapter = images.internal.adapters.BinAdapter(outFolder,'a');
        adapter.writeMetadata(ref.ImageSize, args.BlockSize, size(results{i},3), class(results{i}),[]);
        
        if iscategorical(results{i})
            classnames = categories(results{i});
            varargout{i} = bigimage(adapter,'classnames', classnames, ...
                'PixelLabelIDs', 1:numel(classnames));
            
        else
            varargout{i} = bigimage(adapter);
        end
        varargout{i}.SpatialReferencing = ref;
        varargout{i}.setBlockIntrinsic(1, [blockLocationsIntrinsic{2:3}], trimPadding(results{i},size(requiredblock),args.BorderSize))
    elseif isstruct(results{i})
        results{i}.blockOriginWorld = blockInfo.BlockStartWorld;
        varargout{i} = struct2table(results{i}, 'asarray', true);
    else
        varargout{i} = [results(i), {blockInfo.BlockStartWorld}];
    end
end

% Iteratively execute the function handle on rest of the blocks from the
% datastore.
while hasdata(ds)
    % First, compulsory image
    [requiredblock, blockInfo, blockLocationsIntrinsic] = read(ds);
    requiredblock = requiredblock{1};
    
    % Extra images
    for i = 1:numel(extraImages)
        elevel = levels(i+1);
        extrablocks{i} = extraImages(i).getRegionPadded(elevel, ...
        blockInfo.DataStartWorld, blockInfo.DataEndWorld, args.PadMethod);
    end
    
    [results{:}] = fcn(requiredblock, extrablocks{:});
    
    for i = 1:numel(results)
        if isBigOutput(i)
            if ~compatibleSizes(requiredblock, results{i}, args.BorderSize)
                error(message('images:bigimage:incompatibleResultSizes', ...
                    i, sprintf("%04.4g", blockInfo.BlockStartWorld(1)), ...
                    sprintf("%04.4g", blockInfo.BlockStartWorld(2)), ...
                    size(requiredblock,1), size(requiredblock,2), ...
                    size(results{i},1), size(results{i},2)))
            end
            results{i} = trimPadding(results{i}, size(requiredblock), args.BorderSize);
            varargout{i}.setBlockIntrinsic(1, [blockLocationsIntrinsic{2:3}], results{i})
        elseif isstruct(results{i})
            results{i}.blockOriginWorld = blockInfo.BlockStartWorld;
            varargout{i} = vertcat(varargout{i}, struct2table(results{i}, 'asarray', true));
        else
            varargout{i} = vertcat(varargout{i}, [results(i), {blockInfo.BlockStartWorld}]);
        end
    end
    
    updateWaitBarFcn();
    if isWaitBarAbortedFcn()
        % Aborted via waitbar
        return
    end
end

end

function ds = makeDatastoreForApply(img, level, masks, args)

if isempty(masks)
    masks = bigimage.empty();
end
ds = bigimageDatastore(img, level, ...
    'BlockSize', args.BlockSize, ...
    'PadMethod', args.PadMethod, ...
    'Mask', masks, ...
    'InclusionThreshold', args.InclusionThreshold);
ds.BorderSize = args.BorderSize;

end


function part1 = merge(part1, part2)

if isempty(part2) || isempty(part2{1})
    return
end

for idx = 1:numel(part1)
    item1 = part1{idx};
    item2 = part2{idx};
    
    assert(isempty(item1) || isequal(class(item1), class(item2)))
    
    % All the output bigimages point to the same output folder, so just
    % picking one of them is enough. For the rest, manually merge
    if ~isa(item1, 'bigimage')
        part1{idx} = vertcat(item1, item2);
    end
end
end


function spatialReferencing = createSpatialReferencing(resolutionLevelSizes)

metadata.Height = resolutionLevelSizes(:,1);
metadata.Width  = resolutionLevelSizes(:,2);


resolutionLevelSizes = [metadata.Height(:), metadata.Width(:)];

maxLevel = 1;
maxPixels = metadata.Height(1) * metadata.Width(1);

numLevels = numel(metadata.Height);
for idx = 2:numLevels
    numPixels = metadata.Height(idx) * metadata.Width(idx);
    if numPixels > maxPixels
        maxPixels = numPixels;
        maxLevel = idx;
    end
end

maxSizes = [metadata.Height(maxLevel), metadata.Width(maxLevel)];

XLims = [0.5, maxSizes(2) + 0.5];
YLims = [0.5, maxSizes(1) + 0.5];

numLevels = size(resolutionLevelSizes,1);

for i = numLevels:-1:1  % Iterate backward to preallocate array.
    levelDims = resolutionLevelSizes(i,:);
    spatialReferencing(i) = imref2d(levelDims, XLims, YLims);
end
end

function validatePathsOnWorkers(ppool, fileOrFolders, outputFolder)
if(~ppool.SpmdEnabled)
    error(message('images:bigimage:spmdRequired'));
end

for loc = fileOrFolders
    loc = loc{1}; %#ok<FXSET>
    if strcmp(loc,'variable')
        % in-mem variable, nothing to check
        continue;
    end
    spmd
        isAvailable = isfolder(loc) || isfile(loc);
    end
    if  ~all([isAvailable{:}])
        error(message('images:bigimage:LocationNotAvailableToParallelWorkers', loc));
    end
end

% Verify output folder, use different message
spmd
    isAvailable = isfolder(outputFolder);
end
if  ~all([isAvailable{:}])
    error(message('images:bigimage:OutputNotAvailableToParallelWorkers', outputFolder));
end
end

