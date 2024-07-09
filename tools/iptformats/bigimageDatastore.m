classdef bigimageDatastore < handle & ...
        matlab.mixin.Copyable & ...
        matlab.io.Datastore & ...
        matlab.io.datastore.Partitionable & ...
        matlab.io.datastore.Shuffleable
    %bigimageDatastore Datastore for use with blocks of bigimage objects
    % BIMDS = bigimageDatastore(BIGIMAGES, LEVELS) creates a
    % bigimageDatastore object using a vector of bigimage objects at
    % corresponding resolution levels specified in the numeric array,
    % LEVELS. LEVELS should either be a scalar or a numeric integer valued
    % vector equal to the length of the BIGIMAGES vector.
    %
    % BIMDS = bigimageDatastore(BIGIMAGES, LEVELS, Name, Value) provides
    % additional arguments using one or more name-value pairs. Supported
    % parameters include:
    %
    % 'BlockSize'             Sets the BlockSize property of the
    %                         bigimageDatastore. It is a 1-by-2 ([rows,
    %                         cols]) numeric array specifying the size of
    %                         data to load in each call to read(). Default
    %                         value is the first bigimage's BlockSize
    %                         property.
    %
    % 'BlockOffsets'          Sets the BlockOffsets property of the
    %                         bigimageDatastore. It is a 1-by-2 [rows,
    %                         columns] numeric array specifying the spacing
    %                         between two adjacent blocks. Use this
    %                         property to specify overlapping blocks.
    %                         Default value is the same as the BlockSize
    %                         property resulting in non-overlapping blocks.
    %
    % 'Masks'                 Specifies an array of bigimage objects as the
    %                         masks to be used on BIGIMAGES. It is either a
    %                         scalar or a vector of bigimage masks equal in
    %                         length to BIGIMAGES. Each mask in MASKS
    %                         should have the same spatial extents as the
    %                         corresponding bigimage in BIGIMAGES.
    %
    % 'InclusionThreshold'    Sets the InclusionThreshold property. This
    %                         threshold is as a value in [0, 1] specifying
    %                         the minimum percentage of nonzero/true pixels
    %                         in MASKS that are required to use the
    %                         corresponding block in BIGIMAGES. A value of
    %                         0 implies that at least one mask pixel must
    %                         be nonzero. A value of 1 implies that all
    %                         mask pixels in the region must be nonzero.
    %                         The length of this vector must be the same as
    %                         the length of Masks.
    %
    % 'IncompleteBlocks'      Sets the IncompleteBlocks property. This
    %                         property controls how incomplete blocks are
    %                         handled. METHOD should be a string with one
    %                         of these values:
    %
    %                         'exclude'  - Incomplete blocks are excluded
    %                                      from read().
    %                         'same'     - Incomplete blocks are included,
    %                                      and read() returns the partial
    %                                      data as-is. data will be smaller
    %                                      than specified BlockSize for
    %                                      these incomplete blocks. This is
    %                                      the default value.
    %                         'pad'      - Incomplete blocks are padded
    %                                      based on the PadMethod property.
    %
    % bigimageDatastore properties:
    %   BlockSize          - Block size of a data unit
    %   BlockOffsets       - Space between two blocks
    %   Images             - Source of image blocks, array of bigimage objects
    %   IncompleteBlocks   - Controls incomplete edge blocks
    %   InclusionThreshold - Mask threshold value controlling block inclusion
    %   Levels             - Resolution levels used to read blocks
    %   Masks              - Masks, array of bigimage objects
    %   PadMethod          - Method for padding incomplete blocks
    %   BorderSize         - Padding size required around a block
    %   ReadSize           - Number of blocks to read per read() call
    %
    %  bigimageDatastore methods:
    %    hasdata       - Returns true if there is more data in the datastore
    %    read          - Read the next block
    %    readRelative  - Read a neighboring block with a relative position
    %    reset         - Resets the datastore to the start of the data
    %    partition     - Returns a new datastore that represents a single
    %                    partitioned portion of the original datastore
    %    preview       - Reads the first block
    %    numpartitions - Returns an estimate for a reasonable number of
    %                    partitions to use with the partition function,
    %                    according to the total data size
    %    shuffle       - Shuffles the files of ImageDatastore using randperm
    %    transform     - Create an altered form of the current
    %                    bigimageDatastore by specifying a function handle that
    %                    will execute after read on the current
    %                    bigimageDatastore.
    %    combine       - Create a new bigimageDatastore that horizontally
    %                    concatenates the result of read from two or more input
    %                    bigimageDatastore.
    %
    %   Example 1
    %   ---------
    %   % Create a bigimageDatastore at a specific level and block size
    %     bim = bigimage('tumor_091R.tif');
    %     subplot(2,1,1)
    %     bigimageshow(bim);
    %     subplot(2,1,2)
    %     % Create a datastore at level 2, with blocksize of [512 512]
    %     bimds = bigimageDatastore(bim, 2, 'BlockSize', [512 512]);
    %     % Read four blocks at a time
    %     bimds.ReadSize = 4;
    %     while hasdata(bimds)
    %         blocks = read(bimds);
    %         % Blocks are returned as cell arrays. Partial edge blocks 
    %         % will have smaller size than interior blocks.
    %         disp(blocks);
    %         % Display the blocks
    %         montage(blocks,'Size', [1 4], 'BorderSize',5,'BackgroundColor','b');
    %         title('Press any key to continue');
    %         pause;
    %     end
    %     title('');
    %
    %   Example 2
    %   ---------
    %   % Create a bigimageDatastore with overlapping blocks
    %     bim = bigimage('tumor_091R.tif');
    %     bimds = bigimageDatastore(bim, 3, ...
    %                'BlockSize', [300 300], 'BlockOffsets', [100 100],...
    %                'IncompleteBlocks','exclude');
    %     % Read two blocks at a time
    %     bimds.ReadSize = 2;
    %     while hasdata(bimds)
    %         blocks = read(bimds);
    %         disp(blocks);
    %         % Display the blocks
    %         montage(blocks,'BorderSize',5,'BackgroundColor'   ,'b');
    %         title('Press any key to continue');
    %         pause;
    %     end
    %     title('');
    %
    %   Example 3
    %   ---------
    %   % Create a bigimageDatastore using a coarse level mask
    %     bim = bigimage('tumor_091R.tif');
    %
    %     % Create a mask at the coarsest level
    %     clevel = bim.CoarsestResolutionLevel;
    %     imcoarse = getFullLevel(bim, clevel);
    %     stainMask = ~imbinarize(rgb2gray(imcoarse));
    %     % Retain the original spatial referencing information
    %     bmask = bigimage(stainMask,...
    %                'SpatialReferencing', bim.SpatialReferencing(clevel));    
    %     figure
    %     bigimageshow(bmask); 
    %
    %     % Create a bigimagedatastore for blocks which are fully within
    %     % the stained region as defined by the mask.
    %     bimds = bigimageDatastore(bim, 1, 'BlockSize', [256 256],...
    %                'Masks', bmask, 'InclusionThreshold',1);
    %
    %     bimds.ReadSize = 4;
    %     figure
    %     while hasdata(bimds)
    %         blocks = read(bimds);
    %         disp(blocks);
    %         % Display the blocks
    %         montage(blocks,'BorderSize',5,'BackgroundColor'   ,'b');
    %         title('Press any key to continue');
    %         pause;
    %     end
    %     title('');
    %
    % See also bigimage, imageDatastore
    
    
    % Copyright 2018-2019 The MathWorks, Inc.
    
    properties (Access = public)
        % ReadSize Number of blocks to read per read() call
        %   Default value is 1.
        ReadSize = 1
        
        % BorderSize  Padding size required around a block
        %   Default value is [0 0].
        BorderSize = [0 0]
        
        % PadMethod  Method for padding incomplete blocks
        %   PadMethod is either a numeric scalar or one of the following
        %   strings:
        %     'replicate' - Repeats the border elements
        %     'symmetric' - Pads with mirror reflections
        %   Default value is 0.
        PadMethod
        
        % InclusionThreshold Mask threshold value controlling block inclusion
        %   Default value is 0.5.  Changing this value resets the datastore.
        InclusionThreshold
    end
    
    properties (SetAccess = private)
        % Masks Masks, array of bigimage objects
        %   This is a read-only property
        Masks = bigimage.empty()
        
        % Images Source of image blocks, array of bigimage objects
        %   This is a read-only property
        Images
        
        % Levels Resolution levels used to read blocks
        %   This is a read-only property
        Levels
        
        % BlockSize Block size of a data unit
        %   Default value is set from the first level BlockSize property of
        %   the first bigimage object in Images property. This is a
        %   read-only property
        BlockSize
        
        % BlockOffsets Space between two blocks
        %   Default value is equal to the BlockSize property. Use this
        %   property to obtain overlapping blocks. For example, specifying
        %   a value of [1 1] will result in blocks that slide by 1 pixel
        %   across the rows, and then 1 pixel down the columns and so on.
        %   This is a read-only property
        BlockOffsets
        
        % IncompleteBlocks Controls incomplete edge blocks
        %   Default value is 'same'
        %   This is a read-only property
        IncompleteBlocks = 'same'
    end
    
    properties (Hidden = true, SetAccess = private, Dependent = true)
        Length
    end
    
    properties (Access = private)
        CurrentImageNumber = 1
        CurrentRowIntrinsic = 1
        CurrentColIntrinsic = 1
        CurrentImage
        CurrentLevel
        CurrentImageSize
        
        NumberOfPartitions = 1;
        PartitionNumber = 1;
        
        %%% Even newer Implementation
        BlockLUT
        ReadOrder
        CurrentReadIndex
        
        AtFirstBlock = false
    end
    
    methods (Access = public)
        function obj = bigimageDatastore(bigimages_, levels_, varargin)  % v3.1
            validateattributes(bigimages_, "bigimage", {'nonempty', 'vector'}, mfilename, "bigimages", 1)
            obj.Images = bigimages_;
            
            if isscalar(levels_)
                % make it equal to number of images
                levels_ = repmat(levels_, [1 numel(obj.Images)]);
            else
                % Vector of levels
                if isscalar(obj.Images)
                    % replicate single image to match number of levels
                    obj.Images = repmat(obj.Images, [1 numel(levels_)]);
                end
            end
            
            % At this point, numel(images)==numel(levels)
            numImages = numel(obj.Images);
            validateattributes(levels_, "numeric", {"integer","positive", "vector", "numel", numImages}, mfilename, "levels", 2)
            
            % Each level should be valid for its corresponding bigimage
            for ind = 1:numel(levels_)
                numLevels = numel(obj.Images(ind).SpatialReferencing);
                validateattributes(levels_(ind), "numeric", {'<=', numLevels}, mfilename, "levels", 2);
            end
            
            obj.Levels = double(levels_);
            
            if nargin > 2
                parser = inputParser;
                parser.FunctionName = 'bigimage';
                parser.CaseSensitive = false;
                parser.PartialMatching = true;
                parser.KeepUnmatched = false;
                parser.addParameter('BlockSize', obj.findDefaultBlockSize(), @(blockSize)...
                    validateattributes(blockSize, "numeric", ...
                    {"positive", "integer", "row", "numel", 2}, mfilename, "BlockSize"));
                parser.addParameter('BlockOffsets', [], @(blockOffsets)...
                    validateattributes(blockOffsets, "numeric", ...
                    {"positive", "integer", "row", "numel", 2}, mfilename, "BlockOffsets"));
                parser.addParameter('PadMethod', 0, @validatePadMethod)
                % Allow using []
                parser.addParameter('Masks', bigimage.empty(), @(x) ...
                    validateattributes(x, {'bigimage'}, ...
                    {}, mfilename, "Masks"));
                parser.addParameter('InclusionThreshold', 0.5, ...
                    @(icth) validateInclusionThreshold(icth, numel(obj.Images)))
                parser.addParameter('IncompleteBlocks', 'same');
                parser.addParameter('ReadSize', []);
                parser.parse(varargin{:});
                obj.BlockSize = parser.Results.BlockSize;
                obj.PadMethod = parser.Results.PadMethod;
                
                obj.Masks = parser.Results.Masks;
                if ~isempty(obj.Masks)
                    if isscalar(obj.Masks)
                        obj.Masks = repmat(obj.Masks, [1 numImages]);
                    else
                        validateattributes(obj.Masks, "bigimage", {'numel', numImages}, mfilename, "Mask");
                    end
                    obj.validateMaskExtents();
                end
                % Turn the memory cache on for masks, and validate channels
                for ind = 1:numel(obj.Masks)
                    obj.Masks(ind).Adapter.UseMemoryCache = true;
                    if obj.Masks(ind).Channels > 1
                        error(message('images:bigimage:maskChannels'))
                    end
                end
                
                obj.InclusionThreshold = parser.Results.InclusionThreshold;
                
                str = validatestring(parser.Results.IncompleteBlocks,...
                    {'exclude', 'same', 'pad'}, mfilename, "IncompleteBlocks");
                obj.IncompleteBlocks = str;
                
                if ~contains('ReadSize', parser.UsingDefaults)
                    obj.ReadSize = parser.Results.ReadSize;
                end
            else
                obj.BlockSize = obj.findDefaultBlockSize();
                obj.PadMethod = 0;
            end
            
            if nargin > 2 && ~isempty(parser.Results.BlockOffsets)
                obj.BlockOffsets = parser.Results.BlockOffsets;
            else
                % Default value
                obj.BlockOffsets = obj.BlockSize;
            end
            if any(obj.BlockOffsets < obj.BlockSize)
                % Overlapping blocks, turn on cache
                for ind=1:numel(obj.Images)
                    obj.Images(ind).Adapter.UseMemoryCache = true;
                end
            end
            
            obj.buildLUT()
            
            obj.reset()
        end
        
        function tf = hasdata(obj)
            %hasdata Returns true if more data is available.
            %     TF = hasdata(bimds) returns a logical scalar TF
            %     indicating availability of data. This method should be
            %     called before calling read. hasdata is used in
            %     conjunction with read to read all the data within the
            %     bigimageDatastore.
            %
            %     Example
            %     -------
            %     bim = bigimage('tumor_091R.tif');
            %     bimds = bigimageDatastore(bim,1);
            %     while hasdata(bimds)
            %         [data, info] = read(bimds);
            %     end
            %
            % See also bigimageDatastore, read, reset
            
            % There is more data if...
            %   * (1) within the current image the are more blocks
            %   * ... OR (2) there are more images in the datastore
            %   * AND (3) the masks for the remaining blocks yield true
            %
            
            if obj.AtFirstBlock
                obj.AtFirstBlock = false;
                if ~isempty(obj.Masks) && ~obj.currentBlockSatisfiesMask()
                    % If masks were given and the first block does not
                    % satisfy the mask, keep moving till we get one.
                    obj.switchToNextBlock()
                end
            end
            
            tf = ~obj.readIndexPassedEnd();
        end
        
        function subds = partition(obj, numP, idx)
            % partition Return a partitioned part of the bigimageDatastore.
            %     SUBBIMDS = partition(BIMDS,N,INDEX) partitions DS into N
            %     parts and returns the partitioned Datastore, SUBDS,
            %     corresponding to INDEX. An estimate for a reasonable
            %     value for N can be obtained by using the NUMPARTITIONS
            %     function.
            %
            %     Example
            %     -------
            %     bim = bigimage('tumor_091R.tif');
            %     bimds = bigimageDatastore(bim,2);
            %
            %     bimdsp1 = partition(bimds, 2, 1);
            %     disp('Partition 1');
            %     while hasdata(bimdsp1)
            %         [data, info] = read(bimdsp1);
            %         disp(info);
            %     end           
            %
            %     bimdsp2 = partition(bimds, 2, 2);
            %     disp('Partition 2');
            %     while hasdata(bimdsp2)
            %         [data, info] = read(bimdsp2);
            %         disp(info);
            %     end           
            %
            % See also bigimageDatastore, numpartitions, maxpartitions
            
            % In partitioned bigimageDatastores, the N datastores start at
            % neighboring locations. (i.e., P1 starts at block 1, P2 at 2,
            % etc.) On each read() the datastore's read index skips forward
            % N blocks. (When shuffled, it skips ahead N random read
            % positions.) Incomplete read behavior and masked out blocks is
            % taken into account when moving to the first or next block
            % after reset() or read() is called.
            %
            % This is an implementation detail and should not be counted on
            % between releases.
            
            validateattributes(numP, "numeric", {"scalar", "nonempty", "positive", "integer"}, ...
                "partition", "N");
            validateattributes(idx, "numeric", {"scalar", "nonempty", "positive", "integer", "<=", numP}, ...
                "partition", "INDEX");
            subds = copy(obj);
            
            subds.NumberOfPartitions = numP;
            subds.PartitionNumber = idx;
            subds.reset()
        end
        
        function oneblock = preview(obj)
            %  preview   Reads the first block
            %     b = preview(bimds) Returns the first block from the start
            %     of the datastore.
            %
            %     See also bigimageDatastore, read, hasdata, reset, readall,
            %     progress.
            oneblock = obj.readOneBlock([1 1], 1, obj.BlockLUT.Levels(1));
        end
        
        function amount = progress(obj)
            %  progress  Amount of datastore read.
            %     R = progress(BIMDS) gives the ratio of the datastore
            %     BIMDS that has been read.
            
            possibleReads = numel(obj.ReadOrder);
            amount = (obj.CurrentReadIndex - 1) / possibleReads;
        end
        
        function [multiReadData, infoStruct, varargout] = read(obj)
            % read Read data and information about the extracted data.
            %     BCELL = read(BIMDS) Returns the data extracted from the
            %     bigimageDatastore, BIMDS. BCELL is a cell array of block
            %     data of length ReadSize.
            %
            %     [B, INFO] = read(BIMDS) also returns information about
            %     where the data was extracted from the bigimageDatastore.
            %     INFO is a scalar struct with the following fields. These
            %     fields are arrays if ReadSize>1.
            %
            %       Level           - The level from which this data was
            %                         read.
            %       ImageNumber     - An index into the bimds.Images
            %                         array corresponding to the bigimage
            %                         from which this block was read.
            %       BlockStartWorld - The center world coordinates of the
            %                         top left pixel of the block,
            %                         excluding any padding.
            %       BlockEndWorld   - The center world coordinates of the
            %                         bottom right pixel of the block,
            %                         excluding any padding.
            %       DataStartWorld  - The center world coordinates of the
            %                         top left pixel of the block,
            %                         including padding pixels.
            %       DataEndWorld    - The center world coordinates of the
            %                         bottom right pixel of the block,
            %                         including padding pixels.
            %
            % Note: When a bigimageDatastore is created with a Mask, the
            % read method includes computation time needed to identify a
            % valid image block which satisfies the Mask at the specified
            % InclusionThreshold. This will result in varying run times
            % which depend on the Mask size and sparsity.
            %            
            %     Example
            %     -------
            %     bim = bigimage('tumor_091R.tif');
            %     bimds = bigimageDatastore(bim,1);
            %     while hasdata(bimds)
            %         [data, info] = read(bimds);
            %         disp(info);
            %     end
            %
            %  See also bigimageDatastore, readsize
            
            multiReadData = cell(obj.ReadSize, 1);
            multiReadInfo = {};
            
            for idx = 1:obj.ReadSize
                if ~obj.hasdata()
                    if idx==1 % Nothing could be read
                        error(message('images:bigimage:noMoreData'));
                    end
                    multiReadData(idx:end) = [];
                    break
                end
                
                regionStartIntrinsic = [obj.CurrentRowIntrinsic obj.CurrentColIntrinsic];
                [data, info] = obj.readOneBlock(regionStartIntrinsic, obj.CurrentImageNumber, obj.CurrentLevel);
                
                varargout{1} = {obj.CurrentLevel, obj.CurrentRowIntrinsic, obj.CurrentColIntrinsic};
                
                % Point to the next valid block (Ensures hasdata() works
                % correctly too)
                obj.switchToNextBlock()
                
                multiReadData{idx} = data;
                multiReadInfo{idx} = info; %#ok<AGROW>
            end
            
            % Flip the info struct inside out to make it easier to use
            structArray = [multiReadInfo{:}];
            infoStruct.Level = [structArray.Level];
            infoStruct.ImageNumber = [structArray.ImageNumber];
            infoStruct.BlockStartWorld = reshape([structArray.BlockStartWorld],[], numel(multiReadInfo))';
            infoStruct.BlockEndWorld = reshape([structArray.BlockEndWorld],[], numel(multiReadInfo))';
            infoStruct.DataStartWorld = reshape([structArray.DataStartWorld],[], numel(multiReadInfo))';
            infoStruct.DataEndWorld = reshape([structArray.DataEndWorld],[], numel(multiReadInfo))';
        end
        
        function data = readall(obj) %#ok<MANU,STOUT>
            error(message('images:bigimage:readallNotSupported'))
        end
        
        function [data, info] = readRelative(obj, infoStruct, blockOffsets)
            % readRelative Read neighboring block using relative position
            %     b = readRelative(bimds, sinfo, boffset) reads a
            %     neighboring block. sinfo is a struct (as returned by read
            %     or readRelative) that specifies the source block. boffset
            %     is a 1-by-2 integer valued vector specifying the offset
            %     from the source block. Offset is specified in units of
            %     blocks. b is [] when boffset puts the neighboring block
            %     out of bounds of the corresponding image.
            %
            %     [b, rinfo] = readRelative(bimds, sinfo, boffset) also
            %     returns the info struct of the block that was read.
            %
            %     Input info structs need to have the following fields:
            %       Level           - The level from which this data was
            %                         read.
            %       ImageNumber     - An index into the input BIGIMAGES
            %                         array corresponding to the bigimage
            %                         from which this block was read.
            %       BlockStartWorld - The coordinates of the extreme top
            %                         left of the block.
            %
            %     Returned info structs have the same format as those
            %     returned by the read() function.
            %
            %     NOTE: Masks are not used; the requested blocks are always
            %     read. The Incompleteblock behavior is respected. If the
            %     requested block is incomplete and IncompleteBlocks has a
            %     value of 'exclude', an empty block will be returned.
            %     PadMethod and BorderSize are also honored.
            %
            %     Example
            %     -------
            %     % Read the 4-connected neighbor blocks
            %     bim = bigimage('cameraman.tif');
            %     bimds = bigimageDatastore(bim,1, 'BlockSize', [64 64]);
            %     % Read the first block
            %     [b, sinfo] = read(bimds);
            %     b = b{1};
            %     % Read its four neighbors:
            %     bLeft   = readRelative(bimds, sinfo, [0 -1]);
            %     bTop    = readRelative(bimds, sinfo, [-1 0]);
            %     bRight  = readRelative(bimds, sinfo, [0 1]);
            %     bBottom = readRelative(bimds, sinfo, [1 0]);
            %     % Assemble as a montage to view relative locations
            %     montage({[], bTop, [], bLeft, b, bRight, [], bBottom, []}, ...
            %       'Size', [3 3], 'BorderSize', 5, 'BackgroundColor', 'b')
            %
            %  See also bigimageDatastore, read
            
            validateattributes(infoStruct, "struct", {"scalar", "nonempty"}, "readRelative", "infoStruct") %#ok<*CLARRSTR>
            validateattributes(blockOffsets, "numeric", {"integer", "numel", 2}, "readRelative", "blockOffsets")
            if ~isfield(infoStruct, 'Level') || ~isfield(infoStruct, 'ImageNumber') || ...
                    ~isfield(infoStruct, 'BlockStartWorld') 
                error(message('images:bigimage:missingInfoFields'))
            end
            
            % Validate the info struct
            validateattributes(infoStruct.ImageNumber, "numeric",...
                {"positive", "scalar", "<=", numel(obj.Images)}, ...
                "readRelative", "sinfo.ImageNumber");
            validateattributes(infoStruct.Level, "numeric",...
                {"positive","integer", "scalar", "<=", numel(obj.Images(infoStruct.ImageNumber).SpatialReferencing)},...
                "readRelative", "sinfo.Level");
            validateattributes(infoStruct.BlockStartWorld, "numeric",...
                {"numel", 2}, "readRelative", "sinfo.BlockStartWorld");
            
            % Additional validation for level, it has to correspond to the
            % one the data store was created with
            if ~isequal(infoStruct.Level, obj.Levels(infoStruct.ImageNumber))
                error(message('images:bigimage:readRelativeLevelDoestMatch',num2str(infoStruct.Level), num2str(obj.Levels(infoStruct.ImageNumber))));
            end
            
            % Find intrinsic location of the block to read.
            theImage = obj.Images(infoStruct.ImageNumber);
            ref = theImage.SpatialReferencing(infoStruct.Level);
            bsize = theImage.getBlockSize(infoStruct.Level);
            yxWorld = infoStruct.BlockStartWorld;
            [xInt, yInt] = ref.worldToIntrinsic(yxWorld(1), yxWorld(2));
            yxIntrinsic = round([yInt xInt]);
            
            % Validate that this belong to a block start
            if ~any(all(obj.BlockLUT.StartYX == yxIntrinsic,2))
                error(message('images:bigimage:invalidBlockStart'));
            end
            
            % Compute the start point for output block
            startIntrinsic = yxIntrinsic + obj.BlockOffsets .* blockOffsets(:)';
            
            % Is that out of bounds?
            if any(startIntrinsic > ref.ImageSize) || any((startIntrinsic + bsize - 1) < 1)
                info = struct.empty();
                data = [];
                return
            end
            
            % Or incomplete?
            switch obj.IncompleteBlocks
                case {'exclude', "exclude"}
                    if obj.isIncomplete(startIntrinsic)
                        info = struct.empty();
                        data = [];
                        return
                    end
            end
            
            % Read the neighbor block
            [data, info] = obj.readOneBlock(startIntrinsic, infoStruct.ImageNumber, infoStruct.Level);
        end
        
        function reset(obj)
            % reset  Set next read to begin at start.
            %     reset(bimds) Resets the bigimageDatastore to the state
            %     where no data has been read from it.
            %
            %     Example
            %     -------
            %     bim = bigimage('tumor_091R.tif');
            %     bimds = bigimageDatastore(bim,2);
            %     while hasdata(bimds)
            %         [data, info] = read(bimds);
            %         disp(info);
            %     end
            %
            %     % Reset to read the blocks again
            %     disp('After reset');
            %     reset(bimds);
            %     while hasdata(bimds)
            %         [data, info] = read(bimds);
            %         disp(info);
            %     end           
            %
            %     See also bigimageDatastore, read, hasdata, progress
            
            obj.CurrentReadIndex = 1;
            blockNum = obj.ReadOrder(obj.CurrentReadIndex);
            
            obj.CurrentImageNumber = obj.BlockLUT.ImageNumbers(blockNum);
            obj.CurrentRowIntrinsic = obj.BlockLUT.StartYX(blockNum, 1);
            obj.CurrentColIntrinsic = obj.BlockLUT.StartYX(blockNum, 2);
            obj.CurrentImage = obj.Images(obj.CurrentImageNumber);
            obj.CurrentLevel = obj.Levels(obj.CurrentImageNumber);
            obj.CurrentImageSize = obj.CurrentImage.SpatialReferencing(obj.CurrentLevel).ImageSize;
            
            % Position at first block.
            for ind = 1:(obj.PartitionNumber - 1)
                % For partitioned datastores, skip the required number of
                % blocks.
                obj.moveToNextBlockIndices();
                if obj.readIndexPassedEnd()
                    break
                end
            end
            % Indicate to hasdata() that we are positioned on a valid block.
            obj.AtFirstBlock = true;
        end
        
        function newds = shuffle(obj)
            % shuffle  Permute order blocks will be read.
            %     NEWDS = shuffle(BIMDS) randomly reorders the read order
            %     of the blocks in BIMDS and returns a new
            %     bigimageDatastore NEWDS. The original datastore is
            %     unchanged.
            %
            %     Example
            %     -------
            %     bim = bigimage('tumor_091R.tif');
            %     bimds = bigimageDatastore(bim,2);
            %     while hasdata(bimds)
            %         [data, info] = read(bimds);
            %         disp(info);
            %     end
            %     
            %     sbimds = shuffle(bimds);
            %     disp('Shuffled Order');
            %     while hasdata(sbimds)
            %         [data, info] = read(sbimds);
            %         disp(info);
            %     end
            %     
            %     See also bigimageDatastore, partition
            newds = copy(obj);
            
            newds.ReadOrder = newds.ReadOrder(randperm(numel(newds.ReadOrder)));
            newds.reset()
        end
    end
    
    % Property mutators and accessors
    methods
        function set.BorderSize(obj, newValue)
            validateattributes(newValue, {'numeric'}, {'numel', 2, 'vector', 'nrows', 1, 'real', 'nonnegative', 'integer'}, 'bigimageDatastore', 'BorderSize')
            obj.BorderSize = double(newValue);
        end
        
        function set.ReadSize(obj, newValue)
            validateattributes(newValue, {'numeric'}, {'scalar', 'real', 'positive', 'integer'}, 'bigimageDatastore', 'ReadSize')
            obj.ReadSize = double(newValue);
        end
        
        function set.PadMethod(obj, newValue)
            [~, newValue] = validatePadMethod(newValue);
            obj.PadMethod = newValue;
        end
        
        function set.InclusionThreshold(obj, newValue)
            [~, newValue] = validateInclusionThreshold(newValue, numel(obj.Images)); %#ok<MCSUP>
            obj.InclusionThreshold = newValue;
        end
        
        function L = get.Length(obj)
            L = numel(obj.ReadOrder);  % Assumes unpartitioned.
        end
    end
    
    methods (Access = protected)
        function num = maxpartitions(obj)
            %maxpartitions  Maximum number of possible partitions.
            num = numel(obj.ReadOrder)/obj.NumberOfPartitions;
            if obj.NumberOfPartitions == 1
                % The first one gets any left overs
                num = ceil(num);
            else
                num = floor(num);
            end
        end
    end
    
    methods (Access = private)
        function buildLUT(obj)
            startIntrinsic = zeros(0,2);
            imageNumbers = zeros(0,1);
            levels = zeros(0,1);
            
            for imgIdx = 1:numel(obj.Images)
                thisImage = obj.Images(imgIdx);
                thisLevel = obj.Levels(imgIdx);
                ref = thisImage.SpatialReferencing(thisLevel);
                
                X = 1:obj.BlockOffsets(2):ref.ImageSize(2);
                Y = 1:obj.BlockOffsets(1):ref.ImageSize(1);
                
                [allY, allX] = meshgrid(Y,X);
                thisStartIntrinsic = [allY(:), allX(:)];
                
                startIntrinsic = cat(1, startIntrinsic, thisStartIntrinsic);
                imageNumbers = cat(1, imageNumbers, repmat(imgIdx, [numel(allY), 1]));
                levels = cat(1, levels, repmat(thisLevel, [numel(allY), 1]));
            end
            
            LUT.StartYX = startIntrinsic;
            LUT.ImageNumbers = imageNumbers;
            LUT.Levels = levels;
            obj.BlockLUT = LUT;
            
            obj.ReadOrder = 1:size(startIntrinsic, 1);
        end
        
        function size = findDefaultBlockSize(obj)
            size = obj.Images(1).BlockSize(1,:);
        end
        
        function moveToNextBlockIndices(obj)
            previousImageNumber = obj.CurrentImageNumber;
            previousLevel = obj.CurrentLevel;
            
            obj.CurrentReadIndex = obj.CurrentReadIndex + 1;
            if obj.readIndexPassedEnd()
                % At this point all of the Current* values are invalid.
                % This is okay, since hasdata() is false and read() will
                % return empty.
                return
            end
            blockNum = obj.ReadOrder(obj.CurrentReadIndex);
            
            obj.CurrentRowIntrinsic = obj.BlockLUT.StartYX(blockNum, 1);
            obj.CurrentColIntrinsic = obj.BlockLUT.StartYX(blockNum, 2);
            obj.CurrentImageNumber = obj.BlockLUT.ImageNumbers(blockNum);
            obj.CurrentLevel = obj.BlockLUT.Levels(blockNum);
            
            if previousImageNumber ~= obj.CurrentImageNumber || ...
                    previousLevel ~= obj.CurrentLevel
                
                obj.CurrentImage = obj.Images(obj.CurrentImageNumber);
                obj.CurrentImageSize = obj.CurrentImage.SpatialReferencing(obj.CurrentLevel).ImageSize;
            end
        end
        
        function switchToNextBlock(obj)
            while hasdata(obj)
                % Skip the right number of blocks before looking at whether
                % to use the block. When partitioned, full and incomplete
                % blocks are treated identically while moving to the next
                % block to interrogate.
                for i = 1:obj.NumberOfPartitions
                    obj.moveToNextBlockIndices();
                    passedImage = obj.readIndexPassedEnd();
                    if passedImage
                        break
                    end
                end
                
                % Keep going until a stopping condition is met.
                if passedImage
                    break
                elseif isequal(obj.IncompleteBlocks, 'exclude') && obj.isIncomplete()
                    continue
                elseif isempty(obj.Masks) || obj.currentBlockSatisfiesMask()
                    break
                end
            end
        end
        
        function [data, info] = readOneBlock(obj, regionStartIntrinsic, imageNumber, level)
            info.Level = level;
            info.ImageNumber = imageNumber;
            
            theImage = obj.Images(imageNumber);
            ref = theImage.SpatialReferencing(level);
            if imageNumber == obj.CurrentImageNumber
                imageSize = obj.CurrentImageSize;
            else
                imageSize = ref.ImageSize;
            end
            
            if any(regionStartIntrinsic<1)
                % Out of bounds
                data = []; info = [];
                return
            end
            
            [x,y] = ref.intrinsicToWorld(regionStartIntrinsic(2), regionStartIntrinsic(1));
            info.BlockStartWorld = [y, x];
            
            regionEndIntrinsic = regionStartIntrinsic + obj.BlockSize - 1;
            [x,y] = ref.intrinsicToWorld(regionEndIntrinsic(2), regionEndIntrinsic(1));
            % Clamp
            x = min(x, ref.ImageExtentInWorldX);
            y = min(y, ref.ImageExtentInWorldY);
            info.BlockEndWorld = [y, x];
            
            % Border size
            regionStartIntrinsic = regionStartIntrinsic - obj.BorderSize;
            regionEndIntrinsic = regionEndIntrinsic + obj.BorderSize;
            
            % Start has potentially changed due to bordersize
            [x,y] = ref.intrinsicToWorld(regionStartIntrinsic(2), regionStartIntrinsic(1));
            info.DataStartWorld = [y, x];
            
            % Padding
            paddingNorthAndWest = max(-regionStartIntrinsic + 1, 0);
            paddingSouthAndEast = max(regionEndIntrinsic - obj.CurrentImageSize, 0);
            if strcmp(obj.IncompleteBlocks, 'same')
                paddingNorthAndWest = min(paddingNorthAndWest, obj.BorderSize);
                paddingSouthAndEast = min(paddingSouthAndEast, obj.BorderSize);
            end
            
            % Clamp and read
            regionStartIntrinsic = max(regionStartIntrinsic, [1 1]);
            regionEndIntrinsic = min(regionEndIntrinsic, imageSize);
            data = theImage.getRegionIntrinsic(level, ...
                regionStartIntrinsic, regionEndIntrinsic);
            
            % Pad, if needed
            if any(paddingNorthAndWest | paddingSouthAndEast)
                data = padData(data, paddingNorthAndWest, paddingSouthAndEast, obj.PadMethod);
            end
            
            % Expand data end to account for padding on lower right
            if any(paddingSouthAndEast)
                regionEndIntrinsic = regionEndIntrinsic + paddingSouthAndEast;
            end
            
            % End could have changed to account for border
            [x,y] = ref.intrinsicToWorld(regionEndIntrinsic(2), regionEndIntrinsic(1));
            info.DataEndWorld = [y, x];
            
            % Convert to x,y format
            info.BlockStartWorld = [info.BlockStartWorld(2), info.BlockStartWorld(1)];
            info.BlockEndWorld = [info.BlockEndWorld(2), info.BlockEndWorld(1)];
            info.DataStartWorld = [info.DataStartWorld(2), info.DataStartWorld(1)];
            info.DataEndWorld = [info.DataEndWorld(2), info.DataEndWorld(1)];
        end
               
        function validateMaskExtents(obj)
            % Check that spatial extents should match
            for ind = 1:numel(obj.Images)
                % numel(images) == numel(masks) now
                mask = obj.Masks(ind);                
                mref = mask.SpatialReferencing(mask.FinestResolutionLevel);
                imref = obj.Images(ind).SpatialReferencing(obj.Levels(ind));
                if ~isequal(mref.XWorldLimits, imref.XWorldLimits) ||...
                        ~isequal(mref.YWorldLimits, imref.YWorldLimits)
                    error(message('images:bigimageDatastore:maskImageDoesNotMatch'))
                end
            end
        end
        
        function tf = currentBlockSatisfiesMask(obj)
            regionStartIntrinsic = [obj.CurrentRowIntrinsic obj.CurrentColIntrinsic];
            regionStartIntrinsic = regionStartIntrinsic - obj.BorderSize;
            regionEndIntrinsic = regionStartIntrinsic + obj.BlockSize + 2*obj.BorderSize - 1;
            
            % Convert image intrinsic to world
            imageRef = obj.Images(obj.CurrentImageNumber).SpatialReferencing(obj.CurrentLevel);
            [x,y] = imageRef.intrinsicToWorldAlgo(regionStartIntrinsic(2), regionStartIntrinsic(1));
            maskStartWorld = [x, y];
            [x,y] = imageRef.intrinsicToWorldAlgo(regionEndIntrinsic(2), regionEndIntrinsic(1));
            maskEndWorld = [x, y];
            
            theMask = obj.Masks(obj.CurrentImageNumber);
            pctNNZ = theMask.computeWorldRegionNNZ(theMask.FinestResolutionLevel,...
                maskStartWorld, maskEndWorld);
            
            if pctNNZ== 0
                tf = false;
            else
                tf = pctNNZ>=obj.InclusionThreshold(obj.CurrentImageNumber);
            end
        end
               
        function tf = isIncomplete(obj, varargin)
            assert(nargin <= 2)
            
            if nargin == 1
                startLocation = [obj.CurrentRowIntrinsic obj.CurrentColIntrinsic];
            else
                startLocation = varargin{1};
            end
            expectedEndOfblock = startLocation + obj.BlockSize - 1;
            
            tf = any(expectedEndOfblock > obj.CurrentImageSize);
        end
        
        function tf = readIndexPassedEnd(obj)
            tf = obj.CurrentReadIndex > numel(obj.ReadOrder);
        end
    end
end


function [tf, padMethod] = validatePadMethod(padMethod)

if isstring(padMethod) || ischar(padMethod)
    padMethod = validatestring(padMethod, {'replicate', 'symmetric'}, 'bigimageDatastore');
else
    validateattributes(padMethod, ["numeric", "logical"], {"scalar", "real", "nonsparse"}, 'bigimageDatastore');
end

tf = true;

end


function [tf, icth] = validateInclusionThreshold(icth, numImages)
tf = true;

if numImages==0
    % numImages == 0 on loadobj
    return
end
if isscalar(icth)
    icth = repmat(icth, [1, numImages]);
end
validateattributes(icth, "numeric", {'real', '<=', 1, '>=', 0, 'numel', numImages}, mfilename, "InclusionThreshold")
end


function data = padData(data, paddingNorthAndWest, paddingSouthAndEast, padMethod)

if iscategorical(padMethod)
    % padarray only works with numeric data
    padMethod = double(padMethod);
end

if iscategorical(data)
    convertToCategorical = true;
    origCategories = categories(data);
    data = double(data);
else
    convertToCategorical = false;
end

if paddingNorthAndWest(1)
    data = padarray(data, [paddingNorthAndWest(1) 0], padMethod, 'pre');
end
if paddingNorthAndWest(2)
    data = padarray(data, [0 paddingNorthAndWest(2)], padMethod, 'pre');
end
if paddingSouthAndEast(1)
    data = padarray(data, [paddingSouthAndEast(1) 0], padMethod, 'post');
end
if paddingSouthAndEast(2)
    data = padarray(data, [0 paddingSouthAndEast(2)], padMethod, 'post');
end

if convertToCategorical
    data = categorical(data, 1:numel(origCategories), origCategories);
end

end
