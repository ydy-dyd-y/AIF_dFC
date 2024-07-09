classdef MatrixAdapter < images.internal.adapters.BigImageAdapter
    
    % Copyright 2019 The MathWorks, Inc.
    
    properties
        Data
    end
    
    methods
        function obj = MatrixAdapter(varName, data)
            obj.DataSource = varName;
            obj.Data = data;
            
            obj.Height = size(data,1);
            obj.Width = size(data,2);
            obj.IOBlockSize = [1024 1024];
            obj.PixelDatatype = class(data);
            obj.Channels = size(data,3);
        end
        
        
        function data = readBlock(obj, level, blockStart)
            % Only single level variables are supported
            assert(level==1);
            
            blockEnd = blockStart + obj.IOBlockSize - 1;
            blockEnd = min(blockEnd, [size(obj.Data,1), size(obj.Data,2)]);
            data =  obj.Data(blockStart(1):blockEnd(1),...
                blockStart(2):blockEnd(2),:);
            
            obj.updateNNZ(level, blockStart, data, false);
        end
    end
end

