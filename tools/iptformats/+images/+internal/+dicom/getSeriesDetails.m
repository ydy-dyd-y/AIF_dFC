function [seriesFilenames, spatialDetails, sliceDim] = getSeriesDetails(allFilenames)
%

% Copyright 2016-2017 The MathWorks, Inc.

if ischar(allFilenames)
    allFilenames = {allFilenames};
end

numImages = numel(allFilenames);

if numImages == 1
    filename = allFilenames{1};
    
    try
        d = images.internal.dicom.DICOMFile(filename);
    catch
        error(message('images:dicomread:loadFile', filename))
    end
    
    verifyImg = dicomread(filename);
    if ndims(verifyImg) == 2 %#ok<ISMAT>
        error(message('images:dicomread:notEnoughSlices'))
    end
    
    spatialDetails = images.internal.dicom.getSpatialDetailsForMultiframe(d);
    
    if (isempty(spatialDetails.PatientPositions))
        throwMissingPatientPosition(filename)
    end
    
    seriesFilenames = allFilenames;
    sliceDim = images.internal.dicom.findSortDimension(spatialDetails.PatientPositions);
    return
end

allPatientPositions = zeros(numImages, 3);
allPixelSpacings = zeros(numImages, 2);
allPatientOrientations = zeros(2, 3, numImages);
sliceDim = [];

if isempty(allFilenames)
    return
end

partOfSeries = false(numImages, 1);
previousSeriesInstanceUID = '';
previousPatientOrientation = [];

firstFilename = allFilenames{1};
for idx = 1:numel(allFilenames)
    filename = allFilenames{idx};
    
    try
        d = images.internal.dicom.DICOMFile(filename);
    catch
        if isdicom(filename)
            error(message('images:dicomread:loadFile', filename))
        else
            partOfSeries(idx) = false;
            continue
        end
    end
    
    try
        allPatientPositions(idx, :) = d.getAttribute(32, 50);  % (0020,0032)
    catch
        throwMissingPatientPosition(filename);
    end
    
    thisPatientOrientation = d.getAttribute(32, 55);  % (0020,0037)
    if isempty(thisPatientOrientation)
        throwMissingOrientation(filename);
    end
    
    try
        allPixelSpacings(idx, :) = d.getAttribute(40, 48);
    catch
        throwMissingPixelSpacing(filename);
    end
    
    thisSeriesInstanceUID    = d.getAttribute(32, 14);  % (0020,000E)
    
    partOfSeries(idx) = true;
    
    if isempty(previousSeriesInstanceUID)
        previousSeriesInstanceUID = thisSeriesInstanceUID;
    elseif ~isequal(previousSeriesInstanceUID, thisSeriesInstanceUID)
        error(message('images:dicomread:multipleSeries', ...
            firstFilename, filename))
    end
    
    if isempty(previousPatientOrientation)
        previousPatientOrientation = thisPatientOrientation;
    elseif ~isequal(previousPatientOrientation, thisPatientOrientation)
        error(message('images:dicomread:differentPatientOrientations'))
    end
    
    if ~isempty(thisPatientOrientation)
        allPatientOrientations(:, :, idx) = reshape(thisPatientOrientation, [3 2])';
    end
end

allPatientPositions = allPatientPositions(partOfSeries,:);
allPixelSpacings = allPixelSpacings(partOfSeries,:);

seriesFilenames = allFilenames(partOfSeries);

if isempty(seriesFilenames)
    error(message('images:dicomread:noSuccessfulReads'))
else
    [seriesFilenames, allPatientPositions, allPixelSpacings, allPatientOrientations, sliceDim] = sortSlices(seriesFilenames, allPatientPositions, allPixelSpacings, allPatientOrientations);
end

spatialDetails.PatientPositions = allPatientPositions;
spatialDetails.PixelSpacings = allPixelSpacings;
spatialDetails.PatientOrientations = allPatientOrientations;

end


function [seriesFilenames, allPatientPositions, allPixelSpacings, allPatientOrientations, sortDim] = sortSlices(seriesFilenames, allPatientPositions, allPixelSpacings, allPatientOrientations)

sortDim = images.internal.dicom.findSortDimension(allPatientPositions);
[~, sortIdx] = sort(allPatientPositions(:,sortDim));

seriesFilenames = seriesFilenames(sortIdx);
allPatientPositions = allPatientPositions(sortIdx,:);
allPixelSpacings = allPixelSpacings(sortIdx,:);
allPatientOrientations = allPatientOrientations(:,:,sortIdx);

end


function throwMissingPatientPosition(filename)

msg = message('images:dicomread:missingPatientPositions', filename);
ex = MException('images:dicomread:missingPatientPositions', ...
    strrep(msg.getString(), '\', '\\'));
throw(ex)

end

function throwMissingOrientation(filename)
msg = message('images:dicomread:missingPatientOrientations', filename);
ex = MException('images:dicomread:missingPatientOrientations', ...
    strrep(msg.getString(), '\', '\\'));
throw(ex)

end

function throwMissingPixelSpacing(filename)
msg = message('images:dicomread:missingPixelSpacing', filename);
ex = MException('images:dicomread:missingPixelSpacing', ...
    strrep(msg.getString(), '\', '\\'));
throw(ex)

end
