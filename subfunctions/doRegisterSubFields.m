%recursive locally affine subfield reregistration algorithm to reposition
%ROI locations between recordings of the same population of cells
function matDataPoints = doRegisterSubFields(imRef,im2,sParams)
	%check input
	if nargin < 2,error([mfilename ':NoImages'],'Insufficient input: registration images were not supplied');end
	
	%set params
	if nargin < 3,sParams = struct;end
	if isfield(sParams,'boolDoPlot'),boolDoPlot=sParams.boolDoPlot;else boolDoPlot=false;end
	
	%get images; transform to double
	if boolDoPlot
		figure,imshow(imRef),title('Reference')
		figure,imshow(im2),title(sprintf('Level %d',0))
	end
	
	%set parameters
	intMinSize = 128; %minimum subfield size (has to be sufficiently large for correct registration)
	intFullSize = min(size(imRef,1),size(imRef,2)); %full size is equal to the shorter dimension of the reference image
	if round(log2(intFullSize)) ~= log2(intFullSize), error([mfilename ':ImSizeNotPowerOf2'],'Image size must be a power of two');end %if image size is not a power of two, the factor 2 subfield size regression will fail
	intLevels = log2(intFullSize)-log2(intMinSize)+1; %calculate number of levels to process, based on the full image size and the minimum subfield size
	
	%pre-allocate translation matrices
	matTranslations = zeros(size(imRef,1),size(imRef,2),2,intLevels); %matTranslations(xPix,yPix,intPlane,intLevel);intPlane=1 => x; intPlane=2 => y
	matTranslated = zeros(size(imRef,1),size(imRef,2),intLevels); %matTranslations(xPix,yPix,intPlane,intLevel);intPlane=1 => x; intPlane=2 => y
	
	%start recursive subfield registration algorithm
	for intLevel=1:intLevels %loop through the number of required subfield regressions
		intSize = 2^(log2(intFullSize)-intLevel+1); %get size of current subfields
		intPrevLevel = max(1,intLevel-1); %get index for pre-displacement data
		
		intMaxRows=size(imRef,1)/intSize; %calculate number of rows (y) to process
		intMaxColumns=size(imRef,2)/intSize; %calculate number of columns (x) to process
		
		%pre-allocate output list
		if intLevel == intLevels
			intDataPoints = intMaxRows * intMaxColumns;
			matDataPoints = zeros(intDataPoints,4);
		end
		
		%loop through rows and columns
		intDataPoint = 0;
		for intRow=1:intMaxRows %loop through subfield rows (y)
			intStartY = (intRow-1)*intSize + 1; %get y pixel start location
			intStopY = intStartY + intSize - 1; %get y pixel stop location
			dblCenterY = (intStartY + intStopY)/2;
			for intCol=1:intMaxColumns  %loop through subfield columns (x)
				intDataPoint = intDataPoint + 1;
				intStartX = (intCol-1)*intSize + 1; %get x pixel start location
				intStopX = intStartX + intSize - 1; %get x pixel stop location
				dblCenterX = (intStartX + intStopX)/2;
				
				%get translations from previous level as starting point
				matDisplaceX = matTranslations(intStartY:intStopY,intStartX:intStopX,1,intPrevLevel); %get x pre-diplacement values for the pixels corresponding to current subfield
				intPreDisplaceX = round(mean(matDisplaceX(:))); %calculate mean x pre-displacement for subfield
				matDisplaceY = matTranslations(intStartY:intStopY,intStartX:intStopX,1,intPrevLevel); %get y pre-diplacement values for the pixels corresponding to current subfield
				intPreDisplaceY = round(mean(matDisplaceY(:))); %calculate mean y pre-displacement for subfield
				
				%get subfields from images
				imSub2 = circshift(im2,[intPreDisplaceX intPreDisplaceY]); %x-y translated to-be-registered image with pre-displacement values from previous level
				imSub2 = imSub2(intStartY:intStopY,intStartX:intStopX); %get the subfield from the pre-translated to-be-registered image
				
				imSubRef = imRef(intStartY:intStopY,intStartX:intStopX); %get the unaltered reference subfield from the reference image
				
				%calculate required registration
				[output matImReg] = dftregistration(fft2(imSubRef),fft2(imSub2), 100); %register the two subfields
				
				dblDisplacementX = output(3); %put x displacement in variable
				dblDisplacementY = output(4); %put y displacement in variable
				
				%put in matrix
				matTranslations(intStartY:intStopY,intStartX:intStopX,1,intLevel) = dblDisplacementX + intPreDisplaceX; %put total registered x-displacement in this level's subfield region
				matTranslations(intStartY:intStopY,intStartX:intStopX,2,intLevel) = dblDisplacementY + intPreDisplaceY; %put total registered y-displacement in this level's subfield region
				
				%put data into list
				if intLevel == intLevels
					matDataPoints(intDataPoint,1) = dblCenterX; %x pixel location
					matDataPoints(intDataPoint,2) = dblCenterY; %y pixel location
					matDataPoints(intDataPoint,3) = dblDisplacementX + intPreDisplaceX; %x pixel translation
					matDataPoints(intDataPoint,4) = dblDisplacementY + intPreDisplaceY; %y pixel translation
				end
				
				%create translated matrix
				imSubTranslated = circshift(im2,round([(dblDisplacementX + intPreDisplaceX) (dblDisplacementY + intPreDisplaceY)])); %get translated to-be-registered image from data gathered for this level
				imSubTranslated = imSubTranslated(intStartY:intStopY,intStartX:intStopX); %get the translated subfield
				matTranslated(intStartY:intStopY,intStartX:intStopX,intLevel) = imSubTranslated; %put it into the proper coordinates of the translated matrix
			end
		end
		if boolDoPlot,figure,imshow(matTranslated(:,:,intLevel)),title(sprintf('Level %d',intLevel));end %make figure of translated to-be-registered image
	end
end