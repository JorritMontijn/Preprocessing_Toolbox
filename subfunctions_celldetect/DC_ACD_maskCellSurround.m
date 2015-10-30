function imOut = DC_ACD_maskCellSurround(imI, dblCenterX, dblCenterY, dblPixelSize, dblSizeROI, boolInvert)
	% Returns an image of original size that is black, except for the roi
	
	%% PERFORM STEPS
	% get image size
	[intMaxY, intMaxX] = size(imI);
	intCenterY = round(dblCenterY);
	intCenterX = round(dblCenterX);
	
	% get number of pixels of ROI-window radius
	intRadiusReq = round((dblSizeROI/dblPixelSize) / 2);
	intRadiusROI = min([intRadiusReq;...
		intMaxY-intCenterY-1;...
		intMaxX-intCenterX-1;...
		intCenterY-1;...
		intCenterX-1]);
	
	% set field of view range
	vecSubGrid = -intRadiusROI:intRadiusROI;
		
	%get subregion (defined by x,y, and ROI size)
	if boolInvert
		imI = 1-imI;
	end
	vecSelectY = vecSubGrid + intCenterY;
	vecSelectX = vecSubGrid + intCenterX;
	
	% make field of view mask circular
	[matGridX,matGridY]=meshgrid(vecSubGrid,vecSubGrid);
	
	%create 2D exponential within subregion
	dblGradientWidth = 1;
	matMask = exp(  -( ( ((matGridX/(intRadiusReq*dblGradientWidth)).^2) + ...
		((matGridY/(intRadiusReq*dblGradientWidth)).^2) ) )  );
	
	%get actual pixels from image, then blur & normalize
	imSub = imI(vecSelectY,vecSelectX);
	ptrFilterKernel = fspecial('disk', 3);
	imSub = imadjust(imfilter(imSub, ptrFilterKernel, 'conv', 'replicate'));
	
	%return same size image all black, except ROI; exponential filter
	imOut = zeros(size(imI));
	
	%multiplied by blurred, normalized original image pixels
	imOut(vecSelectY, vecSelectX) = imSub .* matMask;
end