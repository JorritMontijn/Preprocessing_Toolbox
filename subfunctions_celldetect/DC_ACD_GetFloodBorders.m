function [matFloodBorders,matMergedCentroids,vecMergedRadii] = DC_ACD_GetFloodBorders(im1D,dblMergeThreshold)
	%% build borders based on proximity to flood border
	%get step size
	im1D = (im1D + circshift(im1D,[1 0]))/2;
	im1D = imnorm(im1D);
	vecStepSizes = sort(unique(roundi(diff(unique(im1D(:))),10)),'ascend');
	dblStep = vecStepSizes(2); %[first is 0, second is smallest step size]
	
	%build filters
	matFilt = [-1 -1 -1; -1 8 -1; -1 -1 -1];
	
	%check threshold
	if ~exist('dblMergeThreshold','var') || isempty(dblMergeThreshold),dblMergeThreshold=30;end
	
	%% loop
	matFloodBordersLow = zeros(size(im1D));
	matFloodBordersHigh = zeros(size(im1D));
	for dblThreshold=dblStep:dblStep:max(im1D(:))
		%get flood fill
		imBW = im1D <= dblThreshold;
		
		% get boundaries
		imBorder = imfilter(double(imBW),matFilt,'replicate');
		imDilate = imdilate(abs(imBorder),strel('disk',3));
		
		%assign 1s to all pixels within 5 pixels of boundary
		matFloodBordersLow = matFloodBordersLow + imDilate*(dblStep/dblThreshold);
		matFloodBordersHigh = matFloodBordersHigh + imDilate;
	end
	matFloodBorders = imnorm(matFloodBordersLow) + imnorm(matFloodBordersHigh);
	
	%% detect circles
	if nargin > 1
		[matMergedCentroids,vecMergedRadii] = DC_ACD_GetDisks(matFloodBorders,dblMergeThreshold);
	end
	
	%plot
	%figure,imshow(matFloodBorders);hold on;
	%viscircles(matMergedCentroids, vecMergedRadii,'EdgeColor','r');
end