function DC_detectOverlap
	%detects overlap between object masks
	
	global CDG;
	global sRec;
	
	w = waitbar(0,'Please wait while improving cell-boundaries...');
	
	I = CDG.I{5}(:,:,2);
	pixelSize = sRec.xml.sData.dblActualVoxelSizeX / 1000; % micrometers
	roiSize = CDG.expected_cell_size ; % micrometers
	
	cc = 0;
	for c = 1:CDG.numCells
		if CDG.cells(c).type == 1 || CDG.cells(c).type == 2 || CDG.cells(c).type == 4
			cc = cc + 1;
			x_orig(cc) = [CDG.cells(c).Centroids.X];
			y_orig(cc) = [CDG.cells(c).Centroids.Y];
			type_orig(cc) = CDG.cells(c).type;
		end
	end
	
	for c = 1:CDG.nNonDetectedCentroids
		if CDG.NonDetectedCentroids(c).type == 1 || CDG.NonDetectedCentroids(c).type == 2 ...\
				|| CDG.NonDetectedCentroids(c).type == 4
			cc = cc + 1;
			x_orig(cc) = [CDG.NonDetectedCentroids(c).x];
			y_orig(cc) = [CDG.NonDetectedCentroids(c).y];
			type_orig(cc) = CDG.NonDetectedCentroids(c).type;
		end
	end
	
	nObjects = length(x_orig) ;
	ThresholdPercentageNeuron = 0.7 ;
	ThresholdPercentageAstrocyte = 0.78 ;
	ThresholdPercentageBloodvessel = 0.85 ;
	
	% get cell mask
	for c = 1:nObjects
		% mask cell
		Ic = ACD_maskCellSurround( I, x_orig(c), y_orig(c), pixelSize, roiSize, type_orig(c) );
		% get thresholded contour
		if type_orig(c) == 1
			Im{c} = ACD_thresholdContour( Ic, ThresholdPercentageNeuron );
		elseif type_orig(c) == 2
			Im{c} = ACD_thresholdContour( Ic, ThresholdPercentageAstrocyte );
		elseif type_orig(c) == 4
			Im{c} = ACD_thresholdContour( Ic, ThresholdPercentageBloodvessel );
		end
	end
	
	% remove overlap
	[ImRO, oList, toDelete] = ACD_removeOverlap( Im );
	nObjectsRO = length(ImRO);
	
	typelist = type_orig(toDelete==false);
	find(toDelete>0)
	CDG.cells = [];
	for c = 1:nObjectsRO
		[x, y] = ACD_getPerimeter( ImRO{c} );
		
		[labeledObject,nPolygons] = bwlabel(ImRO{c},4) ;  % also 8 possible
		
		% get properties of identified objects
		cellBasic  = regionprops(labeledObject, 'Basic') ;
		cellPixels = regionprops(labeledObject, 'PixelList') ;
		cellPerimeter = regionprops(labeledObject, 'ConvexHull') ;
		
		CDG.cells(c).ObjectLabel = c;
		CDG.cells(c).type = typelist(c);
		CDG.cells(c).Centroids.X = cellBasic(1).Centroid(1) ;
		CDG.cells(c).Centroids.Y = cellBasic(1).Centroid(2) ;
		CDG.cells(c).Size        = sum(sum(ImRO{c}>0)) ;
		CDG.cells(c).Perimeter   = cellPerimeter(1).ConvexHull ;
		CDG.cells(c).Body        = cellPixels(1).PixelList ;
		CDG.cells(c).Perimeter21px   = [];
		CDG.cells(c).Body21px        = [];
		CDG.cells(c).Luminance   = sum(sum( I( ImRO{c}>0 ) )) / sum(sum(ImRO{c}>0));
		CDG.cells(c).drawn       = 0;
		CDG.cells(c).selected       = 0;
		CDG.cells(c).handles.lines = [];
		CDG.cells(c).handles.marker = 0;
	end
	CDG.numCells = nObjectsRO;
	
	CDG.nNonDetectedCentroids = 0;
	CDG.NonDetectedCentroids.ObjectLabel = [];
	CDG.NonDetectedCentroids.type = [];
	CDG.NonDetectedCentroids.x = [];
	CDG.NonDetectedCentroids.y = [];
	CDG.NonDetectedCentroids.drawn = [];
	CDG.NonDetectedCentroids.handle = [];
	
	redraw;
	close(w);
end