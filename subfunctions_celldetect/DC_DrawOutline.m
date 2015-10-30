function DC_DrawOutline(varargin)
	%waits for user to click anywhere on window to zoom in; then creates
	%zoomed-in figure allowing for more precise selection of boundaries
	
	%get globals
	global sDC;
	global sFig;
	global sRec;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig) || isempty(sRec)
		return;
	else
		try
			intType = get(sFig.ptrListCellType, 'Value');
			axesSize = get(sFig.ptrAxesHandle, 'Position');
			intPresence = get(sFig.ptrListPresence, 'Value');
		catch
			return;
		end
	end
	
	% Get margins of screen
	axesSize = get(sFig.ptrAxesHandle, 'Position');
	xZero = axesSize(1);
	yZero = axesSize(2);
	width = axesSize(3);
	height = axesSize(4);
	
	% get margins of plot-window
	sizeI = size(sFig.imCurrent);
	screenWidth = sizeI(2);
	screenHeight = sizeI(1);
	
	%get magnification  & zoom
	dblMagnification = str2double(get(sFig.ptrEditMagnification,'String'))/100;
	dblZoom = str2double(get(sFig.ptrEditZoom,'String'))/100;
	
	
	if nargin == 0
		%put active fig on top
		figure(sFig.ptrWindowHandle);
		
		%wait for mouse click
		while waitforbuttonpress ~= 0, end
		
		
		%get location of the click
		XY = get(sFig.ptrWindowHandle, 'CurrentPoint');
		X = XY(1);%from right
		Y = XY(2);%from bottom
		
		
		% Calculate the coordinates of the click related to the image
		% coordinates
		dblRealX = ( ((X - xZero) * screenWidth ) / width  )/dblMagnification;
		dblRealY = ( ((Y - yZero) * screenHeight) / height )/dblMagnification;
		intRealX = round(dblRealX);
		intRealY = round(dblRealY);
		
		% check if click was in correct Axes
		get(gcf,'CurrentAxes');
		if sFig.ptrAxesHandle ~= get(gcf,'CurrentAxes')
			return;
		elseif isempty(get(gcf,'CurrentAxes'))
			return;
		end
		
		%set object variables
		intObjects = numel(sDC.ROI);
		intNewObject = intObjects+1;
		intNewType = get(sFig.ptrListCellType, 'Value');
		intNewPresence = get(sFig.ptrListPresence, 'Value');
		intNewRespType = get(sFig.ptrListRespType, 'Value');
	else
		%get object nr
		intNewObject = varargin{1};
		
		%get variables
		intNewRespType = sDC.ROI(intNewObject).intRespType;
		intNewPresence = sDC.ROI(intNewObject).intPresence;
		intNewType = sDC.ROI(intNewObject).intType;
		intRealX = round(sDC.ROI(intNewObject).intCenterX);
		intRealY = round(sDC.ROI(intNewObject).intCenterY);
		
		%clear mask
		sDC.ROI(intNewObject).matPerimeter = [];
		sDC.ROI(intNewObject).matMask = [];
	end
	
	%set help text
	cellText{1} = ['Help info on outline drawing'] ;
	cellText{2} = [''] ;
	cellText{3} = ['Now drawing: ' sDC.metaData.cellType{intNewType}] ;
	cellText{4} = [''] ;
	cellText{5} = ['Left-click any point in the figure to set a corner'] ;
	cellText{6} = [''];
	cellText{7} = ['Right-click when finished; or to cancel when less than 3 corners have been set'] ;
	DC_updateTextInformation(cellText)
	
	% calculate coordinates of a zoomed in  area
	widthheight = round(sRec.sProcLib.x / dblZoom);
	x(1) = intRealX - round(widthheight/2);
	x(2) = intRealX + round(widthheight/2);
	y(1) = intRealY - round(widthheight/2);
	y(2) = intRealY + round(widthheight/2);
	
	%move zoomed-in area if outside original image
	if x(2) - x(1) >= sRec.sProcLib.x
		x = [1 sRec.sProcLib.x];
	else
		if x(1) < 1
			x = x - x(1) + 1;
		end
		if x(2) > sRec.sProcLib.x
			x = x - (x(2) - sRec.sProcLib.x);
		end
	end
	if y(2) - y(1) >= sRec.sProcLib.y
		y = [1 sRec.sProcLib.y];
	else
		if y(1) < 1
			y = y - y(1) + 1;
		end
		if y(2) > sRec.sProcLib.y
			y = y - (y(2) - sRec.sProcLib.y);
		end
	end
	
	%calculate required number of pixels
	newzoomX = round((sRec.sProcLib.x / dblZoom) * dblMagnification);
	newzoomY = round((sRec.sProcLib.y / dblZoom) * dblMagnification);
	
	% open window with zoom-in
	ptrZoomFig = figure;
	set( ptrZoomFig, 'units', 'pixels' );
	set( ptrZoomFig, 'position',[50 50 sRec.sProcLib.x*dblMagnification+10 sRec.sProcLib.y*dblMagnification+10] );
	
	ptrZoomAxes = axes;
	set( ptrZoomAxes, 'units', 'pixels' );
	set( ptrZoomAxes, 'position', ...
		[5 5 sRec.sProcLib.x*dblMagnification sRec.sProcLib.y*dblMagnification] );
	
	%create zoomed image
	vecSelectY = y(1):y(2);
	vecSelectX = x(1):x(2);
	if ndims(sFig.imOriginal) == 3
		zoomI = sFig.imOriginal(vecSelectY,vecSelectX,:);
	else
		zoomI = sFig.imOriginal(vecSelectY,vecSelectX);
	end
	newImage = imresize(zoomI, [ newzoomY newzoomX ] );
	imshow(newImage);
	axis xy;
	
	% draw cells in zoom window
	intObjects = numel(sDC.ROI);
	for intObject = 1:intObjects
		
		%check if selected
		if ismember(intObject,sFig.vecSelectedObjects)
			intLineWidth = 4;
			intMarkerSize = 20;
		else
			intLineWidth = 2;
			intMarkerSize = 16;
		end
		
		%get type
		intType = sDC.ROI(intObject).intType;
		strType = sDC.metaData.cellType{sDC.ROI(intObject).intType};
		
		%assign color
		vecColor = sDC.metaData.cellColor{intType};
		
		%if neuron
		if ismember(intType,sDC.metaData.vecNeurons)
			%get presence
			intPresence = sDC.ROI(intObject).intPresence;
			strPresence = sDC.metaData.cellPresence{intPresence};
			
			%change amount of blue
			if strcmp(strPresence,'include')
				vecColor(end) = 0;
			elseif strcmp(strPresence,'present')
				vecColor = vecColor * 0.7;
				vecColor(end) = 0.7;
			elseif strcmp(strPresence,'absent')
				vecColor(end) = 1;
			end
		end
		
		% draw objects
		%{
		if ~isempty(sDC.ROI(intObject).matPerimeter)
			%draw outline
			for p = 1:length(sDC.ROI(intObject).matPerimeter)-1
				xVec = (sDC.ROI(intObject).matPerimeter(p:p+1, 1)- x(1))*dblMagnification;
				yVec = (sDC.ROI(intObject).matPerimeter(p:p+1, 2)- y(1))*dblMagnification;
				sFig.sObject(intObject).handles.lines(p) = ...
					line(xVec,yVec,'Color', vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth);
			end
			sFig.sObject(intObject).handles.lines(p+1) = ...
				line(([sDC.ROI(intObject).matPerimeter(1,   1) ...
				sDC.ROI(intObject).matPerimeter(end, 1) ] - x(1)) *dblMagnification, ...
				([sDC.ROI(intObject).matPerimeter(1,   2) ...
				sDC.ROI(intObject).matPerimeter(end, 2) ] - y(1)) *dblMagnification, ...
				'Color',  vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth );
			
		elseif isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
			%draw center
			xV = (sDC.ROI(intObject).intCenterX - x(1)) *dblMagnification;
			yV = (sDC.ROI(intObject).intCenterY - y(1)) *dblMagnification;
			sFig.sObject(intObject).handles.marker = ...
				line([xV xV], [yV yV], 'color', vecColor, ...
				'Marker', '.',   'MarkerSize', intMarkerSize) ;
		end
		%}
		if isfield(sDC.ROI,'matMask') && ~isempty(sDC.ROI(intObject).matMask)
			sO = bwboundaries(sDC.ROI(intObject).matMask);
			
			for p = 1:length(sO{1})-1
				xVec = (sO{1}(p:p+1, 2)- x(1))*dblMagnification;
				yVec = (sO{1}(p:p+1, 1)- y(1))*dblMagnification;
				line(xVec,yVec,'Color', vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth);
			end
			line(([sO{1}(1,   2) ...
				sO{1}(end, 2) ] - x(1)) *dblMagnification, ...
				([sO{1}(1,   1) ...
				sO{1}(end, 1) ] - y(1)) *dblMagnification, ...
				'Color',  vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth );
			
		elseif isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
			%draw center
			xV = (sDC.ROI(intObject).intCenterX - x(1)) *dblMagnification;
			yV = (sDC.ROI(intObject).intCenterY - y(1)) *dblMagnification;
			line([xV xV], [yV yV], 'color', vecColor, ...
				'Marker', '.',   'MarkerSize', intMarkerSize) ;
		end
	end
	
	% Get margins of roi screen
	axesSize = get(ptrZoomAxes, 'Position');
	xZero = axesSize(1);
	yZero = axesSize(2);
	xMax = axesSize(3);
	yMax = axesSize(4);
	
	% get margins of zoom-window
	sizeI = size(zoomI);
	screenWidth = sizeI(2);
	screenHeight = sizeI(1);
	
	% read in outline of new roi
	vecRealX = [];
	vecRealY = [];
	w = 0;
	count = 0;
	boolCancel = 0;
	while w == 0
		% wait for input
		w = waitforbuttonpress;
		
		s = get( gcf, 'SelectionType' );
		
		% check if click was in correct Axes
		if ptrZoomAxes ~= get(gcf,'CurrentAxes')
			boolCancel = 1;
			break;
		end
		
		% if mouse click
		if strcmp(s, 'normal')
			count = count + 1;
			
			% get location of the click
			XY = get(ptrZoomFig, 'CurrentPoint');
			X = XY(1);
			Y = XY(2);
			
			% Calculate the coordinates of the click related to the image
			% coordinates
			X = ( ((X - xZero) * screenWidth ) / xMax  );
			Y = ( ((Y - yZero) * screenHeight) / yMax );
			roiX(count) = round(X);
			roiY(count) = round(Y);
			
			intRealX = round(X)+x(1)-1;
			intRealY = round(Y)+y(1)-1;
			
			vecRealX(count) = intRealX;
			vecRealY(count) = intRealY;
			
			%assign color
			vecColorROI = sDC.metaData.cellColor{intNewType};

			%if neuron
			if ismember(intType,sDC.metaData.vecNeurons)
				%get presence
				strPresence = sDC.metaData.cellPresence{intNewPresence};

				%change amount of blue
				if strcmp(strPresence,'include')
					vecColorROI(end) = 0;
				elseif strcmp(strPresence,'present')
					vecColorROI = vecColorROI * 0.7;
					vecColorROI(end) = 0.7;
				elseif strcmp(strPresence,'absent')
					vecColorROI(end) = 1;
				end
			end
			
			% draw line surrounding the cell
			dblLineWidthROI = 3;
			if count > 1
				xVec = ([roiX(count-1) roiX(count)])*dblMagnification;
				yVec = ([roiY(count-1) roiY(count)])*dblMagnification;
				
				line(xVec,yVec,'Color', vecColorROI, 'LineStyle', '-', 'LineWidth', dblLineWidthROI );
			end
		elseif strcmp(s, 'alt')
			%on right mouse click, drawing is finished
			w = 1;
			
			%if less than 3 points have been clicked, cancel
			if count < 3
				boolCancel = true;
			end
		end
	end
	%close figure
	close(ptrZoomFig);
	figure(sFig.ptrWindowHandle);
	
	%update information window
	DC_updateTextInformation;
	
	%if cancelled
	if boolCancel == 1
		return;
	end
	
	%get mask
	matMask = poly2mask(vecRealX, vecRealY, sRec.sProcLib.x, sRec.sProcLib.y );
	
	% save roi and close window
	sDC.ROI(intNewObject).intPresence = intNewPresence;
	sDC.ROI(intNewObject).intRespType = intNewRespType;
	sDC.ROI(intNewObject).intType = intNewType;
	sDC.ROI(intNewObject).intCenterX = mean(vecRealX);
	sDC.ROI(intNewObject).intCenterY = mean(vecRealY);
	
	%sDC.ROI(intNewObject).matPerimeter = [vecRealX; vecRealY]';
	sDC.ROI(intNewObject).matMask = matMask;
	
	
	sFig.sObject(intNewObject).drawn = false;
	sFig.sObject(intNewObject).handles.marker = [];
	sFig.sObject(intNewObject).handles.lines = [];
end