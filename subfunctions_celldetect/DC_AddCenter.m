function DC_AddCenter
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
			intPresence = get(sFig.ptrListPresence, 'Value');
			intRespType = get(sFig.ptrListRespType, 'Value');
		catch
			return;
		end
	end
	
	boolRepeat = true;
	while boolRepeat
		
		%set help text
		cellText{1} = ['Help info on center drawing'] ;
		cellText{2} = [''] ;
		cellText{3} = ['Now drawing: ' sDC.metaData.cellType{intType}] ;
		cellText{4} = [''] ;
		cellText{5} = ['Left-click any point in the figure to set a center-point for later boundary detection'] ;
		cellText{6} = ['Right-click when done with center drawing'];
		DC_updateTextInformation(cellText)
		
		%put active fig on top
		figure(sFig.ptrWindowHandle);
		
		%wait for mouse click
		while waitforbuttonpress ~= 0, end
		
		%check what kind of button press
		s = get(gcf, 'SelectionType' );
		
		if strcmp(s, 'alt')
			%right-click: return
			return;
		else
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
			
			%get location of the click
			XY = get(sFig.ptrWindowHandle, 'CurrentPoint');
			X = XY(1);%from right
			Y = XY(2);%from bottom
			
			%get magnification  & zoom
			dblMagnification = str2double(get(sFig.ptrEditMagnification,'String'))/100;
			
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
			
			
			% get cell nr
			intObjects = numel(sDC.ROI);
			intNewObject = intObjects+1;
			
			% save roi
			sDC.ROI(intNewObject).intPresence = intPresence;
			sDC.ROI(intNewObject).intRespType = intRespType;
			sDC.ROI(intNewObject).intType = intType;
			sDC.ROI(intNewObject).intCenterX = intRealX;
			sDC.ROI(intNewObject).intCenterY = intRealY;
			
			sDC.ROI(intNewObject).matPerimeter = [];
			sDC.ROI(intNewObject).matMask = [];
			
			
			sFig.sObject(intNewObject).drawn = 0;
			sFig.sObject(intNewObject).handles.marker = [];
			sFig.sObject(intNewObject).handles.lines = [];
			
			%redraw
			DC_redraw(0);
		end
	end
end