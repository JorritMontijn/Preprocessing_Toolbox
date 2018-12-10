function [intRealX,intRealY] = DC_WaitForClick(sFig)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
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
end

