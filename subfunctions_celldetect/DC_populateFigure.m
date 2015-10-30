function sFig = DC_populateFigure(handles)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%populate button pointers
	cellFields = fieldnames(handles);
	for intField=1:length(cellFields)
		strField = cellFields{intField};
		if length(strField) > 3 && strcmp(strField(1:3),'ptr')
			sFig.(strField) = handles.(strField);
		end
	end
	
	
	%allocate other variables
	sFig.strAnnotations = 'Hide';
	sFig.ButtonGroupAnnotations = handles.ButtonGroupAnnotations;
	sFig.ptrMainGUI = handles.output;
	sFig.vecSelectedObjects = [];
	sFig.sObject = [];
	sFig.ptrWindowHandle = [];
	sFig.ptrAxesHandle = [];
	sFig.imCurrent = [];
	sFig.imOriginal = [];
	sFig.imROI = [];
	
	
	%sFig.sObject(intObject).drawn
	%sFig.sObject(intObject).handles.marker
	%sFig.sObject(intObject).handles.lines
end