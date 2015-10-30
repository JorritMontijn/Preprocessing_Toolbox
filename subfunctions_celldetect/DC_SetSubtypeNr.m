function DC_SetSubtypeNr
	%globals
	global sDC;
	global sFig;
	global sSetSubtypeNr;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig)
		return;
	else
		try
			intType = get(sFig.ptrListCellType, 'Value');
		catch
			return;
		end
	end
	
	%get data on selected object
	if length(sFig.vecSelectedObjects) == 1
		intObject = sFig.vecSelectedObjects;
		intSubTypeNr = getSubTypeNr(sDC,intObject);
	else
		if sum(sFig.vecSelectedObjects) > 1
			cellText{1} = ['number of selected objects: ' num2str(length(sFig.vecSelectedObjects))] ;
			cellText{2} = [''] ;
			cellText{3} = ['more than 1 cell selected, cannot change nr'] ;
		elseif sum(sFig.vecSelectedObjects) < 1
			cellText{1} = ['no cells selected'] ;
			cellText{2} = [''] ;
			cellText{3} = ['cannot change nr of non-selected object'] ;
		end
		%set figure text
		set(sFig.ptrTextInformation, 'string', cellText );
		drawnow;
		return
	end
	
	%set for exit
	sSetSubtypeNr.boolRunning = true;
	sSetSubtypeNr.intNewSubtypeNr = intSubTypeNr;
	
	%put GUI on top
	figure(sFig.ptrMainGUI);
	figure(sFig.ptrWindowHandle);
	
	%run GUI
	DC_SetSubtypeNrGUI;
	
	%wait for GUI to close
	while sSetSubtypeNr.boolRunning,pause(0.1);end
	
	%check if new value is different
	if sSetSubtypeNr.intNewSubtypeNr ~= intSubTypeNr
		%get data
		intSubType = sDC.ROI(intObject).intType;
		intNewSubtypeNr = sSetSubtypeNr.intNewSubtypeNr;
		intSubNr = 0;
		intObjects = numel(sDC.ROI);
		boolRunning = true;
		boolFound = false;
		intSearchObject = 0;
		while boolRunning
			intSearchObject = intSearchObject + 1;
			if intSearchObject > intObjects
				boolRunning = false;
			elseif sDC.ROI(intSearchObject).intType == intSubType
				intSubNr = intSubNr + 1;
				if intSubNr == intNewSubtypeNr
					boolRunning = false;
					boolFound = true;
				end
			end
		end
		
		%check if found
		if ~boolFound
			cellText{1} = ['Subtype number [' num2str(intNewSubtypeNr) '] is out of bounds'] ;
			cellText{2} = [''] ;
			cellText{3} = ['cannot change nr'] ;
			
			%set figure text
			set(sFig.ptrTextInformation, 'string', cellText );
			drawnow;
			return;
		end
		
		%get new object number
		intNewObjectNr = intSearchObject;
		
		%switch objects
		sDC.ROI(intObjects+1) = sDC.ROI(intObject);
		sFig.sObject(intObjects+1) = sFig.sObject(intObject);
		
		sDC.ROI(intObject) = sDC.ROI(intNewObjectNr);
		sFig.sObject(intObject) = sFig.sObject(intNewObjectNr);
		
		sDC.ROI(intNewObjectNr) = sDC.ROI(intObjects+1);
		sFig.sObject(intNewObjectNr) = sFig.sObject(intObjects+1);
		
		sDC.ROI = sDC.ROI(1:intObjects);
		sFig.sObject = sFig.sObject(1:intObjects);
		
		%set msg
		cellText{1} = ['Subtype number [' num2str(intSubTypeNr) '] switched with [' num2str(intNewSubtypeNr) '] '] ;
		cellText{2} = [''] ;
		cellText{3} = [''] ;
		
		%change selection vector
		sFig.vecSelectedObjects = intNewObjectNr;
		
		%set figure text
		DC_updateTextInformation(cellText)
	end
end
function intSubTypeNr = getSubTypeNr(sDC,intObject)
	intType = sDC.ROI(intObject).intType;
	intSubTypeNr = 1;
	for intObjectCounter=1:(intObject-1)
		if intType == sDC.ROI(intObjectCounter).intType
			intSubTypeNr = intSubTypeNr + 1;
		end
	end
end