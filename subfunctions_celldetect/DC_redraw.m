function DC_redraw(varargin)
	%DC_redraw Sets and redraws windows
	%   DC_redraw([redrawImage=true])
	
	% check if only cells have to be redrawn
	redrawImage = true;
	if nargin > 0
		redrawImage = varargin{1};
	end
	
	%get structures
	global sRec;
	global sDC;
	global sFig;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig) || isempty(sRec)
		return;
	else
		try
			%get current image
			intImSelected = get(sFig.ptrListSelectImage,'Value');
		catch %#ok<CTCH>
			return;
		end
	end
	
	%set msg
	cellTextOld = get(sFig.ptrTextInformation, 'string');
	cellText{1} = ['Redrawing...'] ;
	DC_updateTextInformation(cellText)
	
	%intImSource = sFig.vecSource(intImSelected)
	imThis = sFig.cellIm{intImSelected};
	
	% get the current zoom value and adjust the image
	dblImageMagnification = str2double(get(sFig.ptrEditMagnification,'String'))/100;
	if dblImageMagnification > 10, dblImageMagnification = 10;end
	newX = round(sRec.sProcLib.x * dblImageMagnification);
	newY = round(sRec.sProcLib.y * dblImageMagnification);
	newImage = imresize(imThis, [newY newX ]);
	sFig.imCurrent = newImage;
	sFig.imOriginal = imThis;
	
	%get number of objects
	intObjects = numel(sDC.ROI);
	
	% Draw the image if requested
	if redrawImage ~= 0
		%check if figure is still there
		try
			vGet = get(sFig.ptrAxesHandle);
		catch %#ok<CTCH>
			vGet = [];
		end
		
		%make new figure
		if redrawImage == 2 || isempty(vGet)
			%close figure if old one still present
			if ~isempty(vGet)
				close(sFig.ptrWindowHandle);
			end
			
			% create figure
			sFig.ptrWindowHandle = figure;
			set(sFig.ptrWindowHandle, 'Units', 'pixels' );
			set(sFig.ptrWindowHandle, 'position',[0 0 newX+10 newY+10] );
			
			sFig.ptrAxesHandle = axes;
			set(sFig.ptrAxesHandle, 'Units', 'pixels' );
			set(sFig.ptrAxesHandle, 'position',[5 5 newX newY]);
		else
			%set active figure
			figure(sFig.ptrWindowHandle);
		end
		
		% display the image
		imshow(newImage);
		axis xy;
		
		% reset image-drawn flags
		sFig.sObject = [];
		for intObject = 1:intObjects
			sFig.sObject(intObject).drawn = 0;
		end
		
		% connect the PlotClick_Callback function to handle clicks in the image
		set(gcf, 'WindowButtonDownFcn', {@DC_PlotClick_Callback});
	end
	
	%update total objects in figure
	set(sFig.ptrTextTotalObjects, 'String', intObjects);
	set(sFig.ptrTextTotalType, 'String', DC_getTotalType());
	
	%set active figure
	figure(sFig.ptrWindowHandle);
	
	%check draw type
	intDrawType = get(sFig.ptrListSelectDrawType,'Value');
	cellDrawTypes = get(sFig.ptrListSelectDrawType,'String');
	strDrawType = cellDrawTypes{intDrawType};
	boolDrawBorder = strcmpi(strDrawType,'Border');
	
	% draw cells if detected
	for intObject = 1:intObjects
		
		% don't update cells already drawn..
		if sFig.sObject(intObject).drawn == 0
			%check if selected
			if ismember(intObject,sFig.vecSelectedObjects)
				intLineWidth = 4;
				intMarkerSize = 24;
			else
				intLineWidth = 2;
				intMarkerSize = 16;
			end
			
			%get type
			intType = sDC.ROI(intObject).intType;
			strType = sDC.metaData.cellType{intType};
			
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
			if isfield(sDC.ROI(intObject),'matPerimeter') && ~isempty(sDC.ROI(intObject).matPerimeter)
				%draw outline
				for p = 1:length(sDC.ROI(intObject).matPerimeter)-1
					
					sFig.sObject(intObject).handles.lines(p) = ...
						line(sDC.ROI(intObject).matPerimeter(p:p+1, 1) * dblImageMagnification, ...
						sDC.ROI(intObject).matPerimeter(p:p+1, 2) * dblImageMagnification, ...
						'Color', vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth );
				end
				sFig.sObject(intObject).handles.lines(p+1) = ...
					line([sDC.ROI(intObject).matPerimeter(1,   1) ...
					sDC.ROI(intObject).matPerimeter(end, 1) ] * dblImageMagnification, ...
					[sDC.ROI(intObject).matPerimeter(1,   2) ...
					sDC.ROI(intObject).matPerimeter(end, 2) ] * dblImageMagnification, ...
					'Color',  vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth );
			%}
			if isfield(sDC.ROI,'matMask') && ~isempty(sDC.ROI(intObject).matMask) && boolDrawBorder
				sO = bwboundaries(sDC.ROI(intObject).matMask);
				
				%draw outline
				for p = 1:length(sO{1})-1
					sFig.sObject(intObject).handles.lines(p) = ...
						line(sO{1}(p:p+1, 2) * dblImageMagnification, ...
						sO{1}(p:p+1, 1) * dblImageMagnification, ...
						'Color', vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth );
				end
				sFig.sObject(intObject).handles.lines(p+1) = ...
					line([sO{1}(1,   2) ...
					sO{1}(end, 2) ] * dblImageMagnification, ...
					[sO{1}(1,   1) ...
					sO{1}(end, 1) ] * dblImageMagnification, ...
					'Color',  vecColor, 'LineStyle', '-', 'LineWidth', intLineWidth );
				
			elseif isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
				%draw center
				x = sDC.ROI(intObject).intCenterX * dblImageMagnification;
				y = sDC.ROI(intObject).intCenterY * dblImageMagnification;
				sFig.sObject(intObject).handles.marker = ...
					line([x x], [y y], 'color', vecColor, ...
					'Marker', '.',   'MarkerSize', intMarkerSize) ;
			else
				fprintf('Uh-oh... something went wrong with object %d\n',intObject)
			end
			
			%check if annotation need to be drawn
			if strcmp(sFig.strAnnotations,'Show') && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
				x = sDC.ROI(intObject).intCenterX * dblImageMagnification;
				y = sDC.ROI(intObject).intCenterY * dblImageMagnification;
				sFig.sObject(intObject).handles.text = text(x+5*dblImageMagnification,y+5*dblImageMagnification,num2str(getSubTypeNr(sDC,intObject)),'Color',vecColor,'FontSize',12,'FontWeight','demi');
			end
			
			% update drawn-flag
			sFig.sObject(intObject).drawn = 1;
		end
	end
	drawnow;
	
	%set msg
	DC_updateTextInformation(cellTextOld)
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