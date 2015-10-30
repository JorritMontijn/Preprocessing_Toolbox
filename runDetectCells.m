%% starting function
function varargout = runDetectCells(varargin)
	% runDetectCells Detect cells with GUI
	%
	%	Version 2.1 [2014-05-16]
	%	2013-03-14; Created by Jorrit Montijn
	%	2013-10-21; updated cross-recording alignment; incorporated
	%	custom-built recursive locally affine registration algorithm to
	%	automatically realign cell bodies between recordings; removed
	%	matPerimeter field from ROI structure, so everything now works with
	%	masks to reduce superfluous information and avoid possible bugs due
	%	to a two-variable based location storage system [by JM]
	%	2013-11-08; added pre-processing pixel-based selectivity/response
	%	calculation output to blue channel image for increased information
	%	during cell body selection and exclusion process [by JM]
	%	2014-01-27; added 'designation' as extra ROI property, including
	%	'include', 'absent', 'present', 'tuned', 'responsive' [by JM]
	%	2014-02-13; changed 'designation' to two separate properties,
	%	including responsiveness and presence [by JM] 
	%	2014-05-16; made several bug fixes and GUI changes [by JM] 
	
	%set tags
	%#ok<*INUSL>
	%#ok<*INUSD>
	
	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
		'gui_Singleton',  gui_Singleton, ...
		'gui_OpeningFcn', @runDetectCells_OpeningFcn, ...
		'gui_OutputFcn',  @runDetectCells_OutputFcn, ...
		'gui_LayoutFcn',  [] , ...
		'gui_Callback',   []);
	if nargin && ischar(varargin{1})
		gui_State.gui_Callback = str2func(varargin{1});
	end
	
	if nargout
		[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
	else
		gui_mainfcn(gui_State, varargin{:});
	end
	% End initialization code - DO NOT EDIT
	
end
%% these are functions that don't do anything, but are required by matlab
function ptrListCellType_CreateFcn(hObject, eventdata, handles), end %#ok<DEFNU>
function ptrListSelectImage_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrListSelectPixelResponsiveness_CreateFcn(hObject, eventdata, handles), end %#ok<DEFNU>
function ptrEditMagnification_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrEditZoom_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrEditSelect_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrListPresence_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrListRespType_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrListPresence_Callback(hObject, eventdata, handles), end %#ok<DEFNU>
function ptrListRespType_Callback(hObject, eventdata, handles), end %#ok<DEFNU>
function ptrEditSubSelect_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrListSelectDrawType_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ptrListBoundaryType_CreateFcn(hObject, eventdata, handles),end %#ok<DEFNU>
function ButtonGroupAnnotations_CreateFcn(hObject, eventdata, handles),set(hObject,'SelectionChangeFcn','runDetectCells(''ChangeAnnotations'')');end %#ok<DEFNU>
function ptrListBoundaryType_Callback(hObject, eventdata, handles),end %#ok<DEFNU>

%% opening function; initializes output
function runDetectCells_OpeningFcn(hObject, eventdata, handles, varargin)
	%opening actions
	
	%set closing function
	set(hObject,'DeleteFcn','DC_DeleteFcn')
	
	% set rainbow logo
	I = imread('DC_rainbow.tif');
	%I2 = imresize(I, [100 NaN]);
	axes(handles.ptrAxesLogo); %#ok<MAXES>
	imshow(I);
	
	% set default output
	handles.output = hObject;
	guidata(hObject, handles);
	
end
%% defines output variables
function varargout = runDetectCells_OutputFcn(hObject, eventdata, handles)
	%output
	varargout{1} = handles.output;
end
%% exit/save functions
function ptrButtonSave_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%saves data
	
	% make structures global
	global sDC;
	global sRec;
	global sFig;
	
	%save data
	if ~isempty(sDC)
		%set msg
		cellText{1} = 'Saving data... Please wait';
		DC_updateTextInformation(cellText);
		
		%get filename
		strRecFile = sDC.strRecFile;
		strRecPath = sDC.strRecPath;
		if strcmp(getFlankedBy(strRecFile,'_','.mat'),'CD')
			%filename already has cell detection append
			strSaveFile = strRecFile(1:end-4);
		elseif strcmp(getFlankedBy(strRecFile,'_','.mat'),'prepro')
			%filename ends with _prepro.mat
			strSaveFile = [strRecFile(1:end-11) '_CD.mat'];
		else
			strSaveFile = [strRecFile(1:end-4) '_CD.mat'];
		end
		
		%assign to structure
		sRec.sDC = sDC;
		
		%save data
		save([strRecPath strSaveFile], 'sRec' );
		
		%check if figure is still there
		try
			vGet = get(sFig.ptrWindowHandle);
		catch %#ok<CTCH>
			vGet = [];
		end
		
		%save image of figure
		if ~isempty(vGet)
			export_fig([strRecPath 'average' filesep strRecFile(1:end-4) '.tif'],sFig.ptrWindowHandle)
			strFigSave = 'and figure ';
		else
			strFigSave = '';
		end
		
		%set msg
		cellText{1} = sprintf('Saved data %sto:',strFigSave);
		cellText{2} = [strRecPath strSaveFile];
		DC_updateTextInformation(cellText);
	end
end
function ptrButtonSaveExit_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%saves data and exits GUI
	
	% make structures global
	global sDC;
	global sRec;
	global sFig;
	
	%save data
	if ~isempty(sDC)
		%set msg
		cellText{1} = 'Saving data... Please wait';
		DC_updateTextInformation(cellText);
		
		%get filename
		strRecFile = sDC.strRecFile;
		strRecPath = sDC.strRecPath;
		if strcmp(getFlankedBy(strRecFile,'_','.mat'),'CD')
			%filename already has cell detection append
			strSaveFile = strRecFile(1:end-4);
		elseif strcmp(getFlankedBy(strRecFile,'_','.mat'),'prepro')
			%filename ends with _prepro.mat
			strSaveFile = [strRecFile(1:end-11) '_CD.mat'];
		else
			strSaveFile = [strRecFile(1:end-4) '_CD.mat'];
		end
		
		%assign to structure
		sRec.sDC = sDC;
		
		%save data
		save([strRecPath strSaveFile], 'sRec' );
		
		%set msg
		cellText{1} = 'Saved data to:';
		cellText{2} = [strRecPath strSaveFile];
		DC_updateTextInformation(cellText);
	end
	
	%check if figure is still there
	try
		vGet = get(sFig.ptrWindowHandle);
	catch %#ok<CTCH>
		vGet = [];
	end
	
	%close figure
	if ~isempty(vGet)
		%save image of figure
		
		%close figure
		close(sFig.ptrWindowHandle);
	end
	
	% unglobal sRec & CDG
	clear sRec;
	clear sFig;
	clear sDC;
	
	% close main gui
	close;
end
function ptrButtonNoSaveExit_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	% this function exits the cell detection gui without saving
	
	% make structures global
	global sDC; %#ok<NUSED>
	global sRec; %#ok<NUSED>
	global sFig;
	
	%check if figure is still there
	try
		vGet = get(sFig.ptrWindowHandle);
	catch %#ok<CTCH>
		vGet = [];
	end
	
	%close figure
	if ~isempty(vGet)
		close(sFig.ptrWindowHandle);
	end
	
	% unglobal sRec & CDG
	clear sRec;
	clear sFig;
	clear sDC;
	
	% close main gui
	close;
end
%% detect object boundaries
function ptrButtonDrawOutline_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%waits for user to click anywhere on window to zoom in; then creates
	%zoomed-in figure allowing for more precise selection of boundaries
	
	%lock GUI
	DC_lock(handles);
	
	%start drawing function
	DC_DrawOutline;
	
	%redraw cells
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
%% add object center
function ptrButtonAddCenter_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%adds center point for later boundary detection
	
	%lock GUI
	DC_lock(handles);
	
	%add center point
	DC_AddCenter;
	
	%msg
	DC_updateTextInformation;
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
%% remove selected object
function ptrButtonDelSelected_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get global
	global sFig;
	global sDC;
	
	%check if anything is selected
	if isempty(sDC), return;
	elseif ~isempty(sFig.vecSelectedObjects)
		
		%lock GUI
		DC_lock(handles);
		
		%get selected objects
		vecSelected = sFig.vecSelectedObjects;
		sFig.vecSelectedObjects = [];
		
		%remove selected objects
		DC_removeObject(vecSelected);
		
		%redraw
		DC_redraw(0);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
%% remove last object
function ptrButtonDelLast_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get global
	global sDC
	global sFig
	
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
	
	%lock GUI
	DC_lock(handles);
	
	%unselect all objects
	sFig.vecSelectedObjects = [];
	
	%get last object
	intLastObject = numel(sDC.ROI);
	
	%remove objects
	DC_removeObject(intLastObject);
	
	%redraw
	DC_redraw(1);
	
	%unlock GUI
	DC_unlock(handles);
end

%% detect boundaries from centers
function ptrButtonDetectBoundaries_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%detects boundaries
	
	%check if data is loaded
	global sFig;
	try
		%get current image
		intImSelected = get(sFig.ptrListSelectImage,'Value');
	catch %#ok<CTCH>
		rethrow(lasterror)
		return;
	end
	
	%lock GUI
	DC_lock(handles);
	
	%set msg
	cellText{1} = 'Calculating boundaries... Please wait';
	DC_updateTextInformation(cellText);
	
	%detect objects
	[intDetectObjects,indDelete] = DC_detectObjects;
	
	%set msg
	cellText{1} = ['Detected ' num2str(intDetectObjects) ' objects'];
	vecDelete = find(indDelete);
	if ~isempty(vecDelete)
		cellText{2} = '';
		cellText{3} = 'Unable to disentangle following objects:';
		cellText{4} = sprintf('%d; ',vecDelete);
		cellText{5} = '';
		cellText{6} = 'Please check these objects manually';
	end
	DC_updateTextInformation(cellText);
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
%% select which image to display as background
function ptrListSelectImage_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%selected image is automatically queried by drawing function; so no
	%other action is required other than redrawing
	
	%lock GUI
	DC_lock(handles);
	
	%redraw
	DC_redraw(1);
	
	%unlock GUI
	DC_unlock(handles);
end
%% set zoom level
function ptrEditZoom_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%no action required
end
%% set zoom level
function ptrEditMagnification_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%magnification is automatically queried by drawing function; so no
	%other action is required other than redrawing
	
	%lock GUI
	DC_lock(handles);
	
	%redraw
	DC_redraw(2);
	
	%unlock GUI
	DC_unlock(handles);
end
%% this function initializes everything
function ptrButtonLoadRecording_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%This function lets the user select a file and then loads the
	%pre-processed information
	
	%get globals
	global sRec;
	global sFig;
	global sDC;
	
	%lock GUI
	DC_lock(handles);
	
	%switch path
	try
		oldPath = cd('D:\Data\Processed\imagingdata');
	catch
		oldPath = cd();
	end
	
	%get file
	[strRecFile, strRecPath] = uigetfile('*.mat', 'Select Recording file');
	
	%back to old path
	cd(oldPath);
	
	%check whether a file has been selected, otherwise, do nothing
	if ischar(strRecFile)
		%clear old data
		sRec = []; %#ok<NASGU>
		sFig = []; %#ok<NASGU>
		sDC = []; %#ok<NASGU>
		
		%load structure
		sLoad = load([strRecPath strRecFile]);
		if ~isfield(sLoad,'sRec')
			error([mfilename ':NoRecordingStructureInFile'],'Selected file does not contain preprocessing data [%s%s]',strRecPath,strRecFile);
		end
		sRec = sLoad.sRec;
		clear sLoad;
		
		%update figure
		set(handles.ptrTextRecording, 'string', sRec.sProcLib.strRecording);
		set(handles.ptrTextSession, 'string', sRec.strSession);
		
		%get figure location
		intCaCh = sRec.sProcLib.ch; %channel that contains calcium data
		strPath = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording];
		sDir = dir(strPath);
		if isempty(sDir),strPath = strRecPath;end
		strImage = [strPath filesep 'average' filesep 'OverlayProc.tif'];
		
		%load image
		imProc = im2double(imread(strImage));
		if intCaCh ~= 2
			imNew = imProc;
			imNewCh2 = imProc(:,:,intCaCh);
			imNew(:,:,intCaCh) = imNew(:,:,2);
			imNew(:,:,2) = imNewCh2;
			imProc = imNew;
			clear imNew;
		end
		
		%populate cell detection structure
		if isfield(sRec,'sDC')
			sDC = sRec.sDC;
		else
			sDC = DC_populateStructure();
		end
		sDC.strRecFile = strRecFile;
		sDC.strRecPath = strRecPath;
		sDC.metaData.strProcessedPath = strRecPath;
		sRec.sMD.strMasterDir(1) = strRecPath(1);
		
		%populate figure structure
		sFig = DC_populateFigure(handles);
		sFig.imProc = imProc;
		
		%put pixel responsiveness types into list
		if isfield(sRec,'sPixResp') && isfield(sRec.sPixResp,'matPixelSelectivity') && isfield(sRec.sPixResp,'matPixelMaxResponsiveness')
			%for backward compatibility
			sDC.metaData.cellPixRespType{1} = 'None';
			sDC.metaData.cellPixRespType{2} = 'Selectivity';
			sDC.metaData.cellPixRespType{3} = 'Activation';
			
			%set pixel resp types to those defined in DC_populateStructure
			set(sFig.ptrListSelectPixelResponsiveness,'String',sDC.metaData.cellPixRespType)
			
			%add maps to list
			sRec.sPixResp.cellExtraIm{1} = zeros(size(sRec.sPixResp.matPixelSelectivity),class(sRec.sPixResp.matPixelSelectivity));
			sRec.sPixResp.cellExtraIm{2} = sRec.sPixResp.matPixelSelectivity;
			sRec.sPixResp.cellExtraIm{3} = sRec.sPixResp.matPixelMaxResponsiveness;
			
			% put in list
			intType = get(sFig.ptrListSelectPixelResponsiveness, 'Value');
			
			%put channel into image
			matPixResp = sRec.sPixResp.cellExtraIm{intType};
			matNorm = imnorm(matPixResp);
			if min(matNorm(:)) == max(matNorm(:)),matNorm=zeros(size(matNorm),class(matNorm));end
			sFig.imProc(:,:,3) = matNorm;
		end
		
		%put images into list
		[cellIm,cellImName] = DC_createImageList(sFig.imProc);
		sFig.cellIm = cellIm;
		set(sFig.ptrListSelectImage,'String',cellImName);
		sFig.intImSelected = get(sFig.ptrListSelectImage,'Value');
		
		%populate presence list
		%for backward compatibility
		if ~isfield(sDC.metaData,'cellPresence')
			sDC.metaData.cellPresence{1} = 'present';
			sDC.metaData.cellPresence{2} = 'absent';
			sDC.metaData.cellPresence{3} = 'filled';
		end
		set(sFig.ptrListPresence,'String',sDC.metaData.cellPresence);%set designation types to those defined in DC_populateStructure
		
		%populate responsiveness list
		%for backward compatibility
		if ~isfield(sDC.metaData,'cellRespType')
			sDC.metaData.cellRespType{1} = 'tun+resp';
			sDC.metaData.cellRespType{2} = 'tuned';
			sDC.metaData.cellRespType{3} = 'responsive';
			sDC.metaData.cellRespType{4} = 'silent';
		end
		set(sFig.ptrListRespType,'String',sDC.metaData.cellRespType);%set designation types to those defined in DC_populateStructure
		set(sFig.ptrListRespType,'Value',4); %set selected type to be silent by default
		
		%populate draw type list
		%for backward compatibility
		sDC.metaData.cellDrawType{1} = 'Border';
		sDC.metaData.cellDrawType{2} = 'Centroid';
		set(sFig.ptrListSelectDrawType,'String',sDC.metaData.cellDrawType);%set designation types to those defined in DC_populateStructure
		set(sFig.ptrListSelectDrawType,'Value',1); %set selected type to be silent by default
		
		%set cell types to those defined in DC_populateStructure
		set(handles.ptrListCellType,'String',sDC.metaData.cellType)
		
		%set to hide annotations
		set(sFig.ptrButtonRadioAnnotationsHide,'Value',1);
		
		%populate boundary type list
		%for backward compatibility
		sDC.metaData.cellBoundaryType{1} = 'OGB';
		sDC.metaData.cellBoundaryType{2} = 'GCaMP';
		set(sFig.ptrListBoundaryType,'String',sDC.metaData.cellBoundaryType);%set designation types to those defined in DC_populateStructure
		set(sFig.ptrListBoundaryType,'Value',1); %set selected type to be silent by default
		
		%set msg
		cellText{1} = sprintf('Loaded recording %s from %s',strRecFile,strRecPath);
		DC_updateTextInformation(cellText);
		
		%draw image
		DC_redraw(2);
	end
	
	%unlock GUI
	DC_unlock(handles);
end

%% add channel to background image (such as 910nm TdTomato reference image)
function ptrButtonAddImChannel_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%add image as channel to stack
	
	%get globals
	global sFig;
	global sDC;
	global sRec;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig) || isempty(sRec)
		return;
	else
		try
			intType = get(sFig.ptrListCellType, 'Value');
		catch
			return;
		end
	end
	
	%lock GUI
	DC_lock(handles);
	
	%switch path
	try
		oldPath = cd(sDC.metaData.strRawPath);
	catch
		oldPath = cd();
	end
	%get file
	[strImFile, strImPath] = uigetfile('*.*', 'Select image file');
	
	%back to old path
	cd(oldPath);
	
	%check whether a file has been selected, otherwise, do nothing
	if ischar(strImFile)
		%set msg
		cellText{1} = sprintf('Loaded image %s from %s',strImFile,strImPath);
		DC_updateTextInformation(cellText);
		
		%load image
		intMaxType = length(sDC.metaData.cellPixRespType);
		imCh1 = im2double(imread([strImPath strImFile]));
		matNorm = imnorm(imCh1);
		sRec.sPixResp.cellExtraIm{intMaxType+1} = matNorm;
		
		%set pixel resp types to those defined in DC_populateStructure
		sDC.metaData.cellPixRespType{intMaxType+1} = sprintf('Extra map %d',intMaxType+1);
		set(sFig.ptrListSelectPixelResponsiveness,'String',sDC.metaData.cellPixRespType)
		
		% set to current value
		set(sFig.ptrListSelectPixelResponsiveness, 'Value',intMaxType+1);
		
		%put channel into image
		sFig.imProc(:,:,3) = matNorm;
		
		%put images into list
		[cellIm,cellImName] = DC_createImageList(sFig.imProc);
		sFig.cellIm = cellIm;
		set(sFig.ptrListSelectImage,'String',cellImName);
		sFig.intImSelected = get(sFig.ptrListSelectImage,'Value');
		
		%draw image
		DC_redraw(1);
		
		%unlock GUI
		DC_unlock(handles);
	end
	
	%unlock GUI
	DC_unlock(handles);
end
%% add ROIs only
function ptrButtonLoadROI_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%globals
	global sDC
	global sFig
	
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
	
	%lock GUI
	DC_lock(handles);
	
	%switch path
	try
		oldPath = cd(sDC.metaData.strProcessedPath);
	catch
		oldPath = cd();
	end
	
	%get file
	[strRecFile, strRecPath] = uigetfile('*.mat', 'Select Recording file');
	
	%back to old path
	cd(oldPath);
	
	%check whether a file has been selected, otherwise, do nothing
	if ischar(strRecFile)
		%load structure
		sLoad = load([strRecPath strRecFile]);
		if ~isfield(sLoad,'sRec')
			%incorrect; send msg
			warning([mfilename ':NoRecordingStructureInFile'],['Selected file does not contain preprocessing data [%s%s]'],strRecPath,strRecFile);
			
			cellText{1} = sprintf('Selected file does not contain preprocessing data');
			cellText{2} = [''] ;
			cellText{3} = sprintf('[%s%s]',strRecPath,strRecFile) ;
			DC_updateTextInformation(cellText);
			
			%unlock GUI
			DC_unlock(handles);
			
			%return
			return
		end
		sRec = sLoad.sRec;
		clear sLoad;
		
		%load sDC
		if isfield(sRec,'sDC')
			%assign ROIs & metadata
			sDC.ROI = sRec.sDC.ROI;
			sDC.metaData = sRec.sDC.metaData;
			
			%calculate initial guess of ROI displacement
			%get figure location
			intCaCh = sRec.sProcLib.ch; %channel that contains calcium data
			strPath = [sRec.sMD.strMasterDir sRec.sMD.strImgTarget sRec.strSession filesep sRec.sProcLib.strRecording];
			sDir = dir(strPath);
			if isempty(sDir),strPath = strRecPath;end
			strImage = [strPath filesep 'average' filesep 'OverlayProc.tif'];
			
			
			%load image
			imProc = im2double(imread(strImage));
			if intCaCh ~= 2
				imNew = imProc;
				imNewCh2 = imProc(:,:,intCaCh);
				imNew(:,:,intCaCh) = imNew(:,:,2);
				imNew(:,:,2) = imNewCh2;
				imROI = imNew;
				clear imNew;
			else
				imROI = imProc;
			end
			sFig.imROI = imROI;
			
			%calculate required registration
			[output matImReg] = dftregistration( fft2(sFig.imROI(:,:,2)), fft2(sFig.imProc(:,:,2)), 100);
			
			%put output into log
			sDC.metaData.intROIDisplacementX = output(4); %x/y inverted
			sDC.metaData.intROIDisplacementY = output(3); %x/y inverted
			
			%set msg
			cellText{1} = sprintf('Loaded ROI information from %s%s',sRec.strSession,sRec.sProcLib.strRecording);
			cellText{2} = [''] ;
			cellText{3} = sprintf('Probable ROI displacement; x: %.1f; y=%.1f',sDC.metaData.intROIDisplacementX,sDC.metaData.intROIDisplacementY) ;
			DC_updateTextInformation(cellText);
			
			%redraw
			DC_redraw(1);
		end
	end
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonMoveROI_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%moves ROIs
	
	%lock GUI
	%DC_lock(handles);
	
	%move ROIs
	DC_moveROIs;
	
	%unlock GUI
	%DC_unlock(handles);
end
function ptrButtonRecalculateBoundaries_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%globals
	global sDC;
	global sFig;
	
	%check if data has been loaded
	if isempty(sDC) || isempty(sFig)
		return;
	else
		try
			%get current image
			intImSelected = get(sFig.ptrListSelectImage,'Value');
		catch %#ok<CTCH>
			return;
		end
	end
	
	%lock GUI
	DC_lock(handles);
	
	%set msg
	cellText{1} = 'Recalculating boundaries... Please wait';
	DC_updateTextInformation(cellText);
	
	%check selection
	if isempty(sFig.vecSelectedObjects)
		vecObjects = 1:numel(sDC.ROI);
	else
		vecObjects = sFig.vecSelectedObjects;
	end
				
	%loop through selected objects / or through all objects
	for intObject = vecObjects
		%get presence
		intPresence = sDC.ROI(intObject).intPresence;
		strPresence = sDC.metaData.cellPresence{intPresence};
		
		%skip if absent
		if strcmp(strPresence,'absent')
			continue;
		end
		
		% set drawn flag
		sFig.sObject(intObject).drawn = 0;
		
		%remove mask
		sDC.ROI(intObject).matMask = [];
		
		%remove old drawing
		if isfield(sFig.sObject(intObject).handles,'lines')
			for p = 1:length(sFig.sObject(intObject).handles.lines)
				delete(sFig.sObject(intObject).handles.lines(p));
			end
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
			delete(sFig.sObject(intObject).handles.marker);
		end
		if isfield(sFig.sObject(intObject).handles,'text') && ~isempty(sFig.sObject(intObject).handles.text)
			delete(sFig.sObject(intObject).handles.text);
			sFig.sObject(intObject).handles.text = [];
		end
	end
	
	%detect objects
	[intDetectObjects,indDelete] = DC_detectObjects(vecObjects);
	
	%set msg
	cellText{1} = ['Detected ' num2str(intDetectObjects) ' objects'];
	vecDelete = find(indDelete);
	if ~isempty(vecDelete)
		cellText{2} = '';
		cellText{3} = 'Unable to disentangle following objects:';
		cellText{4} = sprintf('%d; ',vecDelete);
		cellText{5} = '';
		cellText{6} = 'Please check these objects manually';
	end
	DC_updateTextInformation(cellText);
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonCheckDisplacements_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%make distance-difference matrix per neuronal pairs for both recordings
	
	%globals
	global sDC
	global sFig
	
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
	
	%lock GUI
	DC_lock(handles);
	
	%switch path
	try
		oldPath = cd(sDC.metaData.strProcessedPath);
	catch
		oldPath = cd();
	end
	
	%get file
	[strRecFile, strRecPath] = uigetfile('*.mat', 'Select Recording file');
	
	%back to old path
	cd(oldPath);
	
	%check whether a file has been selected, otherwise, do nothing
	if ischar(strRecFile)
		%load structure
		sLoad = load([strRecPath strRecFile]);
		if ~isfield(sLoad,'sRec')
			error([mfilename ':NoRecordingStructureInFile'],['Selected file does not contain preprocessing data [' strRecPath strRecFile ']']);
		end
		sRecRef = sLoad.sRec;
		clear sLoad;
		
		%load sDC
		if isfield(sRecRef,'sDC') && ~isempty(sDC)
			sDCRef = sRecRef.sDC;
			
			intObjectsRef = numel(sDCRef.ROI);
			intObjects = numel(sDC.ROI);
			
			if intObjects ~= intObjectsRef
				%set msg
				cellText{1} = sprintf('ROI shift check output; comparing to %s',strRecFile);
				cellText{2} = '' ;
				cellText{3} = 'ERROR:';
				cellText{4} = '# of objects is not the same!';
				cellText{5} = sprintf('This recording; %d;',intObjects);
				cellText{6} = sprintf('Reference recording; %d;',intObjectsRef);
				DC_updateTextInformation(cellText)
			else
				%check distances
				matDistances = zeros(intObjects);
				matDistancesRef = zeros(intObjects);
				
				for intObject1=1:intObjects
					for intObject2=(intObject1+1):intObjects
						%get distance this recording
						dblCell1LocX = sDC.ROI(intObject1).intCenterX;
						dblCell1LocY = sDC.ROI(intObject1).intCenterY;
						dblCell2LocX = sDC.ROI(intObject2).intCenterX;
						dblCell2LocY = sDC.ROI(intObject2).intCenterY;
						xDist = abs(dblCell1LocX - dblCell2LocX);
						yDist = abs(dblCell1LocY - dblCell2LocY);
						dblCellDist = sqrt(xDist^2+yDist^2);
						
						matDistances(intObject1,intObject2) = dblCellDist;
						matDistances(intObject2,intObject1) = dblCellDist;
						
						%get distance other recording
						dblCell1LocX = sDCRef.ROI(intObject1).intCenterX;
						dblCell1LocY = sDCRef.ROI(intObject1).intCenterY;
						dblCell2LocX = sDCRef.ROI(intObject2).intCenterX;
						dblCell2LocY = sDCRef.ROI(intObject2).intCenterY;
						xDist = abs(dblCell1LocX - dblCell2LocX);
						yDist = abs(dblCell1LocY - dblCell2LocY);
						dblCellDist = sqrt(xDist^2+yDist^2);
						
						matDistancesRef(intObject1,intObject2) = dblCellDist;
						matDistancesRef(intObject2,intObject1) = dblCellDist;
					end
				end
				matDistancesDiff = abs(matDistances - matDistancesRef);
				
				%output
				[dblMaxDiff,intIndex] = max(matDistancesDiff(:));
				[intY,intX] = ind2sub(size(matDistancesDiff),intIndex);
				
				%make figure
				h=figure;
				imagesc(matDistancesDiff);
				colorbar;
				title(sprintf('Maximum centroid shift is %.1f pixels [object %d]',dblMaxDiff,intY));
				xlabel('Object 1');
				ylabel('Object 2');
				
				%set msg
				cellText{1} = sprintf('ROI shift check output; comparing to %s',strRecFile);
				cellText{2} = '' ;
				cellText{3} = sprintf('Maximum centroid shift is %.1f pixels [object %d]',dblMaxDiff,intY);
				DC_updateTextInformation(cellText)
			end
		end
	end
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonSetToCurrentType_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%sets selected neurons to current type
	
	%get structures
	global sDC;
	global sFig;
	
	%check if anything is selected
	if isempty(sDC), return; end
	if isempty(sFig.vecSelectedObjects), return; end
	
	%lock GUI
	DC_lock(handles);
	
	%set selected objects to current type
	for intObject = 1:numel(sDC.ROI)
		if ismember(intObject,sFig.vecSelectedObjects)
			sDC.ROI(intObject).intType = get(sFig.ptrListCellType, 'Value');
		end
	end
	
	%update text
	DC_updateTextInformation;
	
	%redraw
	DC_redraw(1);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonRedraw_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%globals
	global sFig;
	global sDC;
	
	%check for data
	if isempty(sDC), return;
	elseif length(sFig.vecSelectedObjects) == 1
		%unlock GUI
		DC_lock(handles);
		
		%remove drawing of selected object
		intObject = sFig.vecSelectedObjects(1);
		if isfield(sDC.ROI,'matMask') && ~isempty(sDC.ROI(intObject).matMask) && sFig.sObject(intObject).drawn
			% remove old drawing and set drawn flag
			if isfield(sFig.sObject(intObject).handles,'lines')
				for p = 1:length(sFig.sObject(intObject).handles.lines)
					try delete(sFig.sObject(intObject).handles.lines(p));catch, end
				end
				sFig.sObject(intObject).handles.lines = [];
			else
				p = [];
			end
			if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX) && sFig.sObject(intObject).drawn == 1
				try delete(sFig.sObject(intObject).handles.marker);catch, end
			end
			sFig.sObject(intObject).drawn = 0;
		end
		
		%redraw
		DC_DrawOutline(intObject);
		
		%msg
		DC_updateTextInformation;
		
		%redraw
		DC_redraw(0);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
function ptrButtonUnselectAll_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%globals
	global sFig;
	global sDC;
	
	%check for data
	if isempty(sDC), return;
	elseif ~isempty(sFig.vecSelectedObjects)
		%lock GUI
		DC_lock(handles);
		
		%set selection & update text
		sFig.vecSelectedObjects = [];
		DC_updateTextInformation;
		
		%redraw
		DC_redraw(1);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
function ptrEditSelect_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%select a single object
	
	%get globals
	global sDC;
	global sFig;
	
	%check for data
	if isempty(sDC), return;end
	
	%get data
	intSelect = round(str2double(get(hObject,'String')));
	intObjects = numel(sDC.ROI);
	vecObjects = 1:intObjects;
	
	if ismember(intSelect,vecObjects)
		%lock GUI
		DC_lock(handles);
		
		%remove all selected objects
		DC_unselectAll();
		
		%set selection & update text
		sFig.vecSelectedObjects = intSelect;
		DC_updateTextInformation;
		
		% remove old drawing and set drawn flag
		if isfield(sFig.sObject(intSelect).handles,'lines')
			for p = 1:length(sFig.sObject(intSelect).handles.lines)
				delete(sFig.sObject(intSelect).handles.lines(p));
			end
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intSelect),'intCenterX') && ~isempty(sDC.ROI(intSelect).intCenterX)
			delete(sFig.sObject(intSelect).handles.marker);
		end
		if isfield(sFig.sObject(intSelect).handles,'text') && ~isempty(sFig.sObject(intSelect).handles.text)
			delete(sFig.sObject(intSelect).handles.text);
			sFig.sObject(intSelect).handles.text = [];
		end
		sFig.sObject(intSelect).drawn = 0;
		
		%redraw
		DC_redraw(0);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
function ptrEditSubSelect_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%select a single object by subtype nr
	
	%get globals
	global sDC;
	global sFig;
	
	%check for data
	if isempty(sDC), return;end
	
	%get data
	intSubType = get(sFig.ptrListCellType, 'Value');
	intSubSelectNr = round(str2double(get(hObject,'String')));
	intSubNr = 0;
	intObjects = numel(sDC.ROI);
	boolRunning = true;
	boolFound = false;
	intObject = 0;
	while boolRunning
		intObject = intObject + 1;
		if intObject > intObjects
			boolRunning = false;
		elseif sDC.ROI(intObject).intType == intSubType
			intSubNr = intSubNr + 1;
			if intSubNr == intSubSelectNr
				boolRunning = false;
				boolFound = true;
			end
		end
	end
	
	if boolFound
		%lock GUI
		DC_lock(handles);
		
		%remove all selected objects
		DC_unselectAll();
		
		%set selection & update text
		sFig.vecSelectedObjects = intObject;
		DC_updateTextInformation;
		
		% remove old drawing and set drawn flag
		if isfield(sFig.sObject(intObject).handles,'lines')
			for p = 1:length(sFig.sObject(intObject).handles.lines)
				delete(sFig.sObject(intObject).handles.lines(p));
			end
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
			delete(sFig.sObject(intObject).handles.marker);
		end
		if isfield(sFig.sObject(intObject).handles,'text') && ~isempty(sFig.sObject(intObject).handles.text)
			delete(sFig.sObject(intObject).handles.text);
			sFig.sObject(intObject).handles.text = [];
		end
		sFig.sObject(intObject).drawn = 0;
		
		%redraw
		DC_redraw(0);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
function ptrButtonRealignROIs_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%globals
	global sFig;
	global sDC;
	
	%check for data
	if isempty(sDC), return;
	else
		%lock GUI
		DC_lock(handles);
		
		%check reference image, otherwise ask for source
		if ~isfield(sFig,'imROI') || isempty(sFig.imROI)
			%set msg
			clear cellText;
			cellText{1} = sprintf('No reference image loaded, please use the ''Load ROIs'' button first');
			DC_updateTextInformation(cellText);
			
			%unlock GUI
			DC_unlock(handles);
			
			%return
			return
		end
		
		%set msg
		clear cellText;
		cellText{1} = sprintf('Recursive locally affine subfield reregistration algorithm is running');
		cellText{2} = '' ;
		cellText{3} = sprintf('Calculating ROI shifts... Please wait');
		DC_updateTextInformation(cellText);
		
		%calculate required registration matrix
		matDataPoints = doRegisterSubFields(mean(sFig.imROI(:,:,1:2),3),mean(sFig.imProc(:,:,1:2),3));
		[matXi,matYi] = meshgrid(1:512);
		[matTranslationsY,matTranslationsX] = doInterpolateFromList(matDataPoints,matXi,matYi);
		
		figure
		subplot(2,2,1),imagesc(matTranslationsX),colorbar, title('X pixel shift')
		subplot(2,2,2),imagesc(matTranslationsY),colorbar, title('Y pixel shift')
		subplot(2,2,3),imagesc(sqrt(matTranslationsY.*matTranslationsY + matTranslationsX.*matTranslationsX)),colorbar, title('Total pixel shift')
		
		
		%set msg
		clear cellText;
		cellText{1} = sprintf('Recursive locally affine subfield reregistration algorithm is finished');
		cellText{2} = '' ;
		cellText{3} = sprintf('Applying translation vector map to all ROIs... Please wait');
		DC_updateTextInformation(cellText);
		
		%apply translations
		for intObject=1:numel(sDC.ROI)
			%get mean translation data
			matTranslateX = matTranslationsX(logical(sDC.ROI(intObject).matMask));
			matTranslateY = matTranslationsY(logical(sDC.ROI(intObject).matMask));
			intTranslateX = round(mean(matTranslateX(:)));
			intTranslateY = round(mean(matTranslateY(:)));
			
			sDC.ROI(intObject).intCenterX = sDC.ROI(intObject).intCenterX - intTranslateX;
			sDC.ROI(intObject).intCenterY = sDC.ROI(intObject).intCenterY - intTranslateY;
			
			%sDC.ROI(intObject).matPerimeter = [sDC.ROI(intObject).matPerimeter(:,1)-intTranslateX sDC.ROI(intObject).matPerimeter(:,2)-intTranslateY];
			sDC.ROI(intObject).matMask = circshift(sDC.ROI(intObject).matMask,[-intTranslateY -intTranslateX]);
			
			% remove old drawing and set drawn flag
			for p = 1:length(sFig.sObject(intObject).handles.lines)
				delete(sFig.sObject(intObject).handles.lines(p));
			end
			if isempty(p) && isfield(sDC.ROI(intObject),'intCenterX') && ~isempty(sDC.ROI(intObject).intCenterX)
				delete(sFig.sObject(intObject).handles.marker);
			end
			if isfield(sFig.sObject(intObject).handles,'text') && ~isempty(sFig.sObject(intObject).handles.text)
				delete(sFig.sObject(intObject).handles.text);
			end
			sFig.sObject(intObject).drawn = 0;
			sFig.sObject(intObject).handles.marker = [];
			sFig.sObject(intObject).handles.lines = [];
			sFig.sObject(intObject).handles.text = [];
		end
		%redraw
		DC_redraw(0);
		
		%set msg
		clear cellText;
		cellText{1} = sprintf('ROIs have been relocated, please check for inconsistencies');
		DC_updateTextInformation(cellText);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
function ptrListCellType_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%globals
	global sFig;
	global sDC;
	
	%check for data
	if isempty(sDC), return;
	else
		set(sFig.ptrTextTotalType, 'String', DC_getTotalType());
	end
end
function ptrListSelectPixelResponsiveness_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%put selected image type in figure list; then redraw
	
	%globals
	global sFig;
	global sDC;
	global sRec;
	
	%check for data
	if isempty(sDC) || ~isfield(sRec,'sPixResp') || ~isfield(sRec.sPixResp,'matPixelSelectivity') || ~isfield(sRec.sPixResp,'matPixelMaxResponsiveness'), return;
	else
		%lock GUI
		DC_lock(handles);
		
		% get selection
		intType = get(sFig.ptrListSelectPixelResponsiveness, 'Value');
		
		%put channel into image
		matPixResp = sRec.sPixResp.cellExtraIm{intType};
		matNorm = imnorm(matPixResp);
		matNorm(isnan(matNorm)) = nanmean(matNorm(:));
		if min(matNorm(:)) == max(matNorm(:)),matNorm=zeros(size(matNorm),class(matNorm));end
		sFig.imProc(:,:,3) = matNorm;
		
		%put images into list
		[cellIm,cellImName] = DC_createImageList(sFig.imProc);
		sFig.cellIm = cellIm;
		set(sFig.ptrListSelectImage,'String',cellImName);
		sFig.intImSelected = get(sFig.ptrListSelectImage,'Value');
		
		%draw image
		DC_redraw(1);
		
		%unlock GUI
		DC_unlock(handles);
	end
end
function ptrButtonSetPresence_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%sets selected neurons to current presence designation
	
	%get structures
	global sDC;
	global sFig;
	
	%check if anything is selected
	if isempty(sDC), return; end
	if isempty(sFig.vecSelectedObjects), return; end
	
	%lock GUI
	DC_lock(handles);
	
	%set selected objects to current type
	for intObject = 1:numel(sDC.ROI)
		if ismember(intObject,sFig.vecSelectedObjects)
			sDC.ROI(intObject).intPresence = get(sFig.ptrListPresence, 'Value');
		end
	end
	
	%update text
	DC_updateTextInformation;
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonSetRespType_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%sets selected neurons to current response type
	
	%get structures
	global sDC;
	global sFig;
	
	%check if anything is selected
	if isempty(sDC), return; end
	if isempty(sFig.vecSelectedObjects), return; end
	
	%lock GUI
	DC_lock(handles);
	
	%set selected objects to current type
	for intObject = 1:numel(sDC.ROI)
		if ismember(intObject,sFig.vecSelectedObjects)
			sDC.ROI(intObject).intRespType = get(sFig.ptrListRespType, 'Value');
		end
	end
	
	%update text
	DC_updateTextInformation;
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrListSelectDrawType_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%lock GUI
	DC_lock(handles);
	
	%draw image
	DC_redraw(1);
	
	%unlock GUI
	DC_unlock(handles);
end
function ChangeAnnotations %#ok<DEFNU>
	%get structures
	global sFig
	
	%check initialization
	if isfield(sFig,'ButtonGroupAnnotations')
		%change annotation type
		sFig.strAnnotations = get(get(sFig.ButtonGroupAnnotations,'SelectedObject'),'String');
		
		%lock GUI
		DC_lock(sFig);
		
		%draw image
		DC_redraw(1);
		
		%unlock GUI
		DC_unlock(sFig);
	end
	
	%SelectedObject
end

function ptrButtonSelectAllOfType_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get structures
	global sDC;
	global sFig;
	
	%check for data
	if isempty(sDC), return;end
	
	%clear selected objects
	DC_unselectAll();
	
	%get data
	intSubType = get(sFig.ptrListCellType, 'Value');
	intObjects = numel(sDC.ROI);
	sFig.vecSelectedObjects = nan(1,intObjects); %over-pre-allocate
	intSelectionObject = 0;
	for intObject=1:intObjects
		if sDC.ROI(intObject).intType == intSubType
			intSelectionObject = intSelectionObject + 1;
			sFig.vecSelectedObjects(intSelectionObject) = intObject;
		end
	end
	
	%remove trailing nans
	sFig.vecSelectedObjects = sFig.vecSelectedObjects(1:(find(isnan(sFig.vecSelectedObjects),1,'first')-1));
	
	%lock GUI
	DC_lock(handles);
	
	%set selection & update text
	DC_updateTextInformation;
	
	% remove old drawings and set drawn flags
	for intSelectedObject=sFig.vecSelectedObjects
		if isfield(sFig.sObject(intSelectedObject).handles,'lines') && sFig.sObject(intObject).drawn == 1
			for p = 1:length(sFig.sObject(intSelectedObject).handles.lines)
				try delete(sFig.sObject(intSelectedObject).handles.lines(p));catch,end
			end
			sFig.sObject(intSelectedObject).handles.lines = [];
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intSelectedObject),'intCenterX') && ~isempty(sDC.ROI(intSelectedObject).intCenterX) && sFig.sObject(intObject).drawn == 1
			try delete(sFig.sObject(intSelectedObject).handles.marker);catch,end
		end
		if isfield(sFig.sObject(intObject).handles,'text') && ~isempty(sFig.sObject(intObject).handles.text)
			try delete(sFig.sObject(intObject).handles.text);catch,end
			sFig.sObject(intObject).handles.text = [];
		end
		sFig.sObject(intSelectedObject).drawn = 0;
	end
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonSelectAllOfPresence_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get structures
	global sDC;
	global sFig;
	
	%check for data
	if isempty(sDC), return;end
	
	%clear selected objects
	DC_unselectAll();
	
	%get data
	intPresence = get(sFig.ptrListPresence, 'Value');
	intObjects = numel(sDC.ROI);
	sFig.vecSelectedObjects = nan(1,intObjects); %over-pre-allocate
	intSelectionObject = 0;
	for intObject=1:intObjects
		if sDC.ROI(intObject).intPresence == intPresence
			intSelectionObject = intSelectionObject + 1;
			sFig.vecSelectedObjects(intSelectionObject) = intObject;
		end
	end
	
	%remove trailing nans
	sFig.vecSelectedObjects = sFig.vecSelectedObjects(1:(find(isnan(sFig.vecSelectedObjects),1,'first')-1));
	
	%lock GUI
	DC_lock(handles);
	
	%set selection & update text
	DC_updateTextInformation;
	
	% remove old drawings and set drawn flags
	for intSelectedObject=sFig.vecSelectedObjects
		if isfield(sFig.sObject(intSelectedObject).handles,'lines') && sFig.sObject(intObject).drawn == 1
			for p = 1:length(sFig.sObject(intSelectedObject).handles.lines)
				try delete(sFig.sObject(intSelectedObject).handles.lines(p));catch,end
			end
			sFig.sObject(intSelectedObject).handles.lines = [];
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intSelectedObject),'intCenterX') && ~isempty(sDC.ROI(intSelectedObject).intCenterX) && sFig.sObject(intObject).drawn == 1
			try delete(sFig.sObject(intSelectedObject).handles.marker);catch,end
		end
		sFig.sObject(intSelectedObject).drawn = 0;
	end
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonSelectAllOfRespType_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get structures
	global sDC;
	global sFig;
	
	%check for data
	if isempty(sDC), return;end
	
	%clear selected objects
	DC_unselectAll();
	
	%get data
	intRespType = get(sFig.ptrListRespType, 'Value');
	intObjects = numel(sDC.ROI);
	sFig.vecSelectedObjects = nan(1,intObjects); %over-pre-allocate
	intSelectionObject = 0;
	for intObject=1:intObjects
		if sDC.ROI(intObject).intRespType == intRespType
			intSelectionObject = intSelectionObject + 1;
			sFig.vecSelectedObjects(intSelectionObject) = intObject;
		end
	end
	
	%remove trailing nans
	sFig.vecSelectedObjects = sFig.vecSelectedObjects(1:(find(isnan(sFig.vecSelectedObjects),1,'first')-1));
	
	%lock GUI
	DC_lock(handles);
	
	%set selection & update text
	DC_updateTextInformation;
	
	% remove old drawings and set drawn flags
	for intSelectedObject=sFig.vecSelectedObjects
		if isfield(sFig.sObject(intSelectedObject).handles,'lines') && sFig.sObject(intObject).drawn == 1
			for p = 1:length(sFig.sObject(intSelectedObject).handles.lines)
				try delete(sFig.sObject(intSelectedObject).handles.lines(p));catch,end
			end
			sFig.sObject(intSelectedObject).handles.lines = [];
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intSelectedObject),'intCenterX') && ~isempty(sDC.ROI(intSelectedObject).intCenterX) && sFig.sObject(intObject).drawn == 1
			try delete(sFig.sObject(intSelectedObject).handles.marker);catch,end
		end
		sFig.sObject(intSelectedObject).drawn = 0;
	end
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonSelectCombo_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get structures
	global sDC;
	global sFig;
	
	%check for data
	if isempty(sDC), return;end
	
	%clear selected objects
	DC_unselectAll();
	
	%get data
	intRespType = get(sFig.ptrListRespType, 'Value');
	intPresence = get(sFig.ptrListPresence, 'Value');
	intSubType = get(sFig.ptrListCellType, 'Value');
	intObjects = numel(sDC.ROI);
	sFig.vecSelectedObjects = nan(1,intObjects); %over-pre-allocate
	intSelectionObject = 0;
	for intObject=1:intObjects
		if sDC.ROI(intObject).intRespType == intRespType && sDC.ROI(intObject).intType == intSubType && sDC.ROI(intObject).intPresence == intPresence
			intSelectionObject = intSelectionObject + 1;
			sFig.vecSelectedObjects(intSelectionObject) = intObject;
		end
	end
	
	%remove trailing nans
	sFig.vecSelectedObjects = sFig.vecSelectedObjects(1:(find(isnan(sFig.vecSelectedObjects),1,'first')-1));
	
	%lock GUI
	DC_lock(handles);
	
	%set selection & update text
	DC_updateTextInformation;
	
	% remove old drawings and set drawn flags
	for intSelectedObject=sFig.vecSelectedObjects
		if isfield(sFig.sObject(intSelectedObject).handles,'lines') && sFig.sObject(intObject).drawn == 1
			for p = 1:length(sFig.sObject(intSelectedObject).handles.lines)
				try delete(sFig.sObject(intSelectedObject).handles.lines(p));catch,end
			end
			sFig.sObject(intSelectedObject).handles.lines = [];
		else
			p = [];
		end
		if isempty(p) && isfield(sDC.ROI(intSelectedObject),'intCenterX') && ~isempty(sDC.ROI(intSelectedObject).intCenterX) && sFig.sObject(intObject).drawn == 1
			try delete(sFig.sObject(intSelectedObject).handles.marker);catch,end
		end
		sFig.sObject(intSelectedObject).drawn = 0;
	end
	
	%redraw
	DC_redraw(0);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrButtonSetSubtypeNr_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%set new subtype number
	
	%lock GUI
	DC_lock(handles);
	
	%move ROIs
	DC_SetSubtypeNr;
	
	%redraw
	DC_redraw(1);
	
	%unlock GUI
	DC_unlock(handles);
end
function ptrPanicButton_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%unlock GUI
	clc;
	DC_unlock(handles);
	
	%set selection & update text
	DC_updateTextInformation;
	
	%show relaxing picture
	try I = imread('DC_KeepCalmAndRelax.tif');catch,return;end
	h=figure;imshow(I)
end
