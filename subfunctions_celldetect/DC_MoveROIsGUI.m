function varargout = DC_MoveROIsGUI(varargin)
	% DC_MOVEROISGUI MATLAB code for DC_MoveROIsGUI.fig
	%      DC_MOVEROISGUI, by itself, creates a new DC_MOVEROISGUI or raises the existing
	%      singleton*.
	%
	%      H = DC_MOVEROISGUI returns the handle to a new DC_MOVEROISGUI or the handle to
	%      the existing singleton*.
	%
	%      DC_MOVEROISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
	%      function named CALLBACK in DC_MOVEROISGUI.M with the given input arguments.
	%
	%      DC_MOVEROISGUI('Property','Value',...) creates a new DC_MOVEROISGUI or raises the
	%      existing singleton*.  Starting from the left, property value pairs are
	%      applied to the GUI before DC_MoveROIsGUI_OpeningFcn gets called.  An
	%      unrecognized property name or invalid value makes property application
	%      stop.  All inputs are passed to DC_MoveROIsGUI_OpeningFcn via varargin.
	%
	%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
	%      instance to run (singleton)".
	%
	% See also: GUIDE, GUIDATA, GUIHANDLES
	
	% Edit the above text to modify the response to help DC_MoveROIsGUI
	
	% Last Modified by GUIDE v2.5 22-Apr-2014 12:34:25
	
	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
		'gui_Singleton',  gui_Singleton, ...
		'gui_OpeningFcn', @DC_MoveROIsGUI_OpeningFcn, ...
		'gui_OutputFcn',  @DC_MoveROIsGUI_OutputFcn, ...
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
% --- Executes just before DC_MoveROIsGUI is made visible.
function DC_MoveROIsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
	% This function has no output args, see OutputFcn.
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	% varargin   command line arguments to DC_MoveROIsGUI (see VARARGIN)
	
	% Choose default command line output for DC_MoveROIsGUI
	handles.output = hObject;
	
	% Update handles structure
	guidata(hObject, handles);
	
	% UIWAIT makes DC_MoveROIsGUI wait for user response (see UIRESUME)
	% uiwait(handles.figure1);
	
end
% --- Outputs from this function are returned to the command line.
function varargout = DC_MoveROIsGUI_OutputFcn(hObject, eventdata, handles)
	% varargout  cell array for returning output args (see VARARGOUT);
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	
	% Get default command line output from handles structure
	varargout{1} = handles.output;
	
end
function ptrButtonDone_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%get globals
	global sMoveROI
	
	%set for exit
	sMoveROI.boolRunning = false;
	
	%done; close GUI
	close;
end
function ptrButtonLeft_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	
	%lock GUI
	DC_lock(handles);
	
	%move left
	intIncrement = getIncrement(handles);
	global sMoveROI;
	sMoveROI.boolUpdated = true;
	sMoveROI.intROIDisplacementX = sMoveROI.intROIDisplacementX - intIncrement;
	
	%lock GUI
	DC_unlock(handles);
end
function ptrButtonUp_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%lock GUI
	DC_lock(handles);
	
	%move up
	intIncrement = getIncrement(handles);
	global sMoveROI;
	sMoveROI.boolUpdated = true;
	sMoveROI.intROIDisplacementY = sMoveROI.intROIDisplacementY + intIncrement;
	
	%lock GUI
	DC_unlock(handles);
end
function ptrButtonDown_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%lock GUI
	DC_lock(handles);
	
	%move down
	intIncrement = getIncrement(handles);
	global sMoveROI;
	sMoveROI.boolUpdated = true;
	sMoveROI.intROIDisplacementY = sMoveROI.intROIDisplacementY - intIncrement;
	
	%lock GUI
	DC_unlock(handles);
end
function ptrButtonRight_Callback(hObject, eventdata, handles) %#ok<DEFNU>
	%lock GUI
	DC_lock(handles);
	
	%move right
	intIncrement = getIncrement(handles);
	global sMoveROI;
	sMoveROI.boolUpdated = true;
	sMoveROI.intROIDisplacementX = sMoveROI.intROIDisplacementX + intIncrement;
	
	%lock GUI
	DC_unlock(handles);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok<DEFNU>
	% hObject    handle to figure1 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	%get globals
	global sMoveROI
	
	%set for exit
	sMoveROI.boolRunning = false;
	
	% Hint: delete(hObject) closes the figure
	delete(hObject);
end



function ptrEditPixelsPerShift_Callback(hObject, eventdata, handles)
% hObject    handle to ptrEditPixelsPerShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	intIncrement = round(str2double(get(handles.ptrEditPixelsPerShift,'String')));
	if isempty(intIncrement) || isnan(intIncrement),intIncrement = 1;end
	set(handles.ptrEditPixelsPerShift,'String',num2str(intIncrement));
end
% --- Executes during object creation, after setting all properties.
function ptrEditPixelsPerShift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptrEditPixelsPerShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function intIncrement = getIncrement(handles)
	intIncrement = round(str2double(get(handles.ptrEditPixelsPerShift,'String')));
	if isempty(intIncrement),intIncrement = 1;end
end
