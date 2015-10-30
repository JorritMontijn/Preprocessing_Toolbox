function varargout = prepro_phaseCorr(varargin)
% PREPRO_PHASECORR MATLAB code for prepro_phaseCorr.fig
%      PREPRO_PHASECORR, by itself, creates a new PREPRO_PHASECORR or raises the existing
%      singleton*.
%
%      H = PREPRO_PHASECORR returns the handle to a new PREPRO_PHASECORR or the handle to
%      the existing singleton*.
%
%      PREPRO_PHASECORR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPRO_PHASECORR.M with the given input arguments.
%
%      PREPRO_PHASECORR('Property','Value',...) creates a new PREPRO_PHASECORR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prepro_phaseCorr_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prepro_phaseCorr_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prepro_phaseCorr

% Last Modified by GUIDE v2.5 11-Dec-2012 17:08:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prepro_phaseCorr_OpeningFcn, ...
                   'gui_OutputFcn',  @prepro_phaseCorr_OutputFcn, ...
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

setappdata(gcf,'imRaw',varargin{1});

% --- Executes just before prepro_phaseCorr is made visible.
function prepro_phaseCorr_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prepro_phaseCorr (see VARARGIN)

% Choose default command line output for prepro_phaseCorr
handles.output = hObject;

imRaw = getappdata(gcf,'imRaw');
h=figure;
%axes(handles.axes1);
imshow(imRaw);
setappdata(gcf,'dblCorr',0);
setappdata(gcf,'sizeX',size(imRaw,2));

minCorr = get(hObject,'Min');
maxCorr = get(hObject,'Max');
set(handles.slider1,'Value',((maxCorr-minCorr)/2 + minCorr));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prepro_phaseCorr wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prepro_phaseCorr_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = getappdata(gcf,'dblCorr');


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

minCorr = get(hObject,'Min');
maxCorr = get(hObject,'Max');
curCorr = get(hObject,'Value');

dblFrac = (curCorr - minCorr) / (maxCorr - minCorr) - 0.5;
sizeX = getappdata(gcf,'sizeX');

dblCorr = round((sizeX-1)*dblFrac) + 1;

imRaw = getappdata(gcf,'imRaw');
%axes(handles.axes1);

imCorr = doPhaseCorrect(imRaw,dblCorr);
imshow(imCorr);
setappdata(gcf,'dblCorr',dblCorr);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
