function varargout = ipcPlot(varargin)
% IPCPLOT MATLAB code for ipcPlot.fig
%      IPCPLOT, by itself, creates a new IPCPLOT or raises the existing
%      singleton*.
%
%      H = IPCPLOT returns the handle to a new IPCPLOT or the handle to
%      the existing singleton*.
%
%      IPCPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPCPLOT.M with the given input arguments.
%
%      IPCPLOT('Property','Value',...) creates a new IPCPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ipcPlot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ipcPlot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ipcPlot

% Last Modified by GUIDE v2.5 10-Oct-2017 03:06:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ipcPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @ipcPlot_OutputFcn, ...
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


% --- Executes just before ipcPlot is made visible.
function ipcPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ipcPlot (see VARARGIN)

% Choose default command line output for ipcPlot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ipcPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ipcPlot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fs_Callback(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs as text
%        str2double(get(hObject,'String')) returns contents of fs as a double


% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seglen_Callback(hObject, eventdata, handles)
% hObject    handle to seglen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seglen as text
%        str2double(get(hObject,'String')) returns contents of seglen as a double


% --- Executes during object creation, after setting all properties.
function seglen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seglen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function overlap_Callback(hObject, eventdata, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of overlap as text
%        str2double(get(hObject,'String')) returns contents of overlap as a double


% --- Executes during object creation, after setting all properties.
function overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PathName = get(handles.cur_dir,'String');

base = get(handles.base,'String');
sampfreq=str2double(get(handles.fs,'String'));
segleng=str2double(get(handles.seglen,'String'));
overlap=str2double(get(handles.overlap,'String'));
segshift=floor((1-overlap)*segleng);
segleng = floor(segleng);
t0=0;
chan1=get(handles.C1_,'String');
chan2=get(handles.C2,'String');
mapfile  = [PathName base '-map.mat'];
mapper1 = load(mapfile);
mapper = mapper1.myMap;
chan1=mapper(chan1);
chan2=mapper(chan2);
tmp = [PathName base '-orig.mat'];
tmp = load(tmp);
CS = abs(imag(tmp.data));
axes(handles.ax1);
[timePointsX,freqPointsY,finalData] = plot_coh_onepair2(CS,sampfreq,segleng,segshift,t0,chan1,chan2);
pan on
datacursormode on

axes(handles.ax2);
tmp = [PathName base '-post.mat'];
tmp = load(tmp);
res = abs(imag(tmp.data.val));
plot(timePointsX,squeeze(res(chan1,:,chan2)));
xlabel('Time [ms]');
ylabel('Average IPC in the band');

axes(handles.ax3);
tmp = [PathName base '-cycle.mat'];
tmp = load(tmp);
res = tmp.data.val;
plot(tmp.data.freq,squeeze(res(chan1,:,chan2)));
xlabel('Freq. [Hz]');
ylabel('Periodicity of absolute IPC');




function C1__Callback(hObject, eventdata, handles)
% hObject    handle to C1_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C1_ as text
%        str2double(get(hObject,'String')) returns contents of C1_ as a double


% --- Executes during object creation, after setting all properties.
function C1__CreateFcn(hObject, eventdata, handles)
% hObject    handle to C1_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C2_Callback(hObject, eventdata, handles)
% hObject    handle to C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C2 as text
%        str2double(get(hObject,'String')) returns contents of C2 as a double


% --- Executes during object creation, after setting all properties.
function C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.vhdr*','Select the IPC file');
base = strsplit(FileName,'.');
base = base{1};
set(handles.base,'String',base);
set(handles.cur_dir,'String',PathName);
