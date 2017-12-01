function varargout = IPC_kamal(varargin)
% IPC_KAMAL MATLAB code for IPC_kamal.fig
%      IPC_KAMAL, by itself, creates a new IPC_KAMAL or raises the existing
%      singleton*.
%
%      H = IPC_KAMAL returns the handle to a new IPC_KAMAL or the handle to
%      the existing singleton*.
%
%      IPC_KAMAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPC_KAMAL.M with the given input arguments.
%
%      IPC_KAMAL('Property','Value',...) creates a new IPC_KAMAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IPC_kamal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IPC_kamal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IPC_kamal

% Last Modified by GUIDE v2.5 08-Oct-2017 19:42:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IPC_kamal_OpeningFcn, ...
                   'gui_OutputFcn',  @IPC_kamal_OutputFcn, ...
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


% --- Executes just before IPC_kamal is made visible.
function IPC_kamal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IPC_kamal (see VARARGIN)

% Choose default command line output for IPC_kamal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IPC_kamal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IPC_kamal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn.
function btn_Callback(hObject, eventdata, handles)
% hObject    handle to btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
postIPC(handles);


function fmin_Callback(hObject, eventdata, handles)
% hObject    handle to fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fmin as text
%        str2double(get(hObject,'String')) returns contents of fmin as a double


% --- Executes during object creation, after setting all properties.
function fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fmax_Callback(hObject, eventdata, handles)
% hObject    handle to fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fmax as text
%        str2double(get(hObject,'String')) returns contents of fmax as a double


% --- Executes during object creation, after setting all properties.
function fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function tmin_Callback(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmin as text
%        str2double(get(hObject,'String')) returns contents of tmin as a double


% --- Executes during object creation, after setting all properties.
function tmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C1_Callback(hObject, eventdata, handles)
% hObject    handle to C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C1 as text
%        str2double(get(hObject,'String')) returns contents of C1 as a double


% --- Executes during object creation, after setting all properties.
function C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tmax_Callback(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmax as text
%        str2double(get(hObject,'String')) returns contents of tmax as a double


% --- Executes during object creation, after setting all properties.
function tmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmax (see GCBO)
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



function dataset_Callback(hObject, eventdata, handles)
% hObject    handle to dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataset as text
%        str2double(get(hObject,'String')) returns contents of dataset as a double
[FileName,PathName] = uigetfile('*.mat','Select the IPC file');
% --- Executes during object creation, after setting all properties.
function dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function postIPC(h)
overlp=0.1;               % segment overlap, should be between 0 and 1. ex: .5 = 50% overlap

sampfreq=str2double(get(h.fs,'String'));            % sampling freq
segleng=768;              % length of fft window  ---- set to 512 to give a resolution of 2000/256=7.8Hz
segshift=overlp*segleng;  % window overlap when doing sliding it across the epoch while doing fft
segstart=0;               % should theoretically be 0
maxfreq=50;                  
t0=str2double(get(h.tmin,'String')); % Stimulus Onset
tmax=str2double(get(h.tmax,'String')); 
art=-1;                   % artifact stuff
freq = 3; 
chan1=str2num(['uint8(' get(h.C1,'String') ')']);
chan2=str2num(['uint8(' get(h.C2,'String') ')']);
set(h.results,'String',class(chan1));
nseg=13;
fmin = str2double(get(h.fmin,'String'));
fmax = str2double(get(h.fmax,'String'));
Fpv = get(h.Fpv,'String');
Fnv = get(h.Fnv,'String');
cur = load([Fpv Fnv]);
data = imag(cur.data);
[timePoints,freqPoints,finalData] = coh_onepair(data,sampfreq,segleng,segshift,t0,chan1,chan2);
[a,b]= size(freqPoints);
f_i1 = -1;
f_i2 = -1;
t_i1 = -1;
t_i2 = -1;
for i=1:b
    if freqPoints(1,i)>=fmin && f_i1<0
        f_i1=i;
    end
    if freqPoints(1,i)>fmax && f_i2<0
        f_i2=i-1;
    end
end
for i=1:a
    if timePoints(i,1)>=t0 && t_i1<0
        t_i1=i;
    end
    if timePoints(i,1)>tmax && t_i2<0
        t_i2=i-1;
    end
end
if timePoints(i,1)<=tmax
    t_i2=i;
end
if t_i1<0 || t_i2<0
    set(h.results,'String','Not enough data in time domain')
    return
end
if f_i1<0 || f_i2<0
    set(h.results,'String','Not enough data in frequency domain')
    return
end
tmp =finalData(t_i1:t_i2,f_i1:f_i2);
[nt,nf] = size(tmp);
res1 = ['#samples(time): ' num2str(nt) char(10) '#samples(freq): ' num2str(nf)];
set(h.results,'String',res1)
time_series = mean(tmp,2);
axes(h.ax1);
plot(timePoints(t_i1:t_i2),time_series);
xlabel('Time(ms)')
ylabel('IPC')
axes(h.ax2);
Fs = 1000*size(data,3)/max(max(timePoints));
N = length(time_series);
xdft = fft(time_series);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

plot(freq,psdx)
grid on
xlabel('Frequency(Hz)')
ylabel('Power')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.mat','Select the IPC file');
if (FileName ~= 0)
set(handles.Fpv,'String',PathName);
set(handles.Fnv,'String',FileName);
set(handles.Fp,'String','IPC dataset file path:');
set(handles.Fn,'String','IPC dataset file:');
end

function [timePoints,freqPoints,finalData] = coh_onepair(CO,sampfreq,segleng,segshift,t0,chan1,chan2);

[highcut,nchan,nsegs,nchan]=size(CO);

T=segleng/sampfreq;
maxfreqexact=(highcut-1)/T;
t0inms=t0*1000/sampfreq;
CO_here=(reshape(CO(:,chan1,:,chan2),highcut,nsegs))';
[timePoints,freqPoints,finalData] = nprl(CO_here,maxfreqexact,t0inms,segleng*1000/sampfreq,segshift*1000/sampfreq,'');

function [timePoints,freqPoints,finalData] = nprl(data,maxfreq,t0inms,segleng,segshift,titlename)
finalData = data;
[a,b]=size(data);
clear x y;
for i=1:a
    for j=1:b
        timePoints(i,j)=(i-1)*segshift+segleng/2-t0inms;
        freqPoints(i,j)=(j-1)*maxfreq/(b-1);
    end
end
