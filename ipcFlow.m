function varargout = ipcFlow(varargin)
% IPCFLOW MATLAB code for ipcFlow.fig
%      IPCFLOW, by itself, creates a new IPCFLOW or raises the existing
%      singleton*.
%
%      H = IPCFLOW returns the handle to a new IPCFLOW or the handle to
%      the existing singleton*.
%
%      IPCFLOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPCFLOW.M with the given input arguments.
%
%      IPCFLOW('Property','Value',...) creates a new IPCFLOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ipcFlow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ipcFlow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ipcFlow

% Last Modified by GUIDE v2.5 30-Nov-2017 04:08:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ipcFlow_OpeningFcn, ...
                   'gui_OutputFcn',  @ipcFlow_OutputFcn, ...
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


% --- Executes just before ipcFlow is made visible.
function ipcFlow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ipcFlow (see VARARGIN)

% Choose default command line output for ipcFlow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ipcFlow wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ipcFlow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function seglen_Callback(hObject, eventdata, handles)

function seglen_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function overlap_Callback(hObject, eventdata, handles)

function overlap_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in configure.
function configure_Callback(hObject, eventdata, handles)
% hObject    handle to configure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.seglen_=str2double(get(handles.seglen,'String'));              % length of fft window  ---- set to 512 to give a resolution of 2000/256=7.8Hz
handles.overlap_=str2double(get(handles.overlap,'String'));          % length of fft window  ---- set to 512 to give a resolution of 2000/256=7.8Hz
handles.fcut_=str2double(get(handles.overlap,'String'));                   % length of fft window  ---- set to 512 to give a resolution of 2000/256=7.8Hz
handles.fs_=str2double(get(handles.fs,'String')); 
handles.segshift_=(1-handles.overlap_)*handles.seglen_;  % window overlap when doing sliding it across the epoch while doing fft
handles.tres_ = 1000*handles.segshift_/handles.fs_;
handles.fres_ = handles.fs_/handles.seglen_;
set(handles.fres,'String',[num2str(handles.fres_,'%.2f') ' hz']);
set(handles.tres,'String',[num2str(handles.tres_,'%.2f') ' ms']);
handles.timeBin=str2double(get(handles.timeBin,'String'));

function ffrom_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ffrom_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fto_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fto_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fileGet.
function fileGet_Callback(hObject, eventdata, handles)
% hObject    handle to fileGet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.vhdr*','Select the IPC file');

handles.mode_ = 0;
handles.filename_ = FileName;
handles.filepath_ = PathName;
set(handles.mode,'String','Single Subject');
set(handles.cur_file,'String',FileName);
set(handles.cur_dir,'String',PathName);


% --- Executes on button press in folderGet.
function folderGet_Callback(hObject, eventdata, handles)
% hObject    handle to folderGet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PathName = uigetdir('*.mat','Select the IPC file');
set(handles.cur_file,'String','*');
handles.filepath = PathName;
set(handles.mode,'String','Batch Processing');
set(handles.cur_dir,'String',PathName);

function fcut_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fcut_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)

function run_Callback(hObject, eventdata, handles)


% [visit,rest] = strtok(filepath,'\');
% while(~strcmpi(rest,'\')) % Keep running until rest is just the final backslash
%     [visit,rest] = strtok(rest,'\');
% end

fullFile = {};
q = 1;
if (~isempty(strfind(handles.mode.String,'Batch')))
    myFolder = handles.cur_dir.String;
    filePattern = fullfile(myFolder , '*.vhdr');
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
          baseFileName = theFiles(k).name;
          if (isempty(strfind(baseFileName,'rest')))
              fullFile{1,q} = {char(theFiles(k).folder),char(theFiles(k).name),2};
              q = q + 1;
          end
    end
else
    myFolder = handles.cur_dir.String;
    filePattern = fullfile(myFolder , handles.cur_file.String);
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
          baseFileName = theFiles(k).name;
          if (isempty(strfind(baseFileName,'rest')))
              fullFile{1,q} = {char(theFiles(k).folder),char(theFiles(k).name),2};
              q = q + 1;
          end
    end
end
for i=1:length(fullFile)
    set(handles.nsub,'String',['Subject ' num2str(i) '/' num2str(length(fullFile))]);
    cur_ = fullFile{i};
    cur_file = cur_{2};
    cur_folder = [cur_{1} '/'];
    st = 'Status: running';
    set(handles.status,'String',[st char(10) 'Starting EEG...']);
  
    eeglab %opens EEGlab

    set(handles.subject,'String',cur_file);
    % If data is going to be saved as a .set file

%     newName = filename{1};
%     newName(end-4:end) = [];
%     newName = [newName '.set'];

    set(handles.status,'String',[st char(10) 'Loading EEG data...']);
    EEG = pop_loadbv(cur_folder,cur_file);
%     EEG1 = pop_loadbv('.','C10_resting_pre.vhdr');
%     EEG2 = pop_loadbv('.','C10_resting_post.vhdr');
%     OUTEEG = pop_mergeset( EEG1, EEG2);
%     EEG = eeg_regepochs(OUTEEG);
    EEG = eeg_checkset( EEG );

    [~, numEvents] = size([EEG.event.latency]);
    boundaryArray = {'boundary'};
    for n = 2:numEvents
        boundaryArray = [boundaryArray {'boundary'}];
    end
    indicies = strcmpi({EEG.event.type}, boundaryArray);
    latencies = [EEG.event.latency];
    EEG.data = epoch(EEG.data,[latencies(~indicies)], [-5000 10000]);
    % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
    [EEG.event.epoch] = deal(1);
    EEG = eeg_checkset( EEG );
    EEG.data = rmbase(EEG.data, 10000,[1:5000]);

    % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname','TE_C13_RH_TCI_epoch_BC','gui','off');  
    %EEG= pop_loadbv(filepath,filename{1});
    EEG = pop_select( EEG,'nochannel',{'CP34' 'Eog'});
     %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 

    data=EEG.data;
    data1=permute(data,[2 3 1]);                                       % May Not be necessary 
    data1=reshape(data1,[size(data1,1)*size(data1,2)  size(data1,3)]); % either this or the one below

    qq = 1;
    myMap = containers.Map();
    for w=EEG.chanlocs
        myMap(getfield(w,'labels')) = qq;
        qq = qq + 1;
    end

    clear data
    data=data1;
    sampfreq=EEG.srate;       % sampling freq
    % segleng = 768;
    epleng=EEG.pnts;          % length of epoch
    %%%%%% 

    overlp=1-str2double(get(handles.overlap,'String'));               % segment overlap, should be between 0 and 1. ex: .5 = 50% overlap
    % overlp = .1;

    segleng=str2double(get(handles.seglen,'String'));              % length of fft window  ---- set to 512 to give a resolution of 2000/256=7.8Hz
    segshift=floor(overlp*segleng);  % window overlap when doing sliding it across the epoch while doing fft
    segstart=0;               % should theoretically be 0
    maxfreq=str2double(get(handles.fcut,'String'));                  
    t0=0;                     % Stimulus Onset
    art=-1;                   % artifact stuff
    % maxfreq = 50;

    nseg=floor((epleng-segstart-segleng)/segshift) +1;
    [n,nchan]=size(data);
    nepoch=floor(n/epleng);

    set(handles.eplen,'String',epleng);
    set(handles.nep,'String',nepoch);
    art=ones(nepoch,nseg,nchan);
    set(handles.status,'String',[st char(10) 'Calculating cross coherency...']);
    [ndat,nchan]=size(data);
    nepoch=floor(ndat/epleng);
    [dim1,dim2]=size(size(art));

    freqres=sampfreq/segleng;
    highcutmax=floor(segleng/2);   
    if maxfreq>0
        highcut=ceil(maxfreq/(sampfreq/2)*(segleng/2))+1;
    end
    if highcut > highcutmax
        highcut=highcutmax;
    end

    nseg=floor((epleng-segstart-segleng)/segshift) +1; % number of FFT analysis segments
    hanwindows=repmat(hanning(segleng),1,nchan);

    CS=zeros(highcut,nchan,nseg,nchan);

    disp(['nchan=',num2str(nchan),' nepoch=',num2str(nepoch),' nseg=',num2str(nseg),' highcut=',num2str(highcut)]);
    for ep=1:nepoch;
        data_ep=data((ep-1)*epleng+1:ep*epleng,:);
        for seg=1:nseg
          artseg=reshape(art(ep,seg,:),1,nchan);
          artseg=repmat(artseg,segleng,1);
          tmp=data_ep((seg-1)*segshift+segstart+1:segleng+(seg-1)*segshift+segstart,:);
          tmp=detrend(tmp).*hanwindows.*artseg; 
          EEG2=fft(tmp);  
          sp=EEG2(1:highcut,:); 
          CSadd=reshape(repmat(sp,1,nchan),highcut,nchan,nchan).*conj(reshape(repmat(sp,nchan,1),highcut,nchan,nchan));
          CS(:,:,seg,:)=CS(:,:,seg,:)+reshape(CSadd,highcut,nchan,1,nchan);
        end;     
    end;
    set(handles.status,'String',[st char(10) 'Calculating IPC...']);
    chanpar=-18;
    [CSbase,COy,CO]=cs_ana(CS,sampfreq,maxfreq,segleng,segshift,chanpar);
    COout=cs2csbase(COy,1);
    set(handles.status,'String',[st char(10) 'done.']);
    data=COout;
    aa = strsplit(cur_file,'.');
    bb = char([cur_folder aa{1} '-orig.mat']);
    save(bb,'data')

    aa = strsplit(cur_file,'.');
    bb = char([cur_folder aa{1} '-map.mat']);
    save(bb,'myMap')

    [b,nc1,a,nc2]=size(data);
    ffrom = str2double(get(handles.ffrom,'String'));
    fto = str2double(get(handles.fto,'String'));
    i1 = -1;
    i2 = 1000000001;
    for i=1:a
        for j=1:b
            timePoints(i,j)=(i-1)*segshift+segleng/2;
            tmp = (j-1)*maxfreq/(b-1);
            freqPoints(i,j)= tmp;
            if (i1<0 && tmp>=ffrom)
                i1=j;
            end
            if (i2>1000000000 && tmp>fto)
                i2=j;
            end

        end
    end
    if (i1<0)
        i1=1;
    end
    if (i2>1000000000)
        i2=b;
    end
    t = timePoints(:,1)/sampfreq;
    Ts = t(2)-t(1);
    Fs = 1/Ts;
    val = squeeze(mean(data(i1:i2,:,:,:),1));
    data = {};
    data.t = t;
    data.val = val;
    aa = strsplit(cur_file,'.');
    bb = char([cur_folder aa{1} '-post.mat']);
    save(bb,'data')

    % Save CSV -timestats
    S = 100/1000;
    npoints = floor(S/Ts);
    nbins = floor(length(t)/npoints);
    fid = fopen(char([cur_folder aa{1} '-timebins.csv']),'wt');
    fprintf(fid, 'timebin,c1,c2,mean,median,peak\n');

    for w=myMap.keys()
        for v=myMap.keys()
            i = myMap(w{1});
            j = myMap(v{1});
            before = 1;
            next = npoints;
            for p=1:nbins
                MEAN = mean(abs(imag(val(i,before:next,j))));
                MED = median(abs(imag(val(i,before:next,j))));
                PEAK = mean(abs(imag(val(i,before:next,j))));
                fprintf(fid, '%d,%s,%s,%d,%d,%d\n', p,w{1},v{1},MEAN,MED,PEAK);
                before = next;
                next = next + npoints;
            end
        end
    end
    fclose(fid);
    k=0;
    set(handles.status,'String',[st char(10) 'Calculating periodicity...']);

    % whole data
    fid = fopen(char([cur_folder aa{1} '-cyclic.csv']),'wt');
    fprintf(fid, 'mode,c1,c2,peak,freq\n');
    for w=myMap.keys()
        for v=myMap.keys()
            i = myMap(w{1});
            j = myMap(v{1});
            cur =  squeeze(abs(imag(val(i,:,j))));
            cur = cur - mean(cur);
            N = length(cur);
            xdft = fft(cur);
            xdft = xdft(1:N/2+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            if k==0
                k=1;
                FF = zeros(nc1,length(psdx),nc2);
                FF(i,:,j) = psdx;
                freq = 0:Fs/N:Fs/2;
            else
            FF(i,:,j) = psdx;
            end
            [max_val, idx] = max(psdx);
            maxf = freq(idx);
            fprintf(fid, 'all,%s,%s,%d,%d\n',w{1},v{1},max_val,maxf);
        end
    end
    data = {};
    data.freq = freq;
    data.val = FF;
    aa = strsplit(cur_file,'.');
    bb = char([cur_folder aa{1} '-cycle.mat']);
    save(bb,'data')
    


    % prestim
    k = 0;
    lcut = floor(size(val,2)/3);

    for w=myMap.keys()
        for v=myMap.keys()
            i = myMap(w{1});
            j = myMap(v{1});
            cur =  squeeze(abs(imag(val(i,1:lcut,j))));
            cur = cur - mean(cur);
            N = length(cur);
            xdft = fft(cur);
            xdft = xdft(1:N/2+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            if k==0
                k=1;
                FF = zeros(nc1,length(psdx),nc2);
                FF(i,:,j) = psdx;
                freq = 0:Fs/N:Fs/2;
            else
            FF(i,:,j) = psdx;
            end
            [max_val, idx] = max(psdx);
            maxf = freq(idx);
            fprintf(fid, 'pre-stim,%s,%s,%d,%d\n',w{1},v{1},max_val,maxf);
        end
    end
    data = {};
    data.freq = freq;
    data.val = FF;



    % Post Stim
    k = 0;
    for w=myMap.keys()
        for v=myMap.keys()
            i = myMap(w{1});
            j = myMap(v{1});
            cur =  squeeze(abs(imag(val(i,lcut:end,j))));
            cur = cur - mean(cur);
            N = length(cur);
            xdft = fft(cur);
            xdft = xdft(1:N/2+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            if k==0
                k=1;
                FF = zeros(nc1,length(psdx),nc2);
                FF(i,:,j) = psdx;
                freq = 0:Fs/N:Fs/2;
            else
            FF(i,:,j) = psdx;
            end
            [max_val, idx] = max(psdx);
            maxf = freq(idx);
            fprintf(fid, 'post-stim,%s,%s,%d,%d\n',w{1},v{1},max_val,maxf);
        end
    end
    data = {};
    data.freq = freq;
    data.val = FF;
    fclose(fid);
end
set(handles.status,'String',[st char(10) 'Done.']);
        




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


% --- Executes on button press in synthOn.
function synthOn_Callback(hObject, eventdata, handles)
% hObject    handle to synthOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of synthOn



function c1f_Callback(hObject, eventdata, handles)
% hObject    handle to c1f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1f as text
%        str2double(get(hObject,'String')) returns contents of c1f as a double


% --- Executes during object creation, after setting all properties.
function c1f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c2f_Callback(hObject, eventdata, handles)
% hObject    handle to c2f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c2f as text
%        str2double(get(hObject,'String')) returns contents of c2f as a double


% --- Executes during object creation, after setting all properties.
function c2f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c2f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in modulate.
function modulate_Callback(hObject, eventdata, handles)
% hObject    handle to modulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modulate



function synthID_Callback(hObject, eventdata, handles)
% hObject    handle to synthID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of synthID as text
%        str2double(get(hObject,'String')) returns contents of synthID as a double


% --- Executes during object creation, after setting all properties.
function synthID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to synthID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function da_Callback(hObject, eventdata, handles)
% hObject    handle to da (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of da as text
%        str2double(get(hObject,'String')) returns contents of da as a double


% --- Executes during object creation, after setting all properties.
function da_CreateFcn(hObject, eventdata, handles)
% hObject    handle to da (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rest.
function rest_Callback(hObject, eventdata, handles)
% hObject    handle to rest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rest



function timeBin_Callback(hObject, eventdata, handles)
% hObject    handle to timeBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeBin as text
%        str2double(get(hObject,'String')) returns contents of timeBin as a double


% --- Executes during object creation, after setting all properties.
function timeBin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
