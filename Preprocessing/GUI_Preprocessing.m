function varargout = GUI_Preprocessing(varargin)
% GUI_PREPROCESSING MATLAB code for GUI_Preprocessing.fig
%      GUI_PREPROCESSING, by itself, creates a new GUI_PREPROCESSING or raises the existing
%      singleton*.
%
%      H = GUI_PREPROCESSING returns the handle to a new GUI_PREPROCESSING or the handle to
%      the existing singleton*.
%
%      GUI_PREPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PREPROCESSING.M with the given input arguments.
%
%      GUI_PREPROCESSING('Property','Value',...) creates a new GUI_PREPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Preprocessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Preprocessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Preprocessing

% Last Modified by GUIDE v2.5 24-Apr-2020 10:00:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Preprocessing_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Preprocessing_OutputFcn, ...
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


% --- Executes just before GUI_Preprocessing is made visible.
function GUI_Preprocessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Preprocessing (see VARARGIN)

% Choose default command line output for GUI_Preprocessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Preprocessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = GUI_Preprocessing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function rat_number_Callback(hObject, eventdata, handles)
% hObject    handle to rat_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_number as text
%        str2double(get(hObject,'String')) returns contents of rat_number as a double

handles = guidata(hObject);
handles.GenParams.Rat = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rat_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.GenParams.Rat = str2double(get(hObject,'String'));
guidata(hObject, handles);

function headstage_Callback(hObject, eventdata, handles)
% hObject    handle to headstage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of headstage as text
%        str2double(get(hObject,'String')) returns contents of headstage as a double
handles = guidata(hObject);
handles.GenParams.HeadStage = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function headstage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to headstage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.GenParams.HeadStage = str2double(get(hObject,'String'));
guidata(hObject, handles);

function session_name_Callback(hObject, eventdata, handles)
% hObject    handle to session_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_name as text
%        str2double(get(hObject,'String')) returns contents of session_name as a double
handles = guidata(hObject);
handles.GenParams.Session = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function session_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.GenParams.Session = get(hObject,'String');
guidata(hObject, handles);


function subsession_name_Callback(hObject, eventdata, handles)
% hObject    handle to subsession_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subsession_name as text
%        str2double(get(hObject,'String')) returns contents of subsession_name as a double
handles = guidata(hObject);
handles.GenParams.SubSession = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function subsession_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subsession_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.GenParams.SubSession = get(hObject,'String');
guidata(hObject, handles)



% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.availability,'String','Busy');

% Retrieved data from Excel file on the viability of the channels 
[handles.ChParams] = retriveElect3Info(handles.GenParams.Rat,handles.GenParams.Session);

disp('Data from channels information are retrieved')

% Where are the data
[file,path] = uigetfile('*.mat',...
    'Select All files that you want to concanenate', ...
    'MultiSelect', 'on');
if isequal(file,0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(path,file)]);
end



% with implemented wait bar
[handles.LFP,handles.ElectParams] = Concatenate_data(path,file,handles.GenParams.HeadStage);

handles.ElectParams.SR = handles.ElectParams.ch_KHz*1000;

disp('Loading done')

handles.LFP=double(handles.LFP)*handles.ElectParams.ch_bitResolution/handles.ElectParams.ch_gain/1000; % Convert bit into mV!

disp('conversion is done')

% set needed variable 
handles.struct.EMG = handles.LFP(4,:);

handles.scale = 2;
handles.ElectParams.filtering.HP = 0.1;
handles.ElectParams.filtering.LP = 500;

% display the path of data saving
handles.pathSav = 'D:\Documents\MatlabData_Elect4';
handles.nameFolder = [num2str(handles.GenParams.Rat),char(handles.GenParams.Session),char(handles.GenParams.SubSession)];
mkdir(handles.pathSav,handles.nameFolder) % folder creation


% create file for storing path and rats
if exist(fullfile(handles.pathSav,'managingRat.mat'),'file')
    %counting the number of folder
    all_files = dir(handles.pathSav);
    all_dir = all_files([all_files(:).isdir]);
    num_dir = numel(all_dir);

    load(fullfile(handles.pathSav,'managingRat.mat'))
    
    % if already exist in the list don't do anything
    for i = 1: length(managingRat(4,:))
        lockfile(i,1) = mean(categorical(cellstr(managingRat{4,i})) == categorical(cellstr(fullfile(handles.pathSav,handles.nameFolder))));
    end
    lenght_managingRat = length(managingRat(1,:))+1;
    if isempty(find(lockfile == 1,1)) == 1
            managingRat{1,lenght_managingRat} = str2double(get(handles.rat_number,'String'));
            managingRat{2,lenght_managingRat} = get(handles.session_name,'String');
            managingRat{3,lenght_managingRat} = get(handles.subsession_name,'String');
            managingRat{4,lenght_managingRat} = fullfile(handles.pathSav,handles.nameFolder);
            namePath= ['managingRat.mat'];
            save(fullfile(handles.pathSav,namePath),'managingRat')
     end
else
    managingRat{1,1} = str2double(get(handles.rat_number,'String'));
    managingRat{2,1} = get(handles.session_name,'String');
    managingRat{3,1} = get(handles.subsession_name,'String');
    managingRat{4,1} = fullfile(handles.pathSav,handles.nameFolder);
    namePath= ['managingRat.mat'];
    save(fullfile(handles.pathSav,namePath),'managingRat')  
end


handles.pathSav= fullfile(handles.pathSav,handles.nameFolder);

if exist(fullfile(handles.pathSav,'temp_local_LFP.mat'),'file')
    disp('The file temp_local_LFP.mat already exist')
    load(fullfile(handles.pathSav,'temp_local_LFP.mat'))
    handles.struct.list = LFP_local;
else
    LFP_local = {'LHb','dHPC','BLA','ACC','PRL'};
    handles.struct.list = LFP_local;
    save(fullfile(handles.pathSav,'temp_local_LFP.mat'),'LFP_local');
end



set(handles.edit_save_path,'String',handles.pathSav);

set(handles.availability,'String','Free');
guidata(hObject, handles);

function scale_Callback(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale as text
%        str2double(get(hObject,'String')) returns contents of scale as a double
handles = guidata(hObject);
handles.scale = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.scale = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes on selection change in select_struct.
function select_struct_Callback(hObject, eventdata, handles)
% hObject    handle to select_struct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_struct contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_struct
handles = guidata(hObject);

switch get(handles.select_struct,'Value')
    case 1
        handles.struct.current = {'LHb'};
        handles.current.lock = find(handles.ChParams.Structures == 'LHb');
        
        set(handles.ch_axe1,'String',handles.ChParams.RatData(handles.current.lock(1)));
        set(handles.ch_axe2,'String',handles.ChParams.RatData(handles.current.lock(2)));
        set(handles.ch_axe3,'String',handles.ChParams.RatData(handles.current.lock(3)));
        disp('you choose LHb')
        
    case 2
        handles.struct.current = {'dHPC'};
        handles.current.lock = find(handles.ChParams.Structures == 'dHPC');
        
        set(handles.ch_axe1,'String',handles.ChParams.RatData(handles.current.lock(1)));
        set(handles.ch_axe2,'String',handles.ChParams.RatData(handles.current.lock(2)));
        set(handles.ch_axe3,'String',handles.ChParams.RatData(handles.current.lock(3)));
        disp('you choose dHPC')


    case 3
        handles.struct.current = {'BLA'};
        handles.current.lock = find(handles.ChParams.Structures == 'BLA');
        
        set(handles.ch_axe1,'String',handles.ChParams.RatData(handles.current.lock(1)));
        set(handles.ch_axe2,'String',handles.ChParams.RatData(handles.current.lock(2)));
        set(handles.ch_axe3,'String',handles.ChParams.RatData(handles.current.lock(3)));
        disp('you choose BLA')


    case 4
        handles.struct.current = {'ACC'};
        handles.current.lock = find(handles.ChParams.Structures == 'ACC');
        
        set(handles.ch_axe1,'String',handles.ChParams.RatData(handles.current.lock(1)));
        set(handles.ch_axe2,'String',handles.ChParams.RatData(handles.current.lock(2)));
        set(handles.ch_axe3,'String',handles.ChParams.RatData(handles.current.lock(3)));
        disp('you choose ACC')


    case 5
        handles.struct.current = {'PRL'};
        handles.current.lock = find(handles.ChParams.Structures == 'PRL');

        set(handles.ch_axe1,'String',handles.ChParams.RatData(handles.current.lock(1)));
        set(handles.ch_axe2,'String',handles.ChParams.RatData(handles.current.lock(2)));
        set(handles.ch_axe3,'String',handles.ChParams.RatData(handles.current.lock(3)));
        disp('you choose PRL')

        
        
    otherwise
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function select_struct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_struct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in select_source.
function select_source_Callback(hObject, eventdata, handles)
% hObject    handle to select_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_source contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_source
handles = guidata(hObject);


switch get(handles.select_source,'Value')
    case 1
        handles.current.source = 'raw LFP';
        handles.current.LFP = handles.LFP(handles.current.lock,:);
        % time calculation
        dt = 1/handles.ElectParams.SR;
        tottime = length(handles.current.LFP(1,:))/handles.ElectParams.SR;
        handles.current.time = 0:dt:tottime-dt;
        handles.current.SR = handles.ElectParams.SR;
        disp(' you choose raw data')
        
    case 2
        handles.current.source = 'filtered LFP';
        handles.current.LFP = handles.LFP(handles.current.lock,:);
        % time calculation
        dt = 1/handles.ElectParams.SR;
        tottime = length(handles.current.LFP(1,:))/handles.ElectParams.SR;
        handles.current.time = 0:dt:tottime-dt;
        handles.current.SR = handles.ElectParams.SR;
        disp(' you choose filtered data')

    case 3
        handles.current.source = 'downsampled LFP';
        handles.current.LFP = handles.LFP(handles.current.lock,:);
        % time calculation
        dt = 1/handles.ElectParams.New_SR;
        tottime = length(handles.current.LFP(1,:))/handles.ElectParams.New_SR;
        handles.current.time = 0:dt:tottime-dt;
        handles.current.SR = handles.ElectParams.New_SR;
        disp(' you choose downsampled data')
  
        
    otherwise
end
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function select_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_track.
function plot_track_Callback(hObject, eventdata, handles)
% hObject    handle to plot_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

cla(handles.axes1)
axes(handles.axes1)
plot(handles.current.time,handles.current.LFP(1,:),'Color',rgb('DarkBlue')); hold on;
title(['Channel n°1 ', cell2mat(handles.struct.current) ,' ',handles.current.source]); ylim([-handles.scale handles.scale]); hold off

cla(handles.axes2)
axes(handles.axes2)
plot(handles.current.time,handles.current.LFP(2,:),'Color',rgb('FireBrick')); hold on;
title(['Channel n°2 ', cell2mat(handles.struct.current) ,' ',handles.current.source]); ylim([-handles.scale handles.scale]); hold off

cla(handles.axes3)
axes(handles.axes3)
plot(handles.current.time,handles.current.LFP(3,:),'Color',rgb('DarkViolet')); hold on;
title(['Channel n°3 ',cell2mat(handles.struct.current) ,' ',handles.current.source]); ylim([-handles.scale handles.scale]); hold off

guidata(hObject, handles);



function hp_filtering_Callback(hObject, eventdata, handles)
% hObject    handle to hp_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hp_filtering as text
%        str2double(get(hObject,'String')) returns contents of hp_filtering as a double
handles = guidata(hObject);
handles.ElectParams.filtering.HP = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function hp_filtering_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hp_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.ElectParams.filtering.HP = str2double(get(hObject,'String'));
guidata(hObject, handles);



function lp_filtering_Callback(hObject, eventdata, handles)
% hObject    handle to lp_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lp_filtering as text
%        str2double(get(hObject,'String')) returns contents of lp_filtering as a double
handles = guidata(hObject);
handles.ElectParams.filtering.LP = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lp_filtering_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lp_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.ElectParams.filtering.LP = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes on button press in apply_filtering.
function apply_filtering_Callback(hObject, eventdata, handles)
% hObject    handle to apply_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
set(handles.availability,'String','Busy');
f = waitbar(0,'Please wait...');

[handles.ElectParams.Nchans, handles.ElectParams.Nsamples] = size(handles.LFP);
handles.LFP_filtered=zeros(handles.ElectParams.Nchans,handles.ElectParams.Nsamples);


% Filtering the Data - Nelson Filter
for k=1:handles.ElectParams.Nchans
    waitbar(k/handles.ElectParams.Nchans,f,['Processins channel number ', num2str(k), ' on ', num2str(handles.ElectParams.Nchans)]);
    handles.LFP_filtered(k,:) = ck_filt(handles.LFP(k,:),1375,[handles.ElectParams.filtering.HP handles.ElectParams.filtering.LP],'band');
end


close(f)
set(handles.availability,'String','Free');
guidata(hObject, handles);


% --- Executes on button press in plot_spectro.
function plot_spectro_Callback(hObject, eventdata, handles)
% hObject    handle to plot_spectro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = guidata(hObject);
set(handles.availability,'String','Busy');


% spectrogram calculation
% params for power calculation 
params.fpass = [0.5 100]; % in Hz (0.5 to 10) - 0.5 is minimum if we want to use a 4 sec window and see one cycle
params.Fs = handles.current.SR;
k = 3; % num tapers (lower for less frequency leakage, higher for more leakage but smoother spectrum) % and it is a tradeoff with how many time points you have. spectrum is less smooth with many time points.
nw = (k+1)/2;
params.tapers=[nw k];

S={};
f={};
S_smooth={};

current_plot =handles.current.LFP;

fw = waitbar(0,'Please wait...');
for j = 1:length(current_plot(:,1)) % selection des structures
    waitbar(j/length(current_plot(:,1)),fw,' Processing ...');
    [S{1,j},f{1,j}]=mtspectrumc(current_plot(j,:),params); %must be samples by channels
    SS{1,j} = 10*log10(S{1,j});
    S_smoothS{1,j} = fastsmooth(SS{1,j},1000,3,1);
end

close(fw)


% Plot the spectrogram
cla(handles.axes5)
axes(handles.axes5)
plot(f{1,1},S_smoothS{1,1},'Color',rgb('DarkBlue')); hold on
title(['Spectro of ', cell2mat(handles.struct.current) ,' ',handles.current.source]);  ylabel('Power (db)'); xlabel('Freq (Hz)'); 

cla(handles.axes6)
axes(handles.axes6)
plot(f{1,2},S_smoothS{1,2},'Color',rgb('FireBrick')); hold on
ylabel('Power (db)'); xlabel('Freq (Hz)');

cla(handles.axes7)
axes(handles.axes7)
plot(f{1,3},S_smoothS{1,3},'Color',rgb('DarkViolet')); hold on
ylabel('Power (db)'); xlabel('Freq (Hz)');


set(handles.availability,'String','Free');
guidata(hObject, handles);




% --- Executes on button press in downsampling.
function downsampling_Callback(hObject, eventdata, handles)
% hObject    handle to downsampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.availability,'String','Busy');


F = 5; % downsampling factor

handles.ElectParams.New_SR = handles.ElectParams.SR/F; % in hz
handles.ElectParams.New_NR = handles.ElectParams.New_SR/2; % in hz
handles.LFP_downsampled = zeros(handles.ElectParams.Nchans,ceil(handles.ElectParams.Nsamples/F));

fw = waitbar(0,'Please wait...');

for k=1:handles.ElectParams.Nchans
    waitbar(k/handles.ElectParams.Nchans,fw,' Processing ...');
    handles.LFP_downsampled(k,:) = decimate(handles.LFP_filtered(k,:),F); % Fonction for downsampling 
end
close(fw)


set(handles.availability,'String','Free');
guidata(hObject, handles);


% --- Executes on selection change in second_source_select.
function second_source_select_Callback(hObject, eventdata, handles)
% hObject    handle to second_source_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns second_source_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from second_source_select
handles = guidata(hObject);


switch get(handles.second_source_select,'Value')
    case 1
        handles.current.source2 = 'raw LFP';
        handles.current.SR2 = handles.ElectParams.SR;
        handles.current.LFP2 = handles.LFP(handles.current.lock,:);
        % time calculation
        dt = 1/handles.current.SR2;
        tottime = length(handles.current.LFP2(1,:))/handles.current.SR2;
        handles.current.time2 = 0:dt:tottime-dt;
        
        disp('you chose raw data for double plot')
        
    case 2
        handles.current.source2 = 'filtered LFP';
        handles.current.SR2 = handles.ElectParams.SR;
        handles.current.LFP2 = handles.LFP_filtered(handles.current.lock,:);
        % time calculation
        dt = 1/handles.current.SR2;
        tottime = length(handles.current.LFP2(1,:))/handles.current.SR2;
        handles.current.time2 = 0:dt:tottime-dt;
        disp('you chose filtered data for double plot')

    case 3
        handles.current.source2 = 'downsampled LFP';
        handles.current.SR2 = handles.ElectParams.New_SR;
        handles.current.LFP2 = handles.LFP_downsampled(handles.current.lock,:);
        % time calculation
        dt = 1/handles.current.SR2;
        tottime = length(handles.current.LFP2(1,:))/handles.current.SR2;
        handles.current.time2 = 0:dt:tottime-dt; 
        disp('you chose downsampling data for double plot')
  
    otherwise
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function second_source_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to second_source_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in double_plot.
function double_plot_Callback(hObject, eventdata, handles)
% hObject    handle to double_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = guidata(hObject);

cla(handles.axes1)
axes(handles.axes1)
plot(handles.current.time,handles.current.LFP(1,:),'Color',rgb('DarkBlue'),'DisplayName',handles.current.source); hold on;
plot(handles.current.time2,handles.current.LFP2(1,:),'Color',rgb('LightBlue'),'DisplayName',handles.current.source2)
legend; title(['Channel n°1 of ',cell2mat(handles.struct.current)]); ylim([-handles.scale handles.scale]); hold off
zoom on

cla(handles.axes2)
axes(handles.axes2)
plot(handles.current.time,handles.current.LFP(2,:),'Color',rgb('FireBrick'),'DisplayName',handles.current.source); hold on;
plot(handles.current.time2,handles.current.LFP2(2,:),'Color',rgb('Salmon'),'DisplayName',handles.current.source2)
legend;
title(['Channel n°2 of ',cell2mat(handles.struct.current)]); ylim([-handles.scale handles.scale]); hold off
zoom on

cla(handles.axes3)
axes(handles.axes3)
plot(handles.current.time,handles.current.LFP(3,:),'Color',rgb('DarkViolet'),'DisplayName',handles.current.source); hold on;
plot(handles.current.time2,handles.current.LFP2(3,:),'Color',rgb('Thistle'),'DisplayName',handles.current.source2)
legend;
title(['Channel n°3 of ',cell2mat(handles.struct.current)]); ylim([-handles.scale handles.scale]); hold off
zoom on

guidata(hObject, handles);


% --- Executes on button press in plot_double_spectro.
function plot_double_spectro_Callback(hObject, eventdata, handles)
% hObject    handle to plot_double_spectro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.availability,'String','Busy');


% spectrogram calculation
% params for power calculation 
params.fpass = [0.5 100]; % in Hz (0.5 to 10) - 0.5 is minimum if we want to use a 4 sec window and see one cycle
params.Fs = handles.current.SR;
k = 3; % num tapers (lower for less frequency leakage, higher for more leakage but smoother spectrum) % and it is a tradeoff with how many time points you have. spectrum is less smooth with many time points.
nw = (k+1)/2;
params.tapers=[nw k];

S={};
f={};
S_smooth={};

current_plot =handles.current.LFP;

fw = waitbar(0,'Please wait...');
for j = 1:length(current_plot(:,1)) % selection des structures
    waitbar(j/6,fw,' Processing ...');
    [S{1,j},f{1,j}]=mtspectrumc(current_plot(j,:),params); %must be samples by channels
    SS{1,j} = 10*log10(S{1,j});
    S_smoothS{1,j} = fastsmooth(SS{1,j},1000,3,1);
end

S2={};
f2={};
S_smooth2={};
current_plot2 =handles.current.LFP2;
params.Fs = handles.current.SR2;

for j = 1:length(current_plot2(:,1)) % selection des structures
    waitbar(j+3/6,fw,' Processing ...');
    [S2{1,j},f2{1,j}]=mtspectrumc(current_plot2(j,:),params); %must be samples by channels
    SS2{1,j} = 10*log10(S2{1,j});
    S_smoothS2{1,j} = fastsmooth(SS2{1,j},1000,3,1);
end

close(fw)


% Plot the spectrogram
cla(handles.axes5)
axes(handles.axes5)
plot(f{1,1},S_smoothS{1,1},'Color',rgb('DarkBlue'),'DisplayName',handles.current.source); hold on
plot(f2{1,1},S_smoothS2{1,1},'Color',rgb('LightBlue'),'DisplayName',handles.current.source2);
legend;
title(['Spectrogram of ',cell2mat(handles.struct.current)]);  ylabel('Power (db)'); xlabel('Freq (Hz)'); 

cla(handles.axes6)
axes(handles.axes6)
plot(f{1,2},S_smoothS{1,2},'Color',rgb('FireBrick'),'DisplayName',handles.current.source); hold on
plot(f2{1,2},S_smoothS2{1,2},'Color',rgb('Salmon'),'DisplayName',handles.current.source2);
legend;
ylabel('Power (db)'); xlabel('Freq (Hz)');

cla(handles.axes7)
axes(handles.axes7)
plot(f{1,3},S_smoothS{1,3},'Color',rgb('DarkViolet'),'DisplayName',handles.current.source); hold on
plot(f2{1,3},S_smoothS2{1,3},'Color',rgb('Thistle'),'DisplayName',handles.current.source2);
legend;
ylabel('Power (db)'); xlabel('Freq (Hz)');


set(handles.availability,'String','Free');
guidata(hObject, handles);


function select_ch2_localref_Callback(hObject, eventdata, handles)
% hObject    handle to select_ch2_localref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of select_ch2_localref as text
%        str2double(get(hObject,'String')) returns contents of select_ch2_localref as a double
handles = guidata(hObject);
handles.localref.ch_21 = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function select_ch2_localref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_ch2_localref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.localref.ch_2 = str2double(get(hObject,'String'));
guidata(hObject, handles);



function select_ch_1_localref_Callback(hObject, eventdata, handles)
% hObject    handle to select_ch_1_localref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of select_ch_1_localref as text
%        str2double(get(hObject,'String')) returns contents of select_ch_1_localref as a double
handles = guidata(hObject);
handles.localref.ch_1 = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function select_ch_1_localref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_ch_1_localref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles = guidata(hObject);
handles.localref.ch_1 = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes on button press in plot_localref.
function plot_localref_Callback(hObject, eventdata, handles)
% hObject    handle to plot_localref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.availability,'String','Busy');

handles.localref.ch_1 = str2double(get(handles.select_ch_1_localref,'String'));
handles.localref.ch_2 = str2double(get(handles.select_ch2_localref,'String'));

handles.local.LFP1 = handles.LFP_downsampled(handles.current.lock(handles.localref.ch_1),:);
handles.local.LFP2 = handles.LFP_downsampled(handles.current.lock(handles.localref.ch_2),:);
handles.local.LFP3 = handles.local.LFP1 - handles.local.LFP2;

dt = 1/handles.ElectParams.New_SR;
tottime = length(handles.local.LFP1(1,:))/handles.ElectParams.New_SR;
time = 0:dt:tottime-dt;

cla(handles.axes4)
axes(handles.axes4)
plot(time,handles.local.LFP1,'Color',rgb('PaleGreen'),'DisplayName',['Ch n°',num2str(handles.localref.ch_1)]); hold on;
plot(time,handles.local.LFP2,'Color',rgb('LightBlue'),'DisplayName',['Ch n°',num2str(handles.localref.ch_2)]);
plot(time,handles.local.LFP3,'Color',rgb('LightSalmon'),'DisplayName',['Local Signal (ch n°',num2str(handles.localref.ch_1),'-',num2str(handles.localref.ch_2),')']);
legend; title(['Comparaison for local ref ',cell2mat(handles.struct.current)]); ylim([-handles.scale handles.scale]); hold off
zoom on

set(handles.availability,'String','Free');
guidata(hObject, handles);



% --- Executes on button press in plot_spectro_localref.
function plot_spectro_localref_Callback(hObject, eventdata, handles)
% hObject    handle to plot_spectro_localref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.availability,'String','Busy');


% spectrogram calculation
% params for power calculation 
params.fpass = [0.5 100]; % in Hz (0.5 to 10) - 0.5 is minimum if we want to use a 4 sec window and see one cycle
params.Fs = handles.ElectParams.New_SR;
k = 3; % num tapers (lower for less frequency leakage, higher for more leakage but smoother spectrum) % and it is a tradeoff with how many time points you have. spectrum is less smooth with many time points.
nw = (k+1)/2;
params.tapers=[nw k];

S={};
f={};
S_smooth={};

current_plot =handles.local.LFP1;

fw = waitbar(0,'Please wait...');
for j = 1:length(current_plot(:,1)) % selection des structures
    waitbar(j/9,fw,' Processing ...');
    [S{1,j},f{1,j}]=mtspectrumc(current_plot(j,:),params); %must be samples by channels
    SS{1,j} = 10*log10(S{1,j});
    S_smoothS{1,j} = fastsmooth(SS{1,j},1000,3,1);
end

S2={};
f2={};
S_smooth2={};
current_plot2 =handles.local.LFP2;

for j = 1:length(current_plot2(:,1)) % selection des structures
    waitbar(j+3/9,fw,' Processing ...');
    [S2{1,j},f2{1,j}]=mtspectrumc(current_plot2(j,:),params); %must be samples by channels
    SS2{1,j} = 10*log10(S2{1,j});
    S_smoothS2{1,j} = fastsmooth(SS2{1,j},1000,3,1);
end

S3={};
f3={};
S_smooth3={};
current_plot3 =handles.local.LFP3;

for j = 1:length(current_plot3(:,1)) % selection des structures
    waitbar(j+6/9,fw,' Processing ...');
    [S3{1,j},f3{1,j}]=mtspectrumc(current_plot3(j,:),params); %must be samples by channels
    SS3{1,j} = 10*log10(S3{1,j});
    S_smoothS3{1,j} = fastsmooth(SS3{1,j},1000,3,1);
end

close(fw)


% Plot the spectrogram
cla(handles.axes8)
axes(handles.axes8)
plot(f{1,1},S_smoothS{1,1},'Color',rgb('PaleGreen'),'DisplayName',['Ch n°',num2str(handles.localref.ch_1)]); hold on
plot(f2{1,1},S_smoothS2{1,1},'Color',rgb('LightBlue'),'DisplayName',['Ch n°',num2str(handles.localref.ch_2)]);
plot(f3{1,1},S_smoothS3{1,1},'Color',rgb('LightSalmon'),'DisplayName',['Local Signal (ch n°',num2str(handles.localref.ch_1),'-',num2str(handles.localref.ch_2),')']);
legend;
title(['Spectro of ',cell2mat(handles.struct.current)]);  ylabel('Power (db)'); xlabel('Freq (Hz)'); 


set(handles.availability,'String','Free');
guidata(hObject, handles);


function edit_save_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_save_path as text
%        str2double(get(hObject,'String')) returns contents of edit_save_path as a double
handles = guidata(hObject);
handles.pathSav = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_save_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in open_windows.
function open_windows_Callback(hObject, eventdata, handles)
% hObject    handle to open_windows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Where are the data

[~,handles.pathSav] = uiputfile('Select where to save the data');
set(handles.edit_save_path,'String',handles.pathSav);

% --- Executes on button press in save_local_ref.
function save_local_ref_Callback(hObject, eventdata, handles)
% hObject    handle to save_local_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.availability,'String','Busy');

namestruct = handles.struct.current{1,1};
handles.struct.cat = categorical(handles.struct.list(1,:));
lock= find(handles.struct.cat == namestruct);
handles.struct.list{2,lock} = handles.local.LFP3;

% find the availability for the ch
ch1_localref_aval = handles.ChParams.RatData(handles.current.lock(handles.localref.ch_1));
ch2_localref_aval = handles.ChParams.RatData(handles.current.lock(handles.localref.ch_2));

% store the availability
meanch= mean([ch1_localref_aval ch2_localref_aval]);
handles.struct.list{3,lock} = meanch;

% stored the ch used for local ref
handles.struct.list{4,lock} = handles.localref.ch_1;
handles.struct.list{5,lock} = handles.localref.ch_2;

% store temporally the data
LFP_local = handles.struct.list;
ChParams = handles.ChParams;
ElectParams = handles.ElectParams;
GenParams = handles.GenParams;
pathSav = handles.pathSav;
nameFolder = handles.nameFolder;
ChParams.localStruc = categorical({'LHb';'dHPC';'BLA';'ACC';'PRL';'EMG'});
save(fullfile(handles.pathSav,'temp_local_LFP.mat'),'ChParams' , 'ElectParams', 'GenParams', 'LFP_local', 'nameFolder', 'pathSav' );


% save figure 
name = ['preprocess_resume',cell2mat(handles.struct.current),'.jpeg'];
outputFileName = fullfile(handles.pathSav,name);


% do figure and save it in JPEG
fignew = figure; 
figure(fignew)

newAxes1 = copyobj(handles.axes1,fignew); % Copy the appropriate axes
newAxes2 = copyobj(handles.axes2,fignew); % Copy the appropriate axes
newAxes3 = copyobj(handles.axes3,fignew); % Copy the appropriate axes
newAxes4 = copyobj(handles.axes4,fignew); % Copy the appropriate axes

newAxes5 = copyobj(handles.axes5,fignew); % Copy the appropriate axes
newAxes6 = copyobj(handles.axes6,fignew); % Copy the appropriate axes
newAxes7 = copyobj(handles.axes7,fignew); % Copy the appropriate axes
newAxes8 = copyobj(handles.axes8,fignew); % Copy the appropriate axes


set(newAxes1,'Position',get(handles.axes1,'InnerPosition')./[2 1 1 1]); 
set(newAxes2,'Position',get(handles.axes2,'InnerPosition')./[2 1 1 1]);
set(newAxes3,'Position',get(handles.axes3,'InnerPosition')./[2 1 1 1]);
set(newAxes4,'Position',get(handles.axes4,'InnerPosition')./[2 1 1 1]);

set(newAxes5,'Position',get(handles.axes5,'InnerPosition')); 
set(newAxes6,'Position',get(handles.axes6,'InnerPosition'));
set(newAxes7,'Position',get(handles.axes7,'InnerPosition'));
set(newAxes8,'Position',get(handles.axes8,'InnerPosition'));
set(gcf, 'Position', get(0, 'Screensize')); % full screen

saveas(fignew,outputFileName);
disp('save big image done')
close(fignew)
set(handles.availability,'String','Free');
guidata(hObject, handles);


% --- Executes on button press in cutting_signal.
function cutting_signal_Callback(hObject, eventdata, handles)
% hObject    handle to cutting_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
set(handles.availability,'String','Busy');

% order the Local signals
handles.LFPlocal = [handles.struct.list{2,1};handles.struct.list{2,2};handles.struct.list{2,3};handles.struct.list{2,4};handles.struct.list{2,5};handles.LFP_downsampled(4,:)];
handles.ChParams.localParams = [handles.struct.list{3,1}; handles.struct.list{3,2}; handles.struct.list{3,3}; handles.struct.list{3,4}; handles.struct.list{3,5}; handles.ChParams.RatData(4,1)];

% cut signal
if handles.GenParams.SubSession == categorical({'S1'}) || handles.GenParams.SubSession == categorical({'S2'})
    Size = length(handles.LFPlocal(1,:))/handles.ElectParams.New_SR/60; % en minute
    message = ['The signal is ',num2str(Size),' minutes long, no need to cut it by hours'];
    disp(message)
    
    handles.DatabyHours = handles.LFPlocal;
else
    CuttingTimes = 3600; % in s
    % Local
    [DatabyHoursLocal] = slicingSignal(handles.LFPlocal,handles.ElectParams.New_SR,CuttingTimes);
    % Not local
    [handles.DatabyHours] = slicingSignal(handles.LFP_downsampled,handles.ElectParams.New_SR,CuttingTimes);
    
    % put the cell arrays together
    handles.DatabyHours(3,:) = DatabyHoursLocal(2,:);
end


% save the variable I need 
ChParams = handles.ChParams;
ElectParams = handles.ElectParams;
GenParams = handles.GenParams;
DatabyHours = handles.DatabyHours;
LFP_local = handles.struct.list;
pathSav = handles.pathSav;
nameFolder = handles.nameFolder;
ChParams.localStruc = categorical({'LHb';'dHPC';'BLA';'ACC';'PRL';'EMG'});


nameData = ['DataPrepro',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.mat'];
save(fullfile(pathSav,nameData),'ChParams' , 'DatabyHours', 'ElectParams', 'GenParams', 'LFP_local', 'nameFolder', 'pathSav' )

% don't forget that now the EMG is in 6 st pposition  - modify iton the
% other scripts.


set(handles.availability,'String','Free');
guidata(hObject, handles);




