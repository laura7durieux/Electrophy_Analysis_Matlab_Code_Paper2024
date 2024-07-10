function varargout = GUI_sleep_V3(varargin)
% GUI_SLEEP_V3 MATLAB code for GUI_sleep_V3.fig
%      GUI_SLEEP_V3, by itself, creates a new GUI_SLEEP_V3 or raises the existing
%      singleton*.
%
%      H = GUI_SLEEP_V3 returns the handle to a new GUI_SLEEP_V3 or the handle to
%      the existing singleton*.
%
%      GUI_SLEEP_V3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SLEEP_V3.M with the given input arguments.
%
%      GUI_SLEEP_V3('Property','Value',...) creates a new GUI_SLEEP_V3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_sleep_V3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_sleep_V3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_sleep_V3

% Last Modified by GUIDE v2.5 20-Apr-2020 10:16:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_sleep_V3_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_sleep_V3_OutputFcn, ...
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


% --- Executes just before GUI_sleep_V3 is made visible.
function GUI_sleep_V3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_sleep_V3 (see VARARGIN)

% Choose default command line output for GUI_sleep_V3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_sleep_V3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_sleep_V3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%% LOADING DATA %%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_load_data.
function pushbutton_load_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
% Load the data
[file1,handles.pathSav] = uigetfile('*.mat',...
   'Select the matlab file C with all the data for the spectrogram', ...
   'MultiSelect', 'off');
load(fullfile(handles.pathSav,file1))
load(fullfile(handles.pathSav,'managingPath.mat'))

handles.managingPath = managingPath;

handles.hours.number = length(PowdHPC(1,:));

set(handles.text5,'String',handles.hours.number);

handles.raw.PowdHPC = PowdHPC;
handles.raw.PowPRL = PowPRL;

handles.raw.EMGthres = EMGthres;
handles.raw.time_EMG = time_EMG;
handles.raw.DLC_Sleep_by_hours = DLC_Sleep_by_hours;
handles.raw.DatabyHours = DatabyHours;

handles.GenParams = GenParams;
handles.ElectParams = ElectParams;
handles.ChParams = ChParams;


% Attribuate the variable needed
handles.hours.selected = 1;

handles.scale.HPC = 0.0001;
handles.scale.PRL = 0.0001;

disp('data loading is done')

guidata(hObject, handles);

%%%%%%%%%%% SELECTION OF THE HOUR TO ANALYSE %%%%%%%%%%%%%%%    
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


handles = guidata(hObject);

handles.hours.selected = str2double(get(hObject,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%% LOADING THE HOUR TO ANALYSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

% % PowdHPC{2,i} = S; % convert to dB;
% % PowdHPC{3,i} = t;
% % PowdHPC{4,i} = f;

handles.HPC.S= handles.raw.PowdHPC{2,handles.hours.selected};
handles.PRL.S= handles.raw.PowPRL{2,handles.hours.selected};
handles.HPC.t= handles.raw.PowdHPC{3,handles.hours.selected};
handles.PRL.t= handles.raw.PowPRL{3,handles.hours.selected};
handles.HPC.f= handles.raw.PowdHPC{4,handles.hours.selected};
handles.PRL.f= handles.raw.PowPRL{4,handles.hours.selected};
handles.HPC.ratio= handles.raw.PowdHPC{5,handles.hours.selected};
handles.PRL.ratio= handles.raw.PowPRL{5,handles.hours.selected};

if isnan(handles.raw.EMGthres{2,1})== 0
    handles.mvt.EMG_bin = handles.raw.EMGthres{3,handles.hours.selected};
    handles.mvt.EMG = handles.raw.EMGthres{2,handles.hours.selected};
end
handles.mvt.EMG_time = handles.raw.time_EMG;
handles.mvt.DLC = handles.raw.DLC_Sleep_by_hours{2,handles.hours.selected}(1,:);
handles.mvt.DLCtime = handles.raw.DLC_Sleep_by_hours{2,1}(2,:);

handles.raw.LFP = handles.raw.DatabyHours{2,handles.hours.selected};
handles.raw.LFP_local = handles.raw.DatabyHours{3,handles.hours.selected};
handles.raw.srate = handles.ElectParams.New_SR;


handles.nameFolder = ['stage_extraction_hour',num2str(handles.hours.selected)];
mkdir(handles.pathSav,handles.nameFolder) % folder creation

disp('Loading the hour to analyse is done')

guidata(hObject, handles);


%%%%%%%%%%%% PLOT SPECTRO FOR SLEEP ASSESSMENT %%%%%%%%%%%%%%%%%%%%% 
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

cla(handles.axes1)
axes(handles.axes1)
plotsig(handles.HPC.S,0,handles.HPC.t,handles.HPC.f); axis xy; colorbar;colormap cool;title('Dorsal Hippocampus')
ylim([0 18])
xlim([0 3600])
caxis([0, handles.scale.HPC]);

cla(handles.axes2)
axes(handles.axes2)
plotsig(handles.PRL.S,0,handles.PRL.t,handles.PRL.f); axis xy; colorbar;colormap cool;title('Prelimbic Cortex')
ylim([0 18])
xlim([0 3600])
caxis([0, handles.scale.PRL]);

cla(handles.axes3)
axes(handles.axes3)
yyaxis left
plot(handles.HPC.t,handles.HPC.ratio,'DisplayName', 'theta delta ratio from HPC')
hold on
yyaxis right
plot(handles.PRL.t,handles.PRL.ratio,'DisplayName', 'delta gamma ratio from PRL')
legend
xlim([0 3750])
title('Power ratio')
hold off

guidata(hObject, handles);


%%%%%%%%%%%%%%%%% PLOT THE MOUVEMENT EMG AND DLC %%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles = guidata(hObject);
switch get(handles.popupmenu1,'Value')
    case 1
        cla(handles.axes4)
        axes(handles.axes4)
        hold on
        plot(handles.mvt.DLCtime,-handles.mvt.DLC,'r','DisplayName','DLC')
        plot(handles.mvt.EMG_time,handles.mvt.EMG_bin,'c','DisplayName','EMG bin')
        plot(handles.mvt.EMG_time,handles.mvt.EMG,'k','DisplayName','EMG')
        legend
        ylim([-2 2])
        xlim([0 3750])
        title('Movements')
        hold off
    case 2
        cla(handles.axes4)
        axes(handles.axes4)
        hold on
        plot(handles.mvt.DLCtime,-handles.mvt.DLC,'r','DisplayName','DLC')
        legend
        ylim([-2 2])
        xlim([0 3750])
        title('Movements')
        hold off
    case 3
        cla(handles.axes4)
        axes(handles.axes4)

        hold on
        plot(handles.mvt.EMG_time,handles.mvt.EMG_bin,'c','DisplayName','EMG bin')
        plot(handles.mvt.EMG_time,handles.mvt.EMG,'k','DisplayName','EMG')
        legend
        ylim([-2 2])
        xlim([0 3750])
        title('Movements')
        hold off
    otherwise
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%% SET HPC SCALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

handles = guidata(hObject);

handles.scale.HPC = str2double(get(hObject,'String'));

guidata(hObject, handles);




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


%%%%%%%%%%%%%% SET PRL SCALE %%%%%%%%%%%%%%%%%%%%%%%%
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles = guidata(hObject);

handles.scale.PRL = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%% STAGES ASSESSMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
% Finding the x from the plot
[x,y,button] = ginput(100); % return Key to stoping the function

% Create sleepTrack
handles.sleepTrack = Position_calculation(x,y,button,(max(handles.HPC.t)));


handles.fSpect = figure;
figure(handles.fSpect)
subplot(211)
hold on
yyaxis left
plotsig(handles.HPC.S,0,handles.HPC.t,handles.HPC.f); axis xy; colorbar;colormap jet;title('Dorsal Hippocampus')
ylim([0 18])
xlim([0 3600])
caxis([0, handles.scale.HPC]);
yyaxis right
plot(handles.sleepTrack,'r','LineWidth',3)
xlim([0 3600]); ylim([0 5]); ylabel('hypnogram'); yticks([1 2 3 4]); yticklabels({'REM','SWS','CW','AW'})
hold off

subplot(212)
hold on
yyaxis left
plotsig(handles.PRL.S,0,handles.PRL.t,handles.PRL.f); axis xy; colorbar;colormap jet;title('Prelimbique Cortex')
ylim([0 18])
xlim([0 3600])
caxis([0, handles.scale.PRL]);
yyaxis right
plot(handles.sleepTrack,'r','LineWidth',3)
xlim([0 3600]); ylim([0 5]); ylabel('hypnogram'); yticks([1 2 3 4]); yticklabels({'REM','SWS','CW','AW'})
hold off

% save figure
savefig(handles.fSpect,fullfile(handles.pathSav,handles.nameFolder,'hypnogram.fig'));

guidata(hObject, handles);


%%%%%%%%%%%%%% SEPARATED THE STAGES DEPENDING OF THE SLEEPTRACK %%%%%%%%
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

[handles.results.raw.AW,handles.results.raw.CW,handles.results.raw.SWS,handles.results.raw.REM] = Separate_statesV2(handles.raw.LFP,handles.sleepTrack,handles.raw.srate);
[handles.results.local.AW,handles.results.local.CW,handles.results.local.SWS,handles.results.local.REM] = Separate_statesV2(handles.raw.LFP_local,handles.sleepTrack,handles.raw.srate);
x=[num2str(size(handles.results.raw.AW,1)),' episode for AW, ', num2str(size(handles.results.raw.CW,1)),' episode for CW, ',num2str(size(handles.results.raw.SWS,1)),' episode for SWS, ',num2str(size(handles.results.raw.REM,1)),' episode for REM'];
disp(x)

guidata(hObject, handles);

%%%%%%%%% STAGES ANALYSIS #LATENCY, OCCURENCIES, TIME %%%%%%%%%%%%%
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

Stage= {handles.results.local.AW; handles.results.local.CW; handles.results.local.SWS; handles.results.local.REM};
[handles.results.SleepSettings,handles.results.SleepSettingsTable] = SleepAnalysis(Stage,handles.sleepTrack,handles.raw.srate);

handles.sleepfig = figure;

figure(handles.sleepfig)
subplot(131)
bar(handles.results.SleepSettingsTable.StartLantencyS,"b")
set(gca,'xticklabel',handles.results.SleepSettingsTable.Stages)
title("Latency")
subplot(132)
bar(handles.results.SleepSettingsTable.OccurenciesNb,"r")
set(gca,'xticklabel',handles.results.SleepSettingsTable.Stages)
title("Occurencies")
subplot(133)
bar(handles.results.SleepSettingsTable.TimeSpendS,"g")
set(gca,'xticklabel',handles.results.SleepSettingsTable.Stages)
title("Total time spend")

% Save the data excel file
writetable(handles.results.SleepSettingsTable,fullfile(handles.pathSav,handles.nameFolder,'sleepAnalysis_results.xlsx'),'Sheet',1)

guidata(hObject, handles);


%%%%%% CHECKING THE FREQUENCIES OF EACH PERIODE %%%%%%%%%%%%%%%%
% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

handles = guidata(hObject);

        % params for power calculation 
        params.fpass = [1 100]; % in Hz (0.5 to 10) - 0.5 is minimum if we want to use a 4 sec window and see one cycle
        params.Fs = handles.raw.srate;
        k = 3; % num tapers (lower for less frequency leakage, higher for more leakage but smoother spectrum) % and it is a tradeoff with how many time points you have. spectrum is less smooth with many time points.
        nw = (k+1)/2;
        params.tapers=[nw k];
        
        

switch get(handles.popupmenu2,'Value')
    case 1 % for AW    
        handles.stages = 'AW';
        disp(handles.stages)
        periodes = size(handles.results.local.AW,1);
        S={};
        f={};
        S_smooth={};
        
        handles.Spectro = figure;
        
        for i = 1:periodes
            for j = [1,3,4,5,6]
                [S{i,1}(:,j),f{i,1}(j,:)]=mtspectrumc(handles.results.local.AW{i,1}(j,:),params); %must be samples by channels
                SS{i,1}(:,j) = 10*log10(S{i,1}(:,j));
                % sprintf('your freq bins are %2.4f Hz',unique(diff(f{i,1}(j,:))))
                S_smoothS{i,1}(:,j) = fastsmooth(SS{i,1}(:,j),1000,3,1);
            end
   
           figure(handles.Spectro);
           subplot(231);
           plot(f{i,1}(1,:),S_smoothS{i,1}(:,1),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for LHb');  ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(232);
           plot(f{i,1}(3,:),S_smoothS{i,1}(:,3),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for dHPC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(233);
           plot(f{i,1}(4,:),S_smoothS{i,1}(:,4),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for BLA'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(234);
           plot(f{i,1}(5,:),S_smoothS{i,1}(:,5),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for ACC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(235);
           plot(f{i,1}(6,:),S_smoothS{i,1}(:,6),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for PRL'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           legend;  
        end

    case 2
        handles.stages = 'CW';
                disp(handles.stages)

        periodes = size(handles.results.local.CW,1);
        S={};
        f={};
        S_smooth={};
        
        handles.Spectro = figure;
        
        for i = 1:periodes
            for j = [1,3,4,5,6]
                [S{i,1}(:,j),f{i,1}(j,:)]=mtspectrumc(handles.results.local.CW{i,1}(j,:),params); %must be samples by channels
                SS{i,1}(:,j) = 10*log10(S{i,1}(:,j));
%                 sprintf('your freq bins are %2.4f Hz',unique(diff(f{i,1}(j,:))))
                S_smoothS{i,1}(:,j) = fastsmooth(SS{i,1}(:,j),1000,3,1);
            end
   
           figure(handles.Spectro);
           subplot(231);
           plot(f{i,1}(1,:),S_smoothS{i,1}(:,1),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for LHb');  ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(232);
           plot(f{i,1}(3,:),S_smoothS{i,1}(:,3),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for dHPC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(233);
           plot(f{i,1}(4,:),S_smoothS{i,1}(:,4),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for BLA'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(234);
           plot(f{i,1}(5,:),S_smoothS{i,1}(:,5),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for ACC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(235);
           plot(f{i,1}(6,:),S_smoothS{i,1}(:,6),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for PRL'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           legend   ;      
        end
    case 3
        handles.stages = 'SWS';
                disp(handles.stages)

        periodes = size(handles.results.local.SWS,1);
        S={};
        f={};
        S_smooth={};
        
        handles.Spectro = figure;
        
        for i = 1:periodes
            for j = [1,3,4,5,6]
                [S{i,1}(:,j),f{i,1}(j,:)]=mtspectrumc(handles.results.local.SWS{i,1}(j,:),params); %must be samples by channels
                SS{i,1}(:,j) = 10*log10(S{i,1}(:,j));
%                 sprintf('your freq bins are %2.4f Hz',unique(diff(f{i,1}(j,:))))
                S_smoothS{i,1}(:,j) = fastsmooth(SS{i,1}(:,j),1000,3,1);
            end
   
           figure(handles.Spectro);
           subplot(231);
           plot(f{i,1}(1,:),S_smoothS{i,1}(:,1),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for LHb');  ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(232);
           plot(f{i,1}(3,:),S_smoothS{i,1}(:,3),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for dHPC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(233);
           plot(f{i,1}(4,:),S_smoothS{i,1}(:,4),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for BLA'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(234);
           plot(f{i,1}(5,:),S_smoothS{i,1}(:,5),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for ACC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(235);
           plot(f{i,1}(6,:),S_smoothS{i,1}(:,6),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for PRL'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           legend;         
        end
        
    case 4
        handles.stages = 'REM';
                disp(handles.stages)

        periodes = size(handles.results.local.REM,1);
        S={};
        f={};
        S_smooth={};
        
        handles.Spectro = figure;
        
        for i = 1:periodes
            for j = [1,3,4,5,6] % selection des structures
                [S{i,1}(:,j),f{i,1}(j,:)]=mtspectrumc(handles.results.local.REM{i,1}(j,:),params); %must be samples by channels
                SS{i,1}(:,j) = 10*log10(S{i,1}(:,j));
%                 sprintf('your freq bins are %2.4f Hz',unique(diff(f{i,1}(j,:))))
                S_smoothS{i,1}(:,j) = fastsmooth(SS{i,1}(:,j),1000,3,1);
            end
   
           figure(handles.Spectro);
           subplot(231);
           plot(f{i,1}(1,:),S_smoothS{i,1}(:,1),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for LHb');  ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(232);
           plot(f{i,1}(3,:),S_smoothS{i,1}(:,3),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for dHPC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(233);
           plot(f{i,1}(4,:),S_smoothS{i,1}(:,4),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for BLA'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(234);
           plot(f{i,1}(5,:),S_smoothS{i,1}(:,5),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for ACC'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           subplot(235);
           plot(f{i,1}(6,:),S_smoothS{i,1}(:,6),'color',rand(1,3),'DisplayName',['period n°', num2str(i)]); hold on
           title('Spectrogram for PRL'); ylabel('Power (db)'); xlabel('Freq (Hz)'); legend;
           legend ;        
        end
    otherwise
end


guidata(hObject, handles);


%%%%%%%%%% SAVE THE SPECTROGRAME OF NEEDED %%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% saving the spectrogram plot
handles = guidata(hObject);

namefig = [char(handles.stages), '_Spectrogram.fig'];
savefig(handles.Spectro,fullfile(handles.pathSav,handles.nameFolder,namefig));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%% SAVING THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Variable management
SleepTrack = handles.sleepTrack;
Results = handles.results;
GenParams = handles.GenParams;
ElectParams = handles.ElectParams;
ChParams = handles.ChParams;


nameData = ['Stages_assessed_hours',num2str(handles.hours.selected),'.mat'];
handles.managingPath{handles.hours.selected+1,:} = fullfile(handles.pathSav,handles.nameFolder,nameData);

% Save the data
save(fullfile(handles.pathSav,handles.nameFolder,nameData),'SleepTrack', 'Results', 'GenParams','ElectParams','ChParams');

managingPath = handles.managingPath;
save((fullfile(handles.pathSav,'managingPath.mat')), 'managingPath');

disp(['Saved in ',fullfile(handles.pathSav,handles.nameFolder,nameData)]);

guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

load(fullfile(handles.pathSav,'managingPath.mat'));

[sleepsettings_total,sleepsettings_total_table] = sleep_analysis_total_time(managingPath);

handles.sleeptotfig = figure;

figure(handles.sleeptotfig)
subplot(131)
bar(sleepsettings_total_table.StartLantencyS,"b")
set(gca,'xticklabel',sleepsettings_total_table.Stages)
title("Latency")
subplot(132)
bar(sleepsettings_total_table.OccurenciesNb,"r")
set(gca,'xticklabel',sleepsettings_total_table.Stages)
title("Occurencies")
subplot(133)
bar(sleepsettings_total_table.TimeSpendS,"g")
set(gca,'xticklabel',sleepsettings_total_table.Stages)
title("Total time spend")

% save figure
savefig(handles.sleeptotfig,fullfile(handles.pathSav,'sleepAnalysis.fig'));

% save data in excel file
writetable(sleepsettings_total_table,fullfile(handles.pathSav,'sleepAnalysis_full.xlsx'),'Sheet',1)

% save into mat file 
save((fullfile(handles.pathSav,'Full_sleep_analysis.mat')), 'sleepsettings_total_table','sleepsettings_total');

disp([" Saved in : ", handles.pathSav])

guidata(hObject, handles);





