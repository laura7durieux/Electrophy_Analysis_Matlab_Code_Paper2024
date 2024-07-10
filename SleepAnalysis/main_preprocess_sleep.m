%% main script for preprocessing before entrering in the GUI for sleep assessement 
clear all
close all

%% Load the data (à mettre en auto dans une prochaine version?)
[file,pathSav] = uigetfile('*.mat',...
    'Select the matlab file that you want to analyze', ...
    'MultiSelect', 'off');
if isequal(file,0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(path,file)]);
end

load(fullfile(pathSav,file))
disp(fullfile(pathSav,file))


%% Select the data for sleep analysis (channels)

ChParams.dHPCExample = LFP_local{4,2}+4; % because the HPCd in raw data beging at the ch 5
ChParams.PRLExample = LFP_local{4,5}+13; % because the PRL in raw data beging at the ch 14


[DatabyHours] = select_ch_for_sleep(DatabyHours,ChParams.localStruc,ChParams.localParams,ChParams.dHPCExample,ChParams.PRLExample);

% plot all the ch for visual checking
checkraw = figure;
NbHour = length(DatabyHours(1,:));
struct = 2;

% Calculate a rep vector as 1 2 1 2 1 2 
repStruc = repmat(1:struct,1,NbHour);

% Calculate a rep vector as 1 1 2 2 3 3
k = [1;1]*(1:NbHour);
repData = k(:)';
 
% do the plot with a loop (depend of how many hours)
for i = 1:(NbHour*struct) %for number of subplot
    figure(checkraw) 
    subplot(NbHour,struct,i)
    plot(DatabyHours{4,repData(i)}(repStruc(i),:))
    x =['Plot of structure n°',num2str(repStruc(i)),' hour n°',num2str(repData(i))];
    title(x)
    ylabel('en mV')
    ylim([-0.5 0.5])
    xlabel('en second')
end

%% Power calculation

[PowdHPC,PowPRL]=power_presleep_chronus(DatabyHours,ElectParams.New_SR);

% checking how is the power for the 1st hour
spect = figure;
figure(spect)
subplot(211)
plotsig(PowdHPC{2,1},0,PowdHPC{3,1},PowdHPC{4,1}); axis xy; colorbar;colormap jet;title('Dorsal Hippocampus')
ylim([0 20])
xlim([0 3600])
caxis([0, 0.0001]);

subplot(212)
plotsig(PowPRL{2,1},0,PowPRL{3,1},PowPRL{4,1}); axis xy; colorbar;colormap jet;title('Prelimbique cortex')
ylim([0 20])
xlim([0 3600])
caxis([0, 0.0001]);

%% Ratio calculation 
% Caculate the power range for HPC

NbHour = length(DatabyHours(1,:));

for i =1:NbHour
    dHPC_S = PowdHPC{2,i};
    dHPC_f = PowdHPC{4,i};
    PRL_S = PowPRL{2,i};
    PRL_f = PowPRL{4,i};
    
    [theta_delta_ratio,delta_gamma_ratio] = ratio_for_sleep(dHPC_S,dHPC_f,PRL_S,PRL_f);
    
    PowdHPC{5,i} = theta_delta_ratio';
    PowPRL{5,i} = delta_gamma_ratio';
end

fratio = figure;
figure(fratio)
plot(PowdHPC{3,1},PowdHPC{5,1},'DisplayName', 'theta delta ratio from HPC')
hold on
plot(PowPRL{3,1},PowPRL{5,1},'DisplayName', 'delta gamma ratio from PRL')
legend
hold off

%% Mouvement scoring managmenet

% EMG Checking (need to be set for all hours)
EMGthres ={};
if mean(isnan(DatabyHours{5,1})) == 0 % not a NaN
    for i =1:NbHour  
        EMG = DatabyHours{5,i};
    
        % calculate mean +sd
        meanEMG = mean(EMG);
        sdEMG = 1.3*std(EMG);
        quantileEMG = quantile(EMG,0.9995);

        % set threshold
        LockValueIn = find(EMG<quantileEMG & EMG> (-quantileEMG));
        EMGbin = ones(1,length(EMG));
        EMGbin(1,LockValueIn)= 0;
        
        %output
        EMGthres{1,i} = i;
        EMGthres{2,i} = EMG;
        EMGthres{3,i} = EMGbin;
    end
    
    % plot (attention montre la troisième heure
    EMG_threshold = figure;
    
    figure(EMG_threshold)
    plot(EMGbin,'b', 'DisplayName','EMG bin')
    hold on
    plot(EMG,'c','DisplayName','EMG data')
    l1 = line([0 length(EMG)],[meanEMG meanEMG],'Color','red','LineWidth',1,'DisplayName','mean');
    l2=line([0 length(EMG)],[-quantileEMG -quantileEMG],'Color','black','LineWidth',1,'DisplayName',' minus sd');
    l3=line([0 length(EMG)],[quantileEMG quantileEMG],'Color','green','LineWidth',1,'DisplayName','plus sd');
    legend
    hold off
else
    disp('No EMG to threshold')
    for i =1:NbHour 
        EMG = DatabyHours{5,i};
        EMGthres{1,i} = i;
        EMGthres{2,i} = EMG;
    end  
end


%% Deeplabcup project
DLC_user = input('Do you want to enter DLC results ? [y/n]','s');

if DLC_user == 'y' 
    [file,path] = uigetfile('*.csv',...
        'Select the cvs file from DLC for sleep', ...
        'MultiSelect', 'off');
    if isequal(file,0)
        disp('User selected Cancel');
    else
        disp(['User selected ', fullfile(path,file)]);
    end



    fid = fopen(fullfile(path,file));
    DLC_Sleep = textscan(fid, '%f', 'Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
    fclose (fid);
    DLC_Sleep = DLC_Sleep{1,1}'; % convert table to array
    div_factor = length(DLC_Sleep)/2;
    DLC_Sleep = [DLC_Sleep(1,1:div_factor); DLC_Sleep(1,div_factor+1:end)];
    
    
    % checking the length of the signal
    lengthSignal = length(DLC_Sleep(1,:))/30; % time in sec
    if  GenParams.Session == "HAB1" | GenParams.Session == "HAB2" | GenParams.Session == "HAB3"
        if lengthSignal < 3600*3
            DLC_User2 = input('Do you want to enter a second DLC results (time is not egal to 3 hours) ? [y/n]','s');
            if DLC_User2 == 'y'
                [file,path] = uigetfile('*.csv',...
                'Select the cvs file from DLC for sleep', ...
                'MultiSelect', 'off');
                if isequal(file,0)
                    disp('User selected Cancel');
                else
                    disp(['User selected ', fullfile(path,file)]);
                end
                
                fid = fopen(fullfile(path,file));
                DLC_Sleep2 = textscan(fid, '%f', 'Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
                fclose (fid);
                DLC_Sleep2 = DLC_Sleep2{1,1}'; % convert table to array
                div_factor = length(DLC_Sleep2)/2;
                DLC_Sleep2 = [DLC_Sleep2(1,1:div_factor); DLC_Sleep2(1,div_factor+1:end)];
                
                DLC_Sleep = [DLC_Sleep DLC_Sleep2];
                
                %time calculation
                len_data = length(DLC_Sleep(1,:));
                tottime = len_data/30;
                dt = 1/30;
                time = 0:dt:tottime;
                DLC_Sleep(2,:)= time(1,1:end-1);
           else
                disp('you will not have the 3 hours')
           end
        end
    elseif (GenParams.Session == "Stress1" | GenParams.Session == "Stress2") & (GenParams.SubSession == "PS1" | GenParams.SubSession == "PS2")
       if lengthSignal < 3600*3
            DLC_User2 = input('Do you want to enter a second DLC results (time is not egal to 3 hours) ? [y/n]','s');
            if DLC_User2 == 'y'
                [file,path] = uigetfile('*.csv',...
                'Select the cvs file from DLC for sleep', ...
                'MultiSelect', 'off');
                if isequal(file,0)
                    disp('User selected Cancel');
                else
                    disp(['User selected ', fullfile(path,file)]);
                end
                
                fid = fopen(fullfile(path,file));
                DLC_Sleep2 = textscan(fid, '%f', 'Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
                fclose (fid);
                DLC_Sleep2 = DLC_Sleep2{1,1}'; % convert table to array
                div_factor = length(DLC_Sleep2)/2;
                DLC_Sleep2 = [DLC_Sleep2(1,1:div_factor); DLC_Sleep2(1,div_factor+1:end)];
                
                DLC_Sleep = [DLC_Sleep DLC_Sleep2];
                
                %time calculation
                len_data = length(DLC_Sleep(1,:));
                tottime = len_data/30;
                dt = 1/30;
                time = 0:dt:tottime;
                DLC_Sleep(2,:)= time(1,1:end-1);
           else
                disp('you will not have the 3 hours')
           end
       end
    elseif GenParams.Session == "StressBL"
       if lengthSignal < 3600*1
            DLC_User2 = input('Do you want to enter a second DLC results (time is not egal to 1 hours) ? [y/n]','s');
            if DLC_User2 == 'y'
                [file,path] = uigetfile('*.csv',...
                'Select the cvs file from DLC for sleep', ...
                'MultiSelect', 'off');
                if isequal(file,0)
                    disp('User selected Cancel');
                else
                    disp(['User selected ', fullfile(path,file)]);
                end
                
                fid = fopen(fullfile(path,file));
                DLC_Sleep2 = textscan(fid, '%f', 'Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
                fclose (fid);
                DLC_Sleep2 = DLC_Sleep2{1,1}'; % convert table to array
                div_factor = length(DLC_Sleep2)/2;
                DLC_Sleep2 = [DLC_Sleep2(1,1:div_factor); DLC_Sleep2(1,div_factor+1:end)];
                
                DLC_Sleep = [DLC_Sleep DLC_Sleep2];
                
                %time calculation
                len_data = length(DLC_Sleep(1,:));
                tottime = len_data/30;
                dt = 1/30;
                time = 0:dt:tottime;
                DLC_Sleep(2,:)= time(1,1:end-1);
           else
                disp('you will not have the 3 hours')
           end
       end
    else
        disp("you don't need to pass you data into the presleep process")
    end
                

    [DLC_Sleep_by_hours] = slicingSignal(DLC_Sleep,30,3600);

    figure 
    plot(DLC_Sleep_by_hours{2,1}(2,:),DLC_Sleep_by_hours{2,1}(1,:))
    ylim([-1 2])
    
    
    if mean(isnan(DatabyHours{5,1})) == 0 % not a NaN
        % calcul of the time 
        time_EMG = [0:1/ElectParams.New_SR:3600];
        time_EMG = time_EMG(:,1:end-1);
        negDLCH1 = -(DLC_Sleep_by_hours{2,1}(1,:));
        %plot
        fDLC=figure;
        figure(fDLC)
        plot(time_EMG,EMGthres{3,1},'b', 'DisplayName','EMG bin')
        hold on
        plot(DLC_Sleep_by_hours{2,1}(2,:),negDLCH1,'r','DisplayName','DLC')
        plot(time_EMG,DatabyHours{5,1},'c','DisplayName','EMG data')
        legend
        hold off
    elseif mean(isnan(DatabyHours{5,1})) == 1 % a NaN
        negDLCH1 = -(DLC_Sleep_by_hours{2,1}(1,:));
        fDLC=figure;
        figure(fDLC)
        plot(DLC_Sleep_by_hours{2,1}(2,:),negDLCH1,'r','DisplayName','DLC')
        
        % creating EMG (avoiding error after)
        time_EMG = [0:1/ElectParams.New_SR:3600];
        time_EMG = time_EMG(:,1:end-1);
    end
    
else
    disp("you didn't load DLC results")
    DLC_Sleep_by_hours = {};
    time_EMG = [0:1/ElectParams.New_SR:3600];
    time_EMG = time_EMG(:,1:end-1);
end


%% Save the signal in parts for the GUI (signal, threshold and ...)

namePlotA = ['checkingdata',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.fig'];
namePlotB = ['spectro',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.fig'];
namePlotC = ['ratioEMG',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.fig'];
namePlotD = ['EMG_threshold',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.fig'];
namePlotE = ['DLC_compare',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.fig'];

saveas(checkraw,fullfile(pathSav,namePlotA));
saveas(spect,fullfile(pathSav,namePlotB));
saveas(fratio,fullfile(pathSav,namePlotC));
saveas(EMG_threshold,fullfile(pathSav,namePlotD));
saveas(fDLC,fullfile(pathSav,namePlotE));

close all
clearvars -except ChParams DatabyHours ElectParams GenParams pathSav nameFolder DLC_Sleep_by_hours EMGthres PowdHPC PowPRL time_EMG

nameData = ['Pre_sleep',num2str(GenParams.Rat),char(GenParams.Session),char(GenParams.SubSession),'.mat'];
save(fullfile(pathSav,nameData))


managingPath{1,:} = fullfile(pathSav,nameFolder);
namePath= ['managingPath.mat'];
save(fullfile(pathSav,namePath),'managingPath')




