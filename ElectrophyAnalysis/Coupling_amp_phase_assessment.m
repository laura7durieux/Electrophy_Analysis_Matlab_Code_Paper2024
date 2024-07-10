%% Coupling analysis
% This script - coupling evaluation part - comes from Nelson Totah (see paper for atributions)
close all
clear

% set up variable for loading 
session_name = 'Stress2';
subsession_name = ['PS2'];
saving_folder= 'AW'
pathway_source = fullfile('D:\Documents\LNCA\Résultats\Electrophy4\Resutats_2021-22',saving_folder,session_name)
file_name = sprintf("AW2min1st_%s_Theta_var.mat", session_name)

load(fullfile(pathway_source, file_name))

name_struc_theta = {'LHb', 'dHPC','BLA','ACC','PRL'};
name_struc_gamma = {'LHb', 'dHPC','BLA','ACC','PRL'};

% Theta gamma coupling
MI_stat_full_rats= [];

for rat =1:length(LFP_cl_full(1,1,:))
fprintf('Rat n° %2.0f \n', list_rat{rat})
MI_stat_full= [];

for struc_theta = 1: length(name_struc_theta)
LFP_Hab_Local =LFP_cl_full(struc_theta,:,rat);

for struc_gamma = 1:length(name_struc_gamma)
LFP_DHipp_Local = LFP_cl_full(struc_gamma,:,rat);
fprintf('%s - Amp & %s - Phase', name_struc_gamma{struc_gamma}, name_struc_theta{struc_theta})

% phase amplitude coupling
% define number of bootstaps
BS = 1000;
n_nodes = 6;

Bin_scalesAmp = 0.01;
Bin_scalesPhase = 0.1;
wavename ='cmor2-0.7958';
fc = centfrq(wavename);
dt = 1/srate;
Fs = srate;

% Phase of DHipp
freqrange = [5 9]; % theta frequency 
scalerange = fc./(freqrange*(1/Fs));
scalesPhase = scalerange(end):Bin_scalesPhase:scalerange(1);
freqPhase = scal2frq(scalesPhase,wavename,dt);
b_size = 0.1;
FreqBins = freqPhase(end):b_size:1;
KeepIdx1 = [];
for n = 1:length(FreqBins)
    KeepIdx1 = [KeepIdx1,find(freqPhase>=FreqBins(n),1,'last')];
end
KeepIdx1 = unique(KeepIdx1);
FreqBins = 1+b_size:b_size:freqPhase(1);
KeepIdx2 = [];
for n = 1:length(FreqBins)
    KeepIdx2 = [KeepIdx2,find(freqPhase>=FreqBins(n),1,'last')];
end
KeepIdx2 = unique(KeepIdx2);
KeepIdx = [];
KeepIdx = sort([KeepIdx1,KeepIdx2]);
scalesPhase = scalesPhase(KeepIdx);
freqPhase = scal2frq(scalesPhase,wavename,dt);

% Amp - Hab
freqrange = [45 65]; % gamma frequency
scalerange = fc./(freqrange*(1/Fs));
scalesAmp = scalerange(end):Bin_scalesAmp:scalerange(1); % I set this to match getting ~ 325 scales from 20 to to 450, which is what we get in the LC MUA triggered spectrogram with a 9kHz signal
freqAmp = scal2frq(scalesAmp,wavename,dt);
b_size = 1;
FreqBins = freqAmp(end):b_size:1;
KeepIdx1 = [];
for n = 1:length(FreqBins)
    KeepIdx1 = [KeepIdx1,find(freqAmp>=FreqBins(n),1,'last')];
end
KeepIdx1 = unique(KeepIdx1);
FreqBins = 1+b_size:b_size:freqAmp(1);
KeepIdx2 = [];
for n = 1:length(FreqBins)
    KeepIdx2 = [KeepIdx2,find(freqAmp>=FreqBins(n),1,'last')];
end
KeepIdx2 = unique(KeepIdx2);
KeepIdx = [];
KeepIdx = sort([KeepIdx1,KeepIdx2]);
scalesAmp = scalesAmp(KeepIdx);
freqAmp = scal2frq(scalesAmp,wavename,dt);

% define phase bins
N = 18;
PhaseBins = -pi:((2*pi)/N):pi-((2*pi)/N);

AmpData = LFP_Hab_Local;
PhaseData = LFP_DHipp_Local;

% Calculate the coeeficients of the wavelets for phase (theta)
PhaseCoefs = cwt(PhaseData,scalesPhase,wavename); % freq X time

% Calculate the coefficients of the wavelet for amplitude (gamma)
AmpCoefs = cwt(AmpData,scalesAmp,wavename); % freq X time

% Get instantaneous phase, p(t)
Phase = angle(PhaseCoefs); % freq X time
clear PhaseCoefs

% Get instantaneous amplitude, A(t)
Amp = abs(AmpCoefs); % freq x time
clear AmpCoefs
clear AmpData PhaseData

mypool = parpool(n_nodes);
pause(0.5)


MI_tmp = ones(length(scalesPhase),length(scalesAmp),'single');
% Bin phases into bins, j
[~,binID] = histc(Phase,PhaseBins); % binID is freq x time
% Select frequency of interest for comodulogram
for fP = 1:size(binID,1)
    fprintf('Calculating MI - Freq for Phase %d of %d \n',fP,size(binID,1))
    selPhase = binID(fP,:);
    parfor fA = 1:size(Amp,1)
        selAmp = Amp(fA,:);
        BinnedAvgAmp = [];
        % Mean amplitude in each phase bin, <A>(j)
        for b = 0:N-1
            BinnedAvgAmp = [BinnedAvgAmp;mean(selAmp(selPhase == b))];
        end
        TotalBins = length(BinnedAvgAmp) - length(find(isnan(BinnedAvgAmp)));
        % Normalize mean power (divide each bin by the sum across all bins)
        NormAmp = BinnedAvgAmp./sum(BinnedAvgAmp); % This is P(j)
        %             bar(1:N,NormAmp)
        %             bar(min(selPhase):max(selPhase),NormAmp)
        % Calculate Shannon Entropy of Amplitude distribution, H(P)
        Entropy = -sum(NormAmp.*log(NormAmp));
        % Calculate KL Distance, D(P,U)
        KLdist = log(TotalBins) - Entropy;
        % Calculate MI
        MI_tmp(fP,fA) = KLdist/log(TotalBins);
    end
end


MI_surrogate = ones(BS,length(scalesPhase),length(scalesAmp),'single');
for bs = 1:BS
    fprintf('BS %d out of %d \n',bs,BS)
    % shuffle amplitude
    surrogate_amp = ones(length(scalesAmp),length(Amp),'single'); % initialize surrogate data
    parfor fA = 1:length(scalesAmp)
        selDat = squeeze(Amp(fA,:))';
        skip = randsample([1:length(selDat)],1);
        % shuffle phase in time (move randomly selected later part of signal to start)
        if size(selDat,1)>1
            surrogate_amp(fA,:) = [selDat(skip:end); selDat(1:skip-1)];
        else
            surrogate_amp(fA,:) = [selDat(skip:end), selDat(1:skip-1)];
        end
    end
    % Get Amplitude Distribution, P(j)
    % Bin phases into bins, j
    [~,binID] = histc(Phase,PhaseBins); % binID is freq x time
    MI_surrogate_tmp = [];
    for fP = 1:size(binID,1)
        
        selPhaseByFreq = binID(fP,:);
        parfor fA = 1:size(surrogate_amp,1)
            selAmpByFreq = surrogate_amp(fA,:);
            BinnedAvgAmp = [];
            % Mean amplitude in each phase bin, <A>(j)
            for b = 0:N-1
                binIdx = find(selPhaseByFreq == b);
                % if there are phases defined in this time
                % window
                if ~isempty(binIdx)
                    BinnedAvgAmp = [BinnedAvgAmp;mean(selAmpByFreq(binIdx))]; % mean of all amplitudes observed at that phase bin
                end
            end
            TotalBins = length(BinnedAvgAmp) - length(find(isnan(BinnedAvgAmp)));
            % Normalize mean power (divide each bin by the sum across all bins)
            NormAmp = BinnedAvgAmp./sum(BinnedAvgAmp); % This is P(j)
            % Calculate Shannon Entropy of Amplitude distribution, H(P)
            % multiply each phase bin by the log of itself. Sum across phase bins. Make negative
            Entropy = -sum(NormAmp.*log(NormAmp));
            % Calculate KL Distance, D(P,U)
            % log(N) - H(P)
            KLdist = log(TotalBins) - Entropy;
            % Calculate MI
            % normalize by log(N)
            MI_surrogate_tmp(fP,fA) = KLdist/log(TotalBins);
        end
    end
    MI_surrogate(bs,:,:) = MI_surrogate_tmp(:,:);
end

clear bs
clear Amp Phase

% Get Treshold-subtracted MI and Z-scored MI - trial averaged
MI_threshold = zeros(length(scalesPhase),length(scalesAmp),'single');
for fP = 1:length(scalesPhase)
    for fA = 1:length(scalesAmp)
        selDat = sort(MI_surrogate(:,fP,fA),'descend');
        selThresh = ceil(0.05 * BS);
        MI_threshold(fP,fA) = selDat(selThresh);
    end
end
MI_stat = MI_tmp - MI_threshold;

clear MI_tmp
clear MI_threshold MI_surrogate

delete(gcp('nocreate'))
pause(0.5)
clear('mypool')

% Change negative values to 0 since negative is not defined. It simply
% means MI was not significant, but that negative variability shouldn't
% effect the normalization.
MI_stat(MI_stat<0)=0;

% mean across rats
MI_stat = MI_stat'; % Amp X Phase
% flip x-axis LR to make frequency for phase increasing
for r = 1:size(MI_stat,1)
    selR = MI_stat(r,:);
    MI_stat(r,:) = fliplr(selR);
    clear selR
end


MI_stat_full(:,:,struc_theta, struc_gamma) = MI_stat;

end
end

MI_stat_full_rats(:,:,:,:,rat) = MI_stat_full;
end

%% Plot on the average of rats
% Do a table with the mean of the frequency data to see if at least one coupling have something significante
% Add the name of the structures on the plot
MI_stat_full_mean = mean(MI_stat_full_rats,5);

% do the plot
fig_coupling = figure('Visible','on');
fig_coupling.Position = [10 10 1800 1000];

count = 1;
average_to_save = zeros(length(name_struc_theta),length(name_struc_gamma));
for struc_theta = 1: length(name_struc_theta)
for struc_gamma = 1:length(name_struc_gamma)
MI_stat = MI_stat_full_mean(:,:, struc_theta,struc_gamma);

% do the maverage of all
average_to_save(struc_theta, struc_gamma)= mean(mean(MI_stat));


subplot(length(name_struc_theta),length(name_struc_gamma),count)
count = count + 1;
wantedFreqAmp = [45:5:65]; % put the gamma frequencies
PlotDat = MI_stat;
PlotDat = double(PlotDat); % Amp X Phase
imagesc(PlotDat);
colormap('hot')
% caxis([0 prctile(PlotDat(:),99)]) % perhaps I should put 95% 
colorbar
LabelFreqAmp = {};
for n = 1:length(freqAmp)
    sel = num2str(freqAmp(n));
    pt = find(sel == '.'); % find decimal pt
    if ~isempty(pt)
        LabelFreqAmp{n} = sel(1,1:pt-1);
    else
        LabelFreqAmp{n} = sel;
    end
end
idx = [];
for a = 1:length(wantedFreqAmp)
    idx(a) = find(freqAmp>=wantedFreqAmp(a),1,'last');
end
set(gca,'YTick',sort(idx)); set(gca,'YTickLabel',fliplr(LabelFreqAmp(idx)));
MinPhase = min(freqPhase);
MaxPhase = max(freqPhase);
LabelFreqPhase = {};
XTickIdx = [];
loopCnt = 0;
for m = floor(MinPhase):1:floor(MaxPhase)
    loopCnt = loopCnt+1;
    XTickIdx = [XTickIdx,find(freqPhase>=m,1,'last')];
    LabelFreqPhase{loopCnt} = m;
end
XTickIdx = length(freqPhase) - XTickIdx;
if XTickIdx(1) == 0
    XTickIdx(1) = 1;
end
set(gca,'XTick',XTickIdx);
set(gca,'XTickLabel',LabelFreqPhase);
set(gca,'Fontsize',12,'fontname','arial')
box off
set(gca, 'ygrid','on')
ygrid = get(gca,'YGridHandle');
set(ygrid,'Color',rgb('purple'), 'GridLineStyle','-', 'LineWidth',1.5);

xlabel(sprintf('%s Phase - freq (Hz)',name_struc_theta{struc_theta}));
ylabel(sprintf('%s Amplitude - freq (Hz)', name_struc_gamma{struc_gamma}));

end
end

% save the figure
name_file = "Coupling.fig";
saveas(fig_coupling,fullfile(pathway_source,"Figure",name_file))
name_file = "Coupling.eps";
saveas(fig_coupling,fullfile(pathway_source,"Figure",name_file))

% save the data 
name_file = ("Coupling_data.mat");
save(fullfile(pathway_source, name_file), 'MI_stat_full_rats', 'MI_stat_full_mean', ...
    'list_rat','name_struc_theta','name_struc_gamma','session_name','subsession_name',...
    'average_to_save')

% save the average of the data
T_average = array2table(average_to_save,'VariableNames',name_struc_gamma, 'RowNames', name_struc_theta);

name_file = sprintf("Average coupling.xlsx")
writetable(T_average,fullfile(pathway_source,name_file));


