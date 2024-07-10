% % Bootstrap on cohence analysis to evaluate the theshold where the coherence can be considered 
% Loading the data
close all
clear
 
% set up variable for loading 
session_name = 'Stress2';
subsession_name = ['PS2'];
saving_folder = "Stress2"

pathway_saving = fullfile('E:\Documents\LNCA\Résultats\Electrophy4\Resutats_2021-22',saving_folder);
pathway_source = fullfile('E:\Documents\LNCA\Résultats\Electrophy4\Resutats_2021-22',saving_folder,"Figures")
file_name = sprintf("Stress_%s_var.mat", session_name)

load(fullfile(pathway_source, file_name))

% clear what is not needed
clear C_smooth f_smooth ftot Stot_Z T_coh_full_gamma T_coh_full_theta T_coh_gamma T_coh_theta T_peak_gamma T_peak_theta T_power_gamma T_power_theta tc ...
    ttot

% Parameters - Coherance avec Chronux. 
band = [1 80];
params.Fs = srate;
params.fpass = band;
params.err = [1 0.05]; % Set the p value we want in the error calculation by chronux
k = 9; % num tapers (lower for less frequency leakage, higher for more leakage but smoother spectrum) % and it is a tradeoff with how many time points you have. spectrum is less smooth with many time points.
nw = (k+1)/2;
params.tapers=[nw k];

[Coh_dat,~,~,~,~,Cohf,~,~] = coherencyc(LFP_cl_full(1,:,1),...
    LFP_cl_full(2,:,1),params);
Coh_smooth = fastsmooth(Coh_dat,1000,3,1);


desiredP = 0.05;
BS = 1000;
pval = desiredP/BS; % This value is bonferoni correction but it is a very conversative one.

Coh_dat_surr = ones(BS,length(Coh_dat(:,1)));
Coh_surr_smooth = ones(BS,length(Coh_dat(:,1)),length(matrixComp1));

T_theta = table();
T_gamma = table();

% % Bootstraping & Smoothing

% This code with 100000 bootstrap takes 24h to run, don't lauch it, if
% you're not sur to have the time - 10 000 takes 15 min
tic
% For theta band - REPLACED BY DELTA FOR SWS
lpf_theta = 9;
hpf_theta = 5;

% For gamma band
lpf_gamma = 65;
hpf_gamma = 45;

% Le gamma et les loops restent à faire
results_theta = struct();
results_gamma = struct();
data_frq =[];

for rat=1:length(list_rat)
fprintf('Rat n° %d (%1.0f / %2.0f)\n',list_rat(rat), rat, length(list_rat))

for s = 1: length(matrixComp1)

indice_results = s + ((rat-1)*length(matrixComp1));

% real coherence
LFP_1 = LFP_cl_full(matrixComp1(s),:,rat);
LFP_2 = LFP_cl_full(matrixComp2(s),:,rat);

[Coh_dat,~,~,~,~,Cohf,~,~] = coherencyc(LFP_1,...
    LFP_2,params);

Coh_smooth(:,s) = fastsmooth(Coh_dat,1000,3,1);



%% Bootstraping
% Creates the surrogate and calculate the coherence on it 
fprintf('%s & %s \n',coherency_list1{s}, coherency_list2{s})
for B = 1:BS
    if mod(B,1000) == 0
    fprintf('Calculating Coherence for BS %d \n',B)
    end
    surrogate_dat = []; % initialize surrogate data
    skip = randsample([1:length(LFP_1)],1);
    surrogate_dat = [LFP_1(skip:end), LFP_1(1:skip-1)]; % shuffle phase in time (move randomly selected later part of signal to start)
    [Coh_dat_surr(B,:),~,~,~,~,~,~,~] = coherencyc(surrogate_dat,LFP_2,params);
    clear surrogate_dat
end

% Smoothing of the surogate coherency
for B= 1:size(Coh_dat_surr,1)
    Coh_surr_smooth(B,:,s) = fastsmooth(Coh_dat_surr(B,:),1000,3,1);
end


% calulate the 95 % Confidence BS over frequencies
% For theta band
Fs_freq = length(Cohf(1,:))/max(Cohf(1,:));
stop_theta = round(lpf_theta*Fs_freq);      %in points (positions)
start_theta = round(hpf_theta*Fs_freq);     %in points (positions)

data_theta = mean(Coh_smooth(start_theta:stop_theta,s));
surrogate_data_theta = mean(Coh_surr_smooth(:,start_theta:stop_theta,s),2); 

[theta] =  calculate_stat_BS(data_theta, surrogate_data_theta, pval);

results_theta.rat(indice_results) = list_rat(rat);
results_theta.structure{indice_results} = sprintf('%s & %s', coherency_list1{s}, coherency_list2{s});

field = fieldnames(theta);
for i = 1: length(field) % find a way to go around that issue
    results_theta.(field{i})(indice_results) = theta.(field{i});
end

% For gamma band
stop_gamma = round(lpf_gamma*Fs_freq);      %in points (positions)
start_gamma = round(hpf_gamma*Fs_freq);     %in points (positions)

data_gamma = mean(Coh_smooth(start_gamma:stop_gamma,s));
surrogate_data_gamma = mean(Coh_surr_smooth(:,start_gamma:stop_gamma,s),2); 

[gamma] =  calculate_stat_BS(data_gamma, surrogate_data_gamma, pval);

results_gamma.rat(indice_results) = list_rat(rat);
results_gamma.structure{indice_results} = sprintf('%s & %s', coherency_list1{s}, coherency_list2{s});

for i = 1: length(field) % find a way to go around that issue
    results_gamma.(field{i})(indice_results) = gamma.(field{i});
end

% calculate the p over frequency
for frq= 1:length(Coh_smooth(:,s))
    data_temps = Coh_smooth(frq,s);
    surrogate_data_temps = Coh_surr_smooth(:,frq,s); 
    temp =  calculate_stat_BS(data_temps, surrogate_data_temps, pval);
    data_frq(frq,s,rat) = temp.pval_eval;
end

end
end

% reorganize the fields inside structures to cerate the table
field_all = fieldnames(results_theta);
for f = 1: length(field_all)
    results_theta.(field_all{f}) = results_theta.(field_all{f})';
    results_gamma.(field_all{f}) = results_gamma.(field_all{f})';
end
toc

% Storing the data into a table
% creating the tables
T_theta = struct2table(results_theta,'ASArray',false);
T_gamma = struct2table(results_gamma,'ASArray',false);

% Saving the data as excel files
name_file = sprintf("Theta_Coherence_bootstraped.xlsx")
writetable(T_theta,fullfile(pathway_source,name_file),'Sheet',sprintf("%s", session_name));

name_file = sprintf("Gamma_Coherence_bootstraped.xlsx")
writetable(T_gamma,fullfile(pathway_source,name_file),'Sheet',sprintf("%s", session_name));

% After nelson totah advice :
% Calculate the p of each frequency DONE
% Plot the % of significante rats across frequency (by duo of structures) DONE
% Do the calculation of the number of rat significante (in %)
% format of the data : data_frq(frq,s,rat)
NBpercent_rat = [];
number_s = length((data_frq(1,:,1)));
n= [];
edges = []; 

for s = 1: number_s
for freq = 1: length(Cohf)
unNum   = unique(data_frq(freq,s,:));
[n(freq,:),edges(freq,:)] = histcounts(data_frq(freq,s,:),2);
end
% put the number in percentage on the total number of rat
NBpercent_rat(:,s)= (n(:,2)*100)/length(data_frq(1,1,:)); % when edges = 1 (either 1 or 0)
end

% Plot
fig_coh = figure();
fig_coh.Position = [10 10 1800 1000]; 
number_s = length(matrixComp1);
for s = 1: number_s
    subplot(2,number_s/2,s)
    bar(Cohf,NBpercent_rat(:,s),'k')
    ylabel('Significante coherence (% of rat)')
    xlabel('Frequency (Hz)')
    title(sprintf('%s & %s', coherency_list1{s}, coherency_list2{s}))
end
% Save the plot
name_file = "coherence_over_rat.fig"
saveas(fig_coh,fullfile(pathway_source,name_file))

name_file = "coherence_over_rat.svg"
saveas(fig_coh,fullfile(pathway_source, name_file))

% save important variables there
save(fullfile(pathway_source,'Variable_coherence_bs.mat'),'Coh_smooth','Coh_surr_smooth',...
    'Cohf','matrixComp1','matrixComp2','NBpercent_rat','results_gamma','results_theta',...
    'T_gamma','T_theta')


% This part is where I stored the function needed for this script
% [results] =  calculate_stat_BS(data, surrogate_data, pval)
% This function calculate the 95 interval confidance of the bootstrap as well as the p value of the original data compared to the bootstrap distribution. Be carrefull it was wrotte only for one frequency value
% INPUT : 
% data = Original coherance data to test (1 x 1)
% surrogate_data = the surrogate results of the coherance (number of BS x 1)
% pval = is the p value wanted (eg 0.05) corrected depending of the number of BS repetition (Bonferoni)
% OUTPUT : 
% results : structures containing CI 95 min and max, as well as the comparison between the real value and the C95. More over it give the estimate p value and the comparison with the p value corrected (cad the significance of our result)
function [results] =  calculate_stat_BS(data, surrogate_data, pval)
    
    % calculate the CI 95 of the bootstrap
    results.CI95min = prctile(surrogate_data,2.5);
    results.CI95max = prctile(surrogate_data,97.5);
    results.CI95 = data < results.CI95min | data > results.CI95max;

    % calculate the p value and the significance of the BS
    [n,xout] = hist(surrogate_data,100);
    col = find(xout >= data);
    BS = length(surrogate_data(:,1));
    if isempty(col)
        results.p = 0; 
        results.pval_eval = 1;
    else 
        prob = sum(n(col));
        results.p = prob/BS; 
        results.pval_eval = results.p < pval;
    end
    results.pval = pval;
end



