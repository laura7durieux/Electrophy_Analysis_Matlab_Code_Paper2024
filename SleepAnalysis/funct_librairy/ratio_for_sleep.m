


function [theta_delta_ratio,delta_gamma_ratio] = ratio_for_sleep(dHPC_S,dHPC_f,PRL_S,PRL_f)

Hand = struct();

% Caculate the power range for HPC
    delta_loc=find(dHPC_f>0.5 & dHPC_f<3); % identify all values of the spectrogram that are within the delta range
    theta_loc=find(dHPC_f>5 & dHPC_f<12); % identify all values of the spectrogram that are within the theta range
 
    Hand.theta_pow_sum=sum(dHPC_S(:,theta_loc),2);% take the sum under the curve for theta accross time
    Hand.theta_pow_mean=nanmean(Hand.theta_pow_sum);
   
    Hand.delta_pow_sum=sum(dHPC_S(:,delta_loc),2);% take the sum under the curve for delta accross time
    Hand.delta_pow_mean=nanmean(Hand.delta_pow_sum);
    
    Hand.theta_delta_ratio=Hand.theta_pow_sum./Hand.delta_pow_sum; % calcul theta/ delta ratio
    Hand.theta_delta_ratio_mean=nanmean(Hand.theta_delta_ratio);
    

% Caculate the power range for PFC
    delta_loc=find(PRL_f>0.5 & PRL_f<3); % identify all values of the spectrogram that are within the delta range
    gamma_loc=find(PRL_f>30 & PRL_f<120); % identify all values of the spectrogram that are within the gamma range

    Hand.gamma_pow_sum=sum(PRL_S(:,gamma_loc),2);% take the sum under the curve for delta accross time
    Hand.gamma_pow_mean=nanmean(Hand.gamma_pow_sum);
   
    Hand.delta_pow_sum=sum(PRL_S(:,delta_loc),2);% take the sum under the curve for delta accross time
    Hand.delta_pow_mean=nanmean(Hand.delta_pow_sum);
    
    Hand.delta_gamma_ratio=Hand.delta_pow_sum./Hand.gamma_pow_sum; % calcul theta/ delta ratio
    Hand.delta_gamma_ratio_mean=nanmean(Hand.delta_gamma_ratio);
    

% Outputs
theta_delta_ratio = Hand.theta_delta_ratio;
delta_gamma_ratio = Hand.delta_gamma_ratio;

end