function [PowdHPC,PowPRL]=power_presleep_chronus(DatabyHours,New_SR)

% params chronux
params.fpass = [1 New_SR/2]; % in Hz (0.5 to 10) - 0.5 is minimum if we want to use a 4 sec window and see one cycle
params.Fs = New_SR;
k = 9; % num tapers (lower for less frequency leakage, higher for more leakage but smoother spectrum) % and it is a tradeoff with how many time points you have. spectrum is less smooth with many time points.
nw = (k+1)/2;
params.tapers=[nw k];
moving_win=5;
win_shift=0.3;

% variable used
NbHour = length(DatabyHours(1,:));
PowdHPC = {};
PowPRL = {};

% for HPC
for i =1:NbHour
    PowdHPC{1,i}=i;
    PowPRL{1,i}=i;
        
    % power calculation
    [S,t,f]=mtspecgramc(DatabyHours{4,i}(1,:),[moving_win win_shift],params);
    
    PowdHPC{2,i} = S; % convert to dB;
    PowdHPC{3,i} = t;
    PowdHPC{4,i} = f;
end

% for PRL
for i =1:NbHour
    PowPRL{1,i}=i;

    % power calculation
    [S,t,f]=mtspecgramc(DatabyHours{4,i}(2,:),[moving_win win_shift],params);
    
    PowPRL{2,i} = S;
    PowPRL{3,i} = t;
    PowPRL{4,i} = f; 
end

sprintf('your freq bins are %2.4f Hz',unique(diff(f)))