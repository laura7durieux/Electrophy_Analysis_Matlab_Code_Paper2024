%%% Separate the waking/sleeping state
% [AW,CW,SWS,REM] = Separate_states(LFP,sleepTrack,srate) allow to separate
% the wake/sleep stage in part from your preprocess signal. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS = 
%       LFP : Should be a matrix (channels x time(in points)) or a
%           vector (1 x time(in points)) containing your LFP signal.
%       sleepTrack : is the sleep stage categorisation and should be a
%           colunm vector (categories x 1). 1 is the REM sleep, 2 is the SWS, 3
%           is calme wake, 4 is active wake.
%       srate : is your sample rate. vector (1 x 1).
% OUTPUTS =
%           AW : cell array {n x 1} with all the active wake part (one part
%           is stock in one cell.
%           CW : cell array {n x 1} with all the calm wake part (one part
%           is stock in one cell.
%           SWS : cell array {n x 1} with all the SWS part (one part
%           is stock in one cell.
%           REM : cell array {n x 1} with all the REM sleep part (one part
%           is stock in one cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [AW,CW,SWS,REM] = Separate_states(LFP,sleepTrack,srate)

%%%%%%%%%%%%%%%%Checking the inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismatrix(LFP)== 0
    error('LFP must be either a matrix or a row vector')
end

if size(LFP) <= [1 ((3599*srate)+10)]
    error('LFP must be 1 hour long !')
end

if iscolumn(sleepTrack)==0
    error('sleepTrack must be column vector')
end

if size(srate) ~= [1 1]
    error('Sampling rate must be a vector (1 x 1)---- single number')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sleepTrack = sleepTrack'; % changing orientation of sleepTrack
[rawLFP, colLFP] = size(LFP);
[rawST, colST] = size(sleepTrack);

% Create a time line for i repetition (loop below) minus 1 second because
if colLFP <= 3601*srate % if LFP is inferior at 1 H
    sleepTrackTime = (colST-1); % minus 1 second at the end of the signal
    disp('1 second at the end of the signal have been lost')
else
    sleepTrackTime = (colST);
    disp('LFP signal have been long enough for process all the signal')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process active wake parts
% Find AW    
lock = find(sleepTrack==4);

% Create a matrix with the position with the wake/sleep stage = 1
T_lock = zeros(rawST, colST);

for i = 1:length(lock)
    T_lock(:,lock(1,i)) = 1;
end

% adding a column at the end 
if T_lock(1,end) == 1
    Temp2 = [1];
    T_lock = [T_lock Temp2];
elseif T_lock(1,end) == 0
    Temp2 = [0];
    T_lock = [T_lock Temp2];
end


% Initialize the matrix needed for the process
j=1;
Temp=[];
Time_Block={};
Time = [1:colST];

% Process the cutting and put them on a cell array
for i = 1:sleepTrackTime
    if T_lock(1,i)-T_lock(1,i+1) == 0 & T_lock(1,i) == 1
        Temp=[Temp Time(:,i)];
    elseif isempty(Temp)==0 % if T is not empty
        j = j+1;
        Time_Block{j,1}=Temp;
        Temp=[];
    end
    if i == sleepTrackTime & T_lock(1,i) == 1% save before the end of the loop
        Temp=[Temp Time(:,i)];
        j = j+1;
        Time_Block{j,1}=Temp;
    end
end

Time_Block = Time_Block(2:end,:);
AW = {};

for i = 1:length(Time_Block)
    Start = srate*Time_Block{i,1}(1,1);
    Stop = srate*Time_Block{i,1}(1,end);
    AW{i,1} = LFP(:,Start:Stop);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process calme wake parts
% Find CW    
lock = find(sleepTrack==3);

% Create a matrix with the position with the wake/sleep stage = 1
T_lock = zeros(rawST, colST);

for i = 1:length(lock)
    T_lock(:,lock(1,i)) = 1;
end

% adding a column at the end 
if T_lock(1,end) == 1
    Temp2 = [1];
    T_lock = [T_lock Temp2];
elseif T_lock(1,end) == 0
    Temp2 = [0];
    T_lock = [T_lock Temp2];
end


% Initialize the matrix needed for the process
j=1;
Temp=[];
Time_Block={};
Time = [1:colST];

% Process the cutting and put them on a cell array
for i = 1:sleepTrackTime
    if T_lock(1,i)-T_lock(1,i+1) == 0 & T_lock(1,i) == 1
        Temp=[Temp Time(:,i)];
    elseif isempty(Temp)==0 % if T is not empty
        j = j+1;
        Time_Block{j,1}=Temp;
        Temp=[];
    end
    if i == sleepTrackTime & T_lock(1,i) == 1% save before the end of the loop
        Temp=[Temp Time(:,i)];
        j = j+1;
        Time_Block{j,1}=Temp;
    end
end

Time_Block = Time_Block(2:end,:);
CW = {};

for i = 1:length(Time_Block)
    Start = srate*Time_Block{i,1}(1,1);
    Stop = srate*Time_Block{i,1}(1,end);
    CW{i,1} = LFP(:,Start:Stop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process slow wave sleep parts
% Find SWS    
lock = find(sleepTrack==2);

% Create a matrix with the position with the wake/sleep stage = 1
T_lock = zeros(rawST, colST);

for i = 1:length(lock)
    T_lock(:,lock(1,i)) = 1;
end

% adding a column at the end 
if T_lock(1,end) == 1
    Temp2 = [1];
    T_lock = [T_lock Temp2];
elseif T_lock(1,end) == 0
    Temp2 = [0];
    T_lock = [T_lock Temp2];
end


% Initialize the matrix needed for the process
j=1;
Temp=[];
Time_Block={};
Time = [1:colST];

% Process the cutting and put them on a cell array
for i = 1:sleepTrackTime
    if T_lock(1,i)-T_lock(1,i+1) == 0 & T_lock(1,i) == 1
        Temp=[Temp Time(:,i)];
    elseif isempty(Temp)==0 % if T is not empty
        j = j+1;
        Time_Block{j,1}=Temp;
        Temp=[];
    end
    if i == sleepTrackTime & T_lock(1,i) == 1% save before the end of the loop
        Temp=[Temp Time(:,i)];
        j = j+1;
        Time_Block{j,1}=Temp;
    end
end

Time_Block = Time_Block(2:end,:);
SWS = {};

for i = 1:length(Time_Block)
    Start = srate*Time_Block{i,1}(1,1);
    Stop = srate*Time_Block{i,1}(1,end);
    SWS{i,1} = LFP(:,Start:Stop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process REM sleep parts
% Find REM    
lock = find(sleepTrack==1);

% Create a matrix with the position with the wake/sleep stage = 1
T_lock = zeros(rawST, colST);

for i = 1:length(lock)
    T_lock(:,lock(1,i)) = 1;
end

% adding a column at the end 
if T_lock(1,end) == 1
    Temp2 = [1];
    T_lock = [T_lock Temp2];
elseif T_lock(1,end) == 0
    Temp2 = [0];
    T_lock = [T_lock Temp2];
end

% Initialize the matrix needed for the process
j=1;
Temp=[];
Time_Block={};
Time = [1:colST];

% Process the cutting and put them on a cell array
for i = 1:sleepTrackTime
    if T_lock(1,i)-T_lock(1,i+1) == 0 & T_lock(1,i) == 1
        Temp=[Temp Time(:,i)];
    elseif isempty(Temp)==0 % if T is not empty
        j = j+1;
        Time_Block{j,1}=Temp;
        Temp=[];
    end
    if i == sleepTrackTime & T_lock(1,i) == 1% save before the end of the loop
        Temp=[Temp Time(:,i)];
        j = j+1;
        Time_Block{j,1}=Temp;
    end
end

Time_Block = Time_Block(2:end,:);
REM = {};

for i = 1:length(Time_Block)
    Start = srate*Time_Block{i,1}(1,1);
    Stop = srate*Time_Block{i,1}(1,end);
    REM{i,1} = LFP(:,Start:Stop);
end

%% End of the process
end

