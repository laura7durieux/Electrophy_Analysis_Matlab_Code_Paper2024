
function [AW,CW,SWS,REM] = Separate_statesV2(LFP,SleepTrack,srate)

% changing the orientation of sleeptrack 
ST = SleepTrack';

% initiation of the matrix
AW = {};
CW = {};
SWS ={};
REM = {};


% for running that we need to add an add on : MinGW-w64 Compiler into
% matlab
% We need also all the toolbox https://www.mathworks.com/matlabcentral/fileexchange/41813-runlength
X = ST;
[B, N, Ind] = RunLength(X);
Ind         = [Ind, length(X)+1];
Multiple    = find(N > 1);
Start       = Ind(Multiple);
Stop        = Ind(Multiple + 1) - 1;

aw =1; cw=1; sws=1; rem=1;

% Create the results cell
for i = 1:length(Start)
    if ST(Start(i))==4
        AW{aw,1} = LFP(:,Start(i)*srate:Stop(i)*srate);
        AW{aw,2} = Start(i);
        AW{aw,3} = Stop(i);
        aw = aw+ 1;
    elseif ST(Start(i))==3
        CW{cw,1} = LFP(:,Start(i)*srate:Stop(i)*srate);
        CW{cw,2} = Start(i);
        CW{cw,3} = Stop(i);
        cw =cw+1;
    elseif ST(Start(i))==2
        SWS{sws,1} = LFP(:,Start(i)*srate:Stop(i)*srate);
        SWS{sws,2} = Start(i);
        SWS{sws,3} = Stop(i);
        sws =sws+1;
    elseif ST(Start(i))==1
        REM{rem,1} = LFP(:,Start(i)*srate:Stop(i)*srate);
        REM{rem,2} = Start(i);
        REM{rem,3} = Stop(i);
        rem=rem+1;
    else
        disp("Something went wrong");
    end
        
end
