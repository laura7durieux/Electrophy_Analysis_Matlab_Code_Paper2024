%% Function allowing to extract informations from sleep analysis

function [Results_stage,Results_stage_table] = SleepAnalysis(Stage,sleepTrack,srate)

%Stage= { AW; CW; SWS; REM};
% Calculate how many parts for each stage

% put away the sleep/wake stage if one of them are empty
sleepTrack = sleepTrack';

nameStageT = {'AW', 'CW', 'SWS', 'REM'};
Results_stage = {};


%%%%%% find the beginning of the sleep/wake stage %%%%%%%%%%%%%%%%%%%%%%
for i = 1: length(nameStageT)
    % put the name at the start of the cell array
    Results_stage{1,i}=nameStageT{i};
    
    % take the start time already calculated before 
    if isempty(Stage{i})==0 % not empty
        Results_stage{2,i}=Stage{i,1}{1,2}; % en second
    else
        Results_stage{2,i}=NaN;
    end
end
  

%%%%%%%% number of time block (occurences) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1: length(nameStageT)
    
    if isempty(Stage{i})==0 % not empty
        Results_stage{3,i}=length(Stage{i,1}(:,1)); % number of periode of the same stage 
    else
        Results_stage{3,i}=NaN;
    end
end

%%%%%%%%%%%% Time of each occurences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1: length(nameStageT) % stages
    TimePeriodTemp =[];
    if isempty(Stage{i})==0 % not empty
        for j = 1: length(Stage{i,1}(:,1)) % periods
            TimePeriodTemp = [TimePeriodTemp length(Stage{i,1}{j,1}(1,:))];
        end
         Results_stage{4,i}=TimePeriodTemp/srate; % en second
    else
         Results_stage{4,i}=NaN;
    end
end

%%%%%%%%%% Time total of occurencies %%%%%%%%%%%%%%%%%%%%%%%

for i = 1: length(nameStageT) % stages
    if isempty(Stage{i})==0 % not empty
         Results_stage{5,i}=sum(Results_stage{4,i}(1,:)); % en second
    else
         Results_stage{5,i}=NaN;
    end
end


%%%%%%%%%%% mean Spend on the stage considered %%%%%%%%%%%%%%%%%%%

for i = 1: length(nameStageT) % stages
    if isempty(Stage{i})==0 % not empty
         Results_stage{6,i}=mean(Results_stage{4,i}(1,:)); % en second
    else
         Results_stage{6,i}=NaN;
    end
end


%%%%%%%%%%%%%%%%%%% Table %%%%%%%%%%%%%%%%%%%%%%%%%

Results_stage_table = [Results_stage(1:3,:) ; Results_stage(5:6,:)]';


varNames = {'Stages','StartLantencyS','OccurenciesNb','TimeSpendS','AverTimeSpendS'};
Results_stage_table = cell2table(Results_stage_table,'VariableNames',varNames);

disp('The analysis on sleep stage is done')

end