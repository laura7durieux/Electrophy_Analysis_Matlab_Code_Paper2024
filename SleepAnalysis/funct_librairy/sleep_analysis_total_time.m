

function [sleepsettings_total,sleepsettings_total_table] = sleep_analysis_total_time(managingPath)

%% load data
for i=2:length(managingPath)
    pathData = managingPath{i,1};
    load(pathData)
    SleepSettings{1,i-1} = Results.SleepSettings; 
end


%% Latency calculation (raw 2)
sleepsettings_total(1,:) = SleepSettings{1,1}(1,:);
TempLat = zeros(length(SleepSettings(1,:)),length(SleepSettings{1,1}(1,:)));
TempLat2 = zeros(length(SleepSettings(1,:)),length(SleepSettings{1,1}(1,:)));

for i = 1:length(SleepSettings(1,:)) % period
    TempLat(i,:) =  cell2mat(SleepSettings{1,i}(2,:));
end

    
for j = 1:length(TempLat(1,:)) % stage
    for i = 1:length(SleepSettings(1,:)) % period

        if isnan(TempLat(i,j)) == 0 % not a NaN
           TempLat2(i,j) = TempLat(i,j);
           if i > 1  
                TempLat2(i,j) = TempLat(i,j)+ TempLat2(i-1,j);
            end
                   
        else % is NaN
            TempLat2(i,j) = 3600; %1 hour
            if i > 1  
                TempLat2(i,j) = TempLat2(i,j)+ TempLat2(i-1,j);
            end
        end

    end    
end


for j = 1:length(TempLat(1,:)) % stage
    v = TempLat2(:,j)/3600;
    v100 = v(mod(v,1)==0);
    if isempty(v100) ==1 
        sleepsettings_total{2,j} = TempLat2(1,j); 
    else
        sleepsettings_total{2,j} = TempLat2(max(v100)+1,j);
    end
end

%% number of occurancies (raw 3)
TempOcc = zeros(length(SleepSettings(1,:)),length(SleepSettings{1,1}(1,:)));

for i = 1:length(SleepSettings(1,:)) % period
    TempOcc(i,:) =  cell2mat(SleepSettings{1,i}(3,:));
end

for j = 1:length(TempOcc(1,:)) % stage
    sleepsettings_total{3,j}=nansum(TempOcc(:,j));
end

%% Sum of time time passed in a stage second (raw4)

TimeTot = zeros(length(SleepSettings(1,:)),length(SleepSettings{1,1}(1,:)));

for i = 1:length(SleepSettings(1,:)) % period
    TimeTot(i,:) =  cell2mat(SleepSettings{1,i}(5,:));
end

for j = 1:length(TimeTot(1,:)) % stage
    sleepsettings_total{4,j}=nansum(TimeTot(:,j));
end

%% Average of time of one episod spend in a stage (raw 5)
TimeAve = zeros(length(SleepSettings(1,:)),length(SleepSettings{1,1}(1,:)));

for i = 1:length(SleepSettings(1,:)) % period
    TimeAve(i,:) =  cell2mat(SleepSettings{1,i}(5,:));
end

for j = 1:length(TimeAve(1,:)) % stage
    sleepsettings_total{5,j}=nanmean(TimeAve(:,j));
end

%% Create table 

varNames = {'Stages','StartLantencyS','OccurenciesNb','TimeSpendS','AverTimeSpendS'};
sleepsettings_total_table = cell2table(sleepsettings_total','VariableNames',varNames);





