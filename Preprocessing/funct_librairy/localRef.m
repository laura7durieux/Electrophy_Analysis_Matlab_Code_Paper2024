%% SCRIPT CALCULATING THE LOCAL REF DEPENDING OF THE EXCEL FILE CHECKING ELECTRODES VIABILITY 
% if writte that : Warning :No LFP local ref for this structure :      EMG 
% it's not true !


function [LFPLocal,localParams, StructuresList]=localRef(ChData,Structures,LFP)

SelectData = {};
localParams = [];

% Find the structures names 
StructuresList = unique(Structures,'stable');

for i = 1:length(StructuresList(:,1))
    SelectData{1,i} = StructuresList(i,1);
    SelectData{2,i} = find(Structures == StructuresList(i,1));
end

for i = 1:length(SelectData(1,:))
    for j=1:length(SelectData{2,i})
        SelectData{3,i}(j,1)= ChData(SelectData{2,i}(j,1),1);
        
        if SelectData{3,i}(j,1)==1
            SelectData{4,i}(j,:)= LFP(SelectData{2,i}(j,1),:);
            SelectData{5,i}(j,:)= zeros(1,length(LFP(SelectData{2,i}(j,1),:)));
        elseif SelectData{3,i}(j,1)==0.5
            SelectData{5,i}(j,:)= LFP(SelectData{2,i}(j,1),:);
            SelectData{4,i}(j,:)= zeros(1,length(LFP(SelectData{2,i}(j,1),:)));
            X = ['Channel n°',num2str(SelectData{2,i}(j,1)),' is not stable of ',...
                SelectData{1,i},' analysis.'];
            disp(X)
        elseif SelectData{3,i}(j,1)== 0
            SelectData{4,i}(j,:)=zeros(1,length(LFP(SelectData{2,i}(j,1),:)));
            SelectData{5,i}(j,:)= zeros(1,length(LFP(SelectData{2,i}(j,1),:)));
            X = ['Channel n°',num2str(SelectData{2,i}(j,1)),' is out of ',...
                SelectData{1,i},' analysis.'];
            disp(X)
        end  
    end 
end

% Sorting the data for putting them one after another 
for i = 1:length(SelectData(1,:))
    k=1;
    for j=1:length(SelectData{2,i})
        if sum(SelectData{4,i}(j,:))~= 0
            SelectData{6,i}(k,:)= SelectData{4,i}(j,:);
            k = k+1;
        end
    end
    
    k=1;
    for j=1:length(SelectData{2,i})
        if sum(SelectData{5,i}(j,:))~= 0
            SelectData{7,i}(k,:)= SelectData{5,i}(j,:);
            k=k+1;
        end
    end
end


% Calculating the local referencing 
for i = 1:length(SelectData(1,:))
    if isempty(SelectData{6,i})==0 % Not empty
        if length(SelectData{6,i}(:,1))>= 2
            SelectData{8,i}= SelectData{6,i}(1,:)-SelectData{6,i}(2,:);
            X = ['Both channels are okay for this structure :',...
                SelectData{1,i}];
            disp(X)
            localParams =[localParams; 1];
        end
    end
    
    if isempty(SelectData{6,i})==0 & isempty(SelectData{7,i})==0 % Not empty
        if (length(SelectData{6,i}(:,1)) == 1) & (length(SelectData{7,i}(:,1))>= 1)
            SelectData{8,i}= SelectData{5,i}(1,:)-SelectData{6,i}(1,:);
            X = ['Warning ! One of the channel for the local ref is not stable for this structure :',...
                SelectData{1,i}];
            disp(X)
            localParams =[localParams; 0.5];
        end
    end
    
    if isempty(SelectData{7,i})==0 % Not empty
        if length(SelectData{7,i}(:,1))>= 2
            SelectData{8,i}= SelectData{7,i}(1,:)-SelectData{7,i}(2,:);
            X = ['Warning : Both channels are unstable for this structure :',...
                SelectData{1,i}];
            disp(X)
            localParams =[localParams; 0.5];
        end
    end
    if isempty(SelectData{8,i})==1 % Empty
        [~,SizeLFP]= size(LFP(1,:));
       SelectData{8,i} = NaN(1,SizeLFP);
       X = ['Warning :No LFP local ref for this structure :',...
                SelectData{1,i}];
            disp(X)
            localParams =[localParams; 0];
    end
end


% for EMG (because needs to be out of the local ref calculation)
if sum(SelectData{4,2}) ~= 0
    SelectData{8,2} = SelectData{4,2}; % it's the second col because it's where EMG is stored
else
    SelectData{8,2} = NaN(1,length(SelectData{4,2}));
end

% for output
LFPLocal=[];
for i = 1:length(SelectData(1,:))
    if isempty(SelectData{8,i})==0 % Not empty
        LFPLocal = [LFPLocal ; SelectData{8,i}(1,:)];
    end
end



