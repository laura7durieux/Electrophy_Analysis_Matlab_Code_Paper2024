function [DatabyHours] = select_ch_for_sleep(DatabyHours,localStruc,localParams,dHPCExample,PRLExample)

% Select the position of the ch needed
LockHPC = find(localStruc == 'dHPC');
LockPRL = find(localStruc == 'PRL');
LockEMG = find(localStruc == 'EMG');


for i = 1:length(DatabyHours(1,:))
    
    % for dHPC (first raw in the matrice into the 4th raw of the cells)
    if localParams(LockHPC,1)== 1
        DatabyHours{4,i}(1,:)=DatabyHours{3,i}(LockHPC,:);
    
    elseif localParams(LockHPC,1)== 0.5
        % plot for chosing what is better
        figure
        subplot(2,1,1)
        plot(DatabyHours{3,i}(LockHPC,:))
        title('Figure 1 : After local ref')
        ylabel('en mV')
        xlabel('en second')
        subplot(2,1,2)
        plot(DatabyHours{2,i}(dHPCExample,:))
        x =['Figure 2 : the channel n°',num2str(dHPCExample),' non-local'];
        title(x)
        ylabel('en mV')
        xlabel('en second')
        
        % choose the better
        Input1 = input('Which figure is better ? (1 or 2)');
        if Input1 == 1
            DatabyHours{4,i}(1,:) = DatabyHours{3,i}(LockHPC,:);
        elseif Input1 == 2
            DatabyHours{4,i}(1,:) = DatabyHours{3,i}(dHPCExample,:);
        else 
            error('You have to choose one of the two possibilities !')
        end
        
    else 
        DatabyHours{4,i}(1,:) = DatabyHours{2,i}(dHPCExample,:);
        x=['Warning : the channel ',num2str(dHPCExample),' had been choose cause no Local ref available'];
        disp(x)   
    end
    
    % for PRL (Second raw in the matrice into the 4th raw of the cells)
    if localParams(LockPRL,1)== 1
        DatabyHours{4,i}(2,:)=DatabyHours{3,i}(LockPRL,:);
    
    elseif localParams(LockPRL,1)== 0.5
        % plot for chosing what is better
        figure
        subplot(2,1,1)
        plot(DatabyHours{3,i}(LockPRL,:))
        title('Figure 1 : After local ref')
        ylabel('en mV')
        xlabel('en second')
        subplot(2,1,2)
        plot(DatabyHours{2,i}(PRLExample,:))
        x =['Figure 2 : the channel n°',num2str(PRLExample),' non-local'];
        title(x)
        ylabel('en mV')
        xlabel('en second')
        
        % choose the better
        Input1 = input('Which figure is better ? (1 or 2)');
        if Input1 == 1
            DatabyHours{4,i}(2,:) = DatabyHours{3,i}(LockPRL,:);
        elseif Input1 == 2
            DatabyHours{4,i}(2,:) = DatabyHours{2,i}(PRLExample,:);
        else 
            error('You have to choose one of the two possibilities !')
        end
        
    else 
        DatabyHours{4,i}(2,:) = DatabyHours{2,i}(PRLExample,:);
        x=['Warning : the channel ',num2str(PRLExample),' had been choose cause no Local ref available'];
        disp(x)   
    end
    
    % for EMG (5th raw of the cells)
    DatabyHours{5,i}(1,:) = DatabyHours{3,i}(LockEMG,:);
    
    % checking if EMG is available
    CheckEMG = localParams(LockEMG);
    if CheckEMG==0
        warning(['The EMG is not avalaible for the hour n°',num2str(i),', you should load DeepLabCut movment analysis !'])
    end
end