%% Breath First Search
clear;
I = [0,0];
G = [0,4];

F={};
Fparent={};
i=1;
P=I;
S(i,:,:) = P(:,:);
Sparent(i) = 0;

while not(isequal(P,G))
    x=P(1);
    y=P(2);
    
    %fill 3L pitcher
    if x==0
        C=action(P,x,y,'fill_3');
        F{end+1} = C;
        Fparent{end+1} = i;
    end

    %Empty 3L pitcher
    if x>0
        C=action(P,x,y,'empty_3');
        F{end+1} = C;
        Fparent{end+1} = i;
    end
    
    %Empty 5L pitcher
    if y>0
        C=action(P,x,y,'empty_5');
        F{end+1} = C;
        Fparent{end+1} = i;
    end
    
    %Empty 3L into 5L pitcher
    if x>0
        C=action(P,x,y,'empty_3to5');
        F{end+1} = C;
        Fparent{end+1} = i;
    end

    %Empty 5L into 3L pitcher
    if y>0
        C=action(P,x,y,'empty_5to3');
        F{end+1} = C;
        Fparent{end+1} = i;
    end

    % visit a new child state in the frontier
    P = F{1}; % breadth-first algorithm
    i = i + 1;
    S(i,:,:) = P(:,:); % memorize the new state 
    Sparent(i) = Fparent{1};
    
    % remove the visited state (first)
    F = F(2:end);
    Fparent = Fparent(2:end);
    
end

% intialize
P(:,:) = S(i,:,:); % restart from the final state (goal)
idp = i;
path = {P};

while not(isequal(P,I))
    idp = Sparent(idp);
    P(:,:) = S(idp,:,:);
    path{end+1} = P;
end
path = flip(path);

disp("PATH")
for j=1:length(path)
    disp(path{j})
end

disp("Breath Frist Algo")
disp("Number of nodes explored to reach the goal")
disp(i)

%% Depth First Search
clear;
I = [0,0];
G = [0,4];
F={[0,0]};
Fparent={};
i=1;
P=I;
S(i,:) = P(:);
steps = 0;

while not(isequal(P,G))
    
    x=P(1);
    y=P(2);
    n = length(F); 
    %fill 3L pitcher
    if x==0
        C=action(P,x,y,'fill_3');
        F = {C, F{:}};
    end
    %Empty 3L pitcher
    if x>0
        C=action(P,x,y,'empty_3');
        F = {C, F{:}};
    end
    %Empty 5L pitcher
    if y>0
        C=action(P,x,y,'empty_5');
        F = {C, F{:}};
    end
    %Empty 3L into 5L pitcher
    if x>0
        C=action(P,x,y,'empty_3to5');
        F = {C, F{:}};
    end
    %Empty 5L into 3L pitcher
    if y>0
        C=action(P,x,y,'empty_5to3');
        F = {C, F{:}};
    end
    
    % visit a new child state in the frontier
    P_check = F{1}; 
    
    if ismember(P_check,S,'rows')
        %already expansed node -- so choose a different node that is not
        %expansed from the frontier
        check=1;
        while ismember(P_check,S,'rows')
            check = check+1;
            P_check= F{check};
            steps = steps+1;
        end
        P = P_check;
        i = i + 1;
        S(i,:) = P(:);
    else 
        steps = steps+1;
        P = F{1}; 
        i = i + 1;
        S(i,:) = P(:);
    end    

    %remove the next visited state
    F(end+1-n) = [];

end

disp("Depth Frist Algo")
disp("PATH")
disp(S)
disp("Number of nodes explored to reach the goal")
disp(steps)

%% Greedy Algorithm

clear;
xlsname = '/Volumes/GoogleDrive/My Drive/Academics/2020-2021/Fall/ML_AI/Week2/dataHW1_Fall2020.xlsx';
data = xlsread(xlsname,'Sheet1','B4:C103');

I = data(1,:);
G = data(end,:);
P={I};
plot(data(:,1),data(:,2),'o');
hold on;

while not(isequal(P{end},G))
    F={};
    data_index=[];
    lastP = P{end};
    %Add points to frontier that are 1.3 units away
    for i = 1:length(data)
        d = distance(data(i,:),lastP);
        if d<=1.3
            F{end+1}=data(i,:);
            data_index(end+1)=i;
        end
    end
    
    %find the point that is least away from goal
    dis_min=sqrt(200);
    min_pos=[];
    remove_data_index=[];
    for i = 1:length(data_index)
        d = distance(data(data_index(i),:),G);
        if d<dis_min
            min_pos=i;
            dis_min=d;
            remove_data_index=data_index(i);
        end
    end
    P{end+1}=F{min_pos};
    data(remove_data_index,:)=[];
end

l = length(P);
path = zeros(l,2);

for i=1:l
    path(i,:)= P{i};
end

plot(path(:,1),path(:,2),'r', 'linewidth',2)

%% Functions 

function d_ = distance(x,y)
    d_ = sqrt(sum((x-y).^2));
end

function C_ = action(P,x,y,s)
    
    if strcmp(s,'empty_3')==1
        C_=P;
        C_(1)=0;
    end
    if strcmp(s,'empty_5')==1
        C_=P;
        C_(2)=0;
    end
    if strcmp(s,'empty_3to5')==1
        C_=P;
        if x+y>5
            C_(1)= x+y-5;
            C_(2)=5;
        else
            C_(1)=0;
            C_(2)=x+y;
        end
    end
    if strcmp(s,'empty_5to3')==1
        C_=P;
        if x+y>3
            C_(2)= x+y-3;
            C_(1)=3;
        else
            C_(2)=0;
            C_(1)=x+y;
        end
    end
    if strcmp(s,'fill_3')==1
        C_=P;
        C_(1)=3;
    end
end
%%
