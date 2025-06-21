function [sol]=EAassignment

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

clear all
close all

%1. Load the data 
% =================================================
load('IL27_NewData');  % loading data  % 1pt
tt=DataT(end); % tt is the total duration of Data % 1pt
tspan=0:0.1:tt;  % unit is min  % 1pt

%2. Set Parameters setting
% =================================================
NOP=8;  % number of old para: muC, muT, mug, CM , kc, kg, k27, I27  % 1pt
NNP=7;  % number of new para: lambdaC, lambdaTC,lambdaTIg, lambdaTI27, lambdaIgT, eta, Sg % 1pt
TPara=NOP+NNP;
VN=3; % variable number % 1pt
AllNum=[NOP,NNP,TPara,VN];

M=500; % Generation times % 1pt
PN=10; % Parent number % 1pt
CN=5; % Children number % 1pt
Parents=zeros(PN,TPara); % 1pt
Children=zeros(PN*CN,TPara); % 1pt

mu=0.1;     % percentage of Parent to have macromutation % 1pt
lambda=0.05; % max possible precent change of a parameter % 1pt
beta=0.02;     % effective energy % 1pt

AllParents=zeros(M,PN,TPara);

%3. Define the first generation
% =================================================
Para=zeros(1,4); % parameter set  % 1pt for all parameter values
Para(1)=0.519; %muC
Para(2)=0.9; %muT
Para(3)=11.04; %mug
Para(4)=1.8127; % CM
Para(5:8)=[1,1,1,1]; %kc, kg, k27, I27

% Para = [muC, muT, mug, CM , kc, kg, k27, I27]

for pn=1:PN  % Give parameters
    Parents(pn,1:NOP)=Para(1:NOP); % Assign muC, muT, mug, CM , kc, kg, k27, I27 % 1pt         
    Parents(pn,NOP+1)=2*rand(1)+3; % Initial guest of lambdaC % 2 pts for NOP-NOP+MNNP
    Parents(pn,NOP+2)=rand(1); % Initial guest of lambdaTC 
    Parents(pn,NOP+3)=rand(1)+1; % Initial guest of lambdaTIg 
    Parents(pn,NOP+4)=rand(1); % Initial guest of lambdaTI27 
    Parents(pn,NOP+5)=3*rand(1)+3; % Initial guest of lambdaIgT 
    Parents(pn,NOP+6)=rand(1); % Initial guest of eta 
    Parents(pn,NOP+7)=2*rand(1)+3; % Initial guest of Sg 
end

% Give IC for (T,I,V)

   init_IC=[1,1,1]; % 1pt
%4. Start Generation
% =================================================
for Gen=1:M
fprintf('Gen=%d \n', Gen)
  % Macromutation for the 1st Parents
    PopuIndex=[1:PN];
    SamInd=randsample(PopuIndex, mu*PN);
    for k=1:mu*PN
        MutUpRate=1.05;
        MutLowRate=0.95; 
        for Set=NOP+1:NOP+NNP % 1pt
            Parents(SamInd(k),Set)=rand(1)*(Parents(SamInd(k),Set)*MutUpRate-Parents(SamInd(k),Set)*MutLowRate)+Parents(SamInd(k),Set)*MutLowRate;
        end 
    end

    % Generate & score Children
    for pn=1:PN
        for cn=1:CN
            for Set=1:NOP 
                Children((pn-1)*CN+cn,Set)=Parents(pn,Set);
            end 
            for Set=NOP+1:NOP+NNP % 1pt
                delta=rand(1);
                Children((pn-1)*CN+cn,Set)=Parents(pn,Set)*(1+lambda*(delta-0.5*1));
            end
            % NewC_IC=CTIgmodel_restIC(Children((pn-1)*CN+cn,:),init_IC);
            C_sol = CTIgmodel_min(AllNum,Children((pn-1)*CN+cn,:),init_IC,tspan);% estimate LL % 1pt
            ChildrenIndex((pn-1)*CN+cn)=(pn-1)*CN+cn;
            ChilSco((pn-1)*CN+cn)=C_sol;
            theta((pn-1)*CN+cn)=exp(-beta*ChilSco((pn-1)*CN+cn)); % 1pt
        end

    end

    % Score & reorder Parents 
    for pn=1:PN
        % NewP_IC=CTIgmodel_restIC(Parents(pn,:),init_IC);  
        P_sol = CTIgmodel_min(AllNum,Parents(pn,:),init_IC,tspan); % 1pt
        PareSco(pn)=P_sol;
    end
    [Order_PareSco,P_INDEX]=sort(PareSco); % 1pt

    for pn=1:PN
        ReOderParents(pn,:)=Parents(P_INDEX(pn),:);    
    end

    % Select the Parents for next genation
     % 30% are chosen from the old Parents with top score
     for pn=1:PN*0.3 
        Parents(pn,:)=ReOderParents(pn,:);   % 1pt
        ParentsScore(pn)=Order_PareSco(pn);
     end

     % 70% are chose from children by theta
     New_theta=theta;
     for pn=1:PN*0.7
         if New_theta~=zeros(size(New_theta)) % 5pts
            FindIndex(pn)=randsample(ChildrenIndex, 1, true, New_theta);
         else
             FindIndex(pn)=randsample(ChildrenIndex,1);
         end
         New_theta(FindIndex(pn))=0; % remove from the current index
     end

     for pn=1:PN*0.7 
        Parents(pn+PN*0.3,:)=Children(FindIndex(pn),:); % 1pt
        ParentsScore(pn+PN*0.3)=ChilSco(FindIndex(pn));
     end

     % Recore the min score at each generation
     % ParentsScore
     MinScore(Gen)=min(ParentsScore) % 1pt
     %pause
     save('EvolutionaryAssignmentData','Gen','Parents','ParentsScore','MinScore')


end


save('EvolutionaryAssignmentData','Gen','Parents','ParentsScore','MinScore')

params=Parents(1,:);
options = odeset('RelTol',1e-6);% SET THE TOLERENCE OF THE SOLVER
duration=tspan(end); 
[sol.x,sol.y] = ode45(CTIgmodel(params),tspan,init_IC,options);

figure(1)
plot(1:length(MinScore), MinScore,'r','LineWidth',1), hold on
xlabel('Generation')
ylabel('Total Error')
legend('MinScore')
figure(2)
plot(DataT',Data(:,1),'g*'), hold on
plot(sol.x(:),sol.y(:,1),'g','LineWidth',1), hold on
plot(DataT',Data(:,2),'r*'), hold on 
plot(sol.x(:),sol.y(:,2),'r','LineWidth',1), hold on
plot(DataT',Data(:,3),'b*'), hold on 
plot(sol.x(:),sol.y(:,3),'b','LineWidth',1), hold on 
legend('Data C','fit C','Data T','fit T','Data IFN','fit IFN')


% Total = 47