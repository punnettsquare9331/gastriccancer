function [sol]=EvaluationAlgorithm1

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

clear all
close all

%1. Load the data 
% =================================================
load('Data');  % loading data 
tt=Data(1,end); % tt is the total duration of Data 
tspan=1:tt*1;  % unit is min
 
%2. Set Parameters setting
% =================================================
NOP=2;  % number of old para: R and fq
NNP=4;  % number of new para: beta, gamma
TPara=NOP+NNP;
VN=3; % variable number 
AllNum=[NOP,NNP,TPara,VN];

M=12; % Generation times
PN=10; % Parent number
CN=5; % Children number
Parents=zeros(PN,TPara);
Children=zeros(PN*CN,TPara);

mu=0.4;     % percentage of Parent to have macromutation
lambda=0.5; % max possible precent change of a parameter
beta=2;     % effective energy

AllParents=zeros(M,PN,TPara);

%3. Define the first generation
% =================================================
Para=zeros(1,4); % parameter set 
Para(1)=10; %A
Para(2)=2; %mu1
Para(3)=4; %mu2
Para(4)=5; % K

% para = [A,  mu1, mu2,k, beta,gamma]

for pn=1:PN  % Give parameters
    Parents(pn,1:NOP)=Para(1:NOP); % Assign A, mu1, mu2, k  
    for Set=NOP+1:NOP+NNP          
        Parents(pn,Set)=1000*rand(1); % Initial guest of beta, gamma
    end    
end


% Give IC for (T,I,V)

   init_IC=[0.1,0.1,0.1];
%4. Start Generation
% =================================================
for Gen=1:M
fprintf('Gen=%d \n', Gen)

  % Macromutation for the 1st Parents
    PopuIndex=[1:PN];
    SamInd=randsample(PopuIndex, mu*PN);
    for k=1:mu*PN
        MutUpRate=100;
        MutLowRate=0; 
        for Set=NOP+1:NOP+NNP 
            Parents(SamInd(k),Set)=rand(1)*(Parents(SamInd(k),Set)*MutUpRate-Parents(SamInd(k),Set)*MutLowRate)+Parents(SamInd(k),Set)*MutLowRate;
        end 
    end

    % Generate & score Children
    for pn=1:PN
        for cn=1:CN
            for Set=1:NOP 
                Children((pn-1)*CN+cn,Set)=Parents(pn,Set);
            end 
            for Set=NOP+1:NOP+NNP 
                delta=rand(1);
                Children((pn-1)*CN+cn,Set)=Parents(pn,Set)*(1+lambda*(delta-0.5*1));
            end
            C_sol = EpisOdes1_min(AllNum,Children((pn-1)*CN+cn,:),init_IC,tspan);% estimate LL
            
            ChildrenIndex((pn-1)*CN+cn)=(pn-1)*CN+cn;
            ChilSco((pn-1)*CN+cn)=C_sol;
            theta((pn-1)*CN+cn)=exp(-beta*ChilSco((pn-1)*CN+cn));
        end
    end
    
    % Score & reorder Parents 
    for pn=1:PN
        P_sol = EpisOdes1_min(AllNum,Parents(pn,:),init_IC,tspan);
        PareSco(pn)=P_sol;
    end
    [Order_PareSco,P_INDEX]=sort(PareSco);

    for pn=1:PN
        ReOderParents(pn,:)=Parents(P_INDEX(pn),:);     
    end

    % Select the Parents for next genation
     % 10% are chosen from the old Parents with top score
     for pn=1:PN*0.1 
        Parents(pn,:)=ReOderParents(pn,:);
        ParentsScore(pn)=Order_PareSco(pn);
     end
     
     % 90% are chose from children by theta
     New_theta=theta;
     for pn=1:PN*0.9
          FindIndex(pn)=randsample(ChildrenIndex, 1, true, New_theta);
          New_theta(FindIndex(pn))=0; % remove from the current index
     end
      
     for pn=1:PN*0.9 
        Parents(pn+PN*0.1,:)=Children(FindIndex(pn),:);
        ParentsScore(pn+PN*0.1)=ChilSco(FindIndex(pn));
     end
    
     % Recore the min score at each generation
     % ParentsScore
     MinScore(Gen)=min(ParentsScore)
     %pause
    
     fpath1=['EvoluAlgor/EA_EpisModel1_P1.mat'];
     save(fpath1,'Gen','Parents','ParentsScore','MinScore')
    

end

fpath1=['EvoluAlgor/EA_EpisModel1_P1.mat'];
save(fpath1,'Gen','Parents','ParentsScore','MinScore')








