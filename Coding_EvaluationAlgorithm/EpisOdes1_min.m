function AllError = EpisOdes1_min(AllNum,ParaStart,IC,tspan)

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

NOP  = AllNum(1); % number of old para
NNP  = AllNum(2);  % number of new para
TPara= AllNum(3);
VN   = AllNum(4);

params = ParaStart;
load('Data');

% ======  Part I fitting =============
options = odeset('RelTol',1e-6);% SET THE TOLERENCE OF THE SOLVER
duration=tspan(end);               
sol = ode15s(EpisOdes_fit1(params),tspan,IC,options);% SOLVE THE ODE

Error=zeros(1,3);
AllError=0;  
New=zeros(VN,length(tspan));
for j=1:VN
    New(j,1:length(tspan))=spline(sol.x(:), sol.y(j,:),tspan);   
    Error(j)=0;
    for q=1:length(Data(1,:))
        Error(j)=Error(j)+(New(j,Data(1,q))-Data(j+1,q))^2;
    end
    
    AllError=AllError+Error(j);
end

        
            
    


