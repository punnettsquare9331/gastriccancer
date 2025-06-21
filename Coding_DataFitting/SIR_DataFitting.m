% function []=SIR_DataFitting
% S'=A-beta*S*I+gamma*R-mu*S;
% I'=beta*S*I-nu*I-mu*I;
% R'=nu*I-gamma*R-mu*R;

close all
clear all

load('SIR_Data')
tspan=[0,10];
A=3; % 2
nu=1; % 1
gamma=11; % 1
mu=0.2; % 1
beta=4; % 3

xt0=zeros(1,8);
xt0(1)=A;
xt0(2)=nu;
xt0(3)=gamma;
xt0(4)=mu;
xt0(5)=beta;
xt0(6:8)=v(1,:);
[pbest,presnorm,presidual,exitflag,output] = lsqcurvefit(@paramfun,xt0,t,v);

fprintf('New parameters:')
pbest


%A=pbest(1);
%nu=pbest(2);
%gamma=pbest(3);
%mu=pbest(4);
%beta=pbest(5);
init = pbest(6:8);
tspan=[0,10];
[newt,newv]=ode45(@(t,v) fun_SIR(t,v,pbest), tspan, init);

figure(1)
plot(t,v(:,1),'g*'), hold on
plot(newt,newv(:,1),'g','LineWidth',1), hold on
plot(t,v(:,2),'r*'), hold on 
plot(newt,newv(:,2),'r','LineWidth',1), hold on
plot(t,v(:,3),'b*'), hold on 
plot(newt,newv(:,3),'b','LineWidth',1), hold on 
legend('Data S','fit S','Data I','fit I','Data R','fit R')

function pos=paramfun(x,tspan)
A=x(1);
nu=x(2);
gamma=x(3);
mu=x(4);
beta=x(5);
xt0=x(6:8);
f=@(t,y) [A-beta*y(1)*y(2)+gamma*y(3)-mu*y(1);...
          beta*y(1)*y(2)-nu*y(2)-mu*y(2);...
          nu*y(2)-gamma*y(3)-mu*y(3)];
      
[~,pos]=ode45(f,tspan,xt0);

end

function dv=fun_SIR(t,v,pbest)
global A nu gamma mu beta
S=v(1);
I=v(2);
R=v(3);
dv=zeros(3,1);

A=pbest(1);
nu=pbest(2);
gamma=pbest(3);
mu=pbest(4);
beta=pbest(5);

dv(1)=A-beta*S*I+gamma*R-mu*S;
dv(2)=beta*S*I-nu*I-mu*I;
dv(3)=nu*I-gamma*R-mu*R;
end

