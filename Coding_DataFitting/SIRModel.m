function []=SIRModel



close all
clear all

global A nu gamma mu beta

A=2;
nu=1;
gamma=1;
mu=1;
beta=3;
% R0=beta/(nu+mu)

S0=1.8;
I0=0.2;
R0=0;

init = [S0,I0, R0];
tspan=[0,10];

[t,v]=ode45(@(t,v) fun_SIR(t,v), tspan, init);

save('SIR_Data','t','v')

figure(1)
plot(t,v(:,1),'g*'), hold on  
plot(t,v(:,2),'r*'), hold on 
plot(t,v(:,3),'b*'), hold on 
legend('S','I','R')

function dv=fun_SIR(t,v)
global A nu gamma mu beta
S=v(1);
I=v(2);
R=v(3);
dv=zeros(3,1);

% S> I > R -> S
dv(1)=A-beta*S*I+gamma*R-mu*S;
dv(2)=beta*S*I-nu*I-mu*I;
dv(3)=nu*I-gamma*R-mu*R;

