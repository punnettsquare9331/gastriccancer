function []=DesignData

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

Data=zeros(4,5);
Data(1,:)=[1,10,20,30,40]; % time series
Data(2,:)=[10, 20, 40, 80, 160]; % value of T
Data(3,:)=[11, 21, 41, 81, 161]; % value of I
Data(4,:)=[12, 22, 42, 82, 162]; % value of V
save('Data','Data')

