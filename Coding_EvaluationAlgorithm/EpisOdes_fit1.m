function h = EpisOdes_fit1(para)

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

h = @fx;
    
A=para(1);
mu1=para(2);
mu2=para(3);
k=para(4);
beta=para(5);
gamma=para(6);

    function [yout] = fx(t,y)       
        yout = zeros(3,1);      
        yout(1) = A-beta*y(1)*y(3)-mu1*y(1); % T
        yout(2) = beta*y(1)*y(3)-mu2*y(2); % I
        yout(3) = gamma*mu2*y(2)-k*y(3); % V
              
    end
end


