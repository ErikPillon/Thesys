function [new_state] = f_switch(pos,delta,old_state,new_configuration)
% 
% 
% 
% 
% 
% 
% 
% 
n = sqrt(length(old_state));

x = pos(1); y = pos(2);
xn = 1+mod(x+delta(1)-1,n); yn = 1+mod(y+delta(2)-1,n);

Matrix = reshape(old_state,n,n);

if Matrix(x,y) == 0 
    error('f_switch:the desired state is empty');
elseif Matrix(xn,yn) == 1  
    error('f_switch:the desired state is already occupied');
else
    Matrix(x,y) = 0;
    Matrix(xn,yn) = 1;
end

one_move = reshape(Matrix,1,n^2);
new_state= [one_move,[xn,yn],new_configuration];


