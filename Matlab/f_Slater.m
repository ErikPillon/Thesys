function [sign] = f_Slater(State_unperturbed,site,new_State,displacement)
% 
% 
% 
% 
% 
% 
% 
% 
% n = sqrt(length(State_unperturbed));
% 
% x = site(1); y = site(2);
% xn = 1+mod(x+displacement(1)-1,n); 
% yn = 1+mod(y+displacement(2)-1,n);
% 
% count = sum(State_unperturbed(1:(y-1)*n+x),2);
% count = count+sum(new_State(1:(yn-1)*n+xn),2);
% 
% sign = (-1)^count;
sign = 1;