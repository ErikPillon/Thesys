function [energy] = f_phi(p,a)
% [energy] = f_phi(p,a)
% 
% this function is intended to evaluate the following integral
%                 inf
%                 /
% phi(p, a) =     | exp(-at) t^p dt
%                 /
%                t=1

f = @(t) exp(-a.*t).*y^p;

energy = integral(f,1,inf)