function [energy] = f_phi(p,a)
% [energy] = f_phi(p,a)
% 
% this function is intended to evaluate the following integral
%              inf
%              /
% phi(p, a) =  | exp(-at) t^p dt
%              /
%             t=1
%
% NOTE: in the particular case in which the function is evaluated with a=0
% the default value is 0 for all the function

f = @(t) exp(-a.*t).*t.^p;

energy = 0; % default value

if a == 0
    return
else
    energy = integral(f,1,inf);
end

end
