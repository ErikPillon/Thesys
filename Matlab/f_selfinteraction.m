function [A] = f_selfinteraction(delta, c, alpha)
% [A] = selfinteraction(delta)
% [A] = f_selfinteraction(delta, c, alpha)
% 
% --- INPUT  ---
% delta : argument of the function of selfinteraction
% c     : convergence coefficient
% alpha : exponent of the potential
%
% --- OUTPUT  ---
% A : energy of selfinteraction

if delta == 0
    A = 2*c^(alpha/2)/alpha;
else
    A = gamma(alpha/2)*gammainc(c*delta^2,alpha/2)/(delta^alpha);
end

end

    