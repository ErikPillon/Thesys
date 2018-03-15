function [energy_gap] = f_one_excitation_energy_gap(delta, ratio, theta, alpha, c) 
% 
% [energy_gap] = f_one_excitation_energy_gap(delta) 
% [energy_gap] = f_one_excitation_energy_gap(delta, ratio, theta) 
% [energy_gap] = f_one_excitation_energy_gap(delta, ratio, theta, alpha)
% [energy_gap] = f_one_excitation_energy_gap(delta, ratio, theta, alpha, c) 
%
% function for the implementation of the ewald summation following 
% notes on pages 3 and following
% suggested reference "Ewald-Rastelli"
%
% --- INPUT  ---
% 
% c : convergence coefficient; in general the algorithm should not depend
%     on this coefficient
%
% --- OUTPUT ---
% energy_gap : difference in term of potential energy due to the
%              displacement delta

% Check inputs
if nargin < 5
  c = 1;
  if nargin < 4
    alpha = 1;
    if nargin < 3
        ratio = 1;
        theta = pi/2
      if nargin < 1
        error(message('f_one_excitation_energy_gap:not enough input arguments'));
      end  
    end
  end
end

%% initialize all the variables
vol = ratio*sin(theta);

G_cut = 20;
R_cut = 20;

% the final summation is given by [sum1-sum2+sum3-sum4]/gammma(alpha/2);
sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;

%% evaluation of all the terms
% sum1
for Gx = -G_cut:G_cut
    for Gy = -G_cut:G_cut
        G = [Gx, Gy];
        sum1 = sum1+((cos(G*delta')-1)*f_phi(-alpha/2,G*G'/4/c));
    end
end
sum1 = sum1*pi*c^(alpha/2-1)/vol;

% sum2
sum2 = f_selfinteraction(norm(delta), c, alpha)-f_selfinteraction(0, c, alpha);

% sum3
for Rx = -R_cut:R_cut
    for Ry = -R_cut:R_cut
        R = [Rx, Ry];
        sum3 = sum3+(f_phi(alpha/2-1,norm(R+delta)^2*c)-f_phi(alpha/2-1,norm(R)^2*c));
    end
end
sum3 = sum3*c^(alpha/2);

% sum4
sum4 = f_phi(alpha/2-1,norm(delta)^2*c)-f_phi(alpha/2-1,0);
sum4 = c^(alpha/2)*sum4;

energy_gap = (sum1-sum2+sum3-sum4)/gamma(alpha/2);

end
