function [gap1, gap2, gap3] = f_3_excitations(delta, delta2, site1, delta3, site2, ratio, theta, alpha, c) 
% 
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site1, delta3) 
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site1, delta3, site2, ratio, theta) 
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site1, delta3, site2, ratio, theta, alpha)
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site1, delta3, site2, ratio, theta, alpha, c) 
%
% function for the implementation of the ewald summation following 
% notes on pages 3 and following
% suggested reference "Ewald-Rastelli"
%
% --- INPUT  ---
% delta : displacement (2D) 
% delta2: second displacement (2D)
% site1 : cell (i,j) in which we apply the second displacement 
% delta3: third displacement (2D)
% site2 : cell (i,j) in which we apply the third displacement 
% ratio : ratio between the two basis vectors
% theta : angle between the two basis vectors
% alpha : exponent of convergence 
% c     : convergence coefficient;  the algorithm should not depend on this 
%         parameter
% 
% --- OUTPUT ---
% gap1 : energy gap due to the first displacement
% gap2 : energy gap due to the second displacement
% gap3 : energy gap due to the third displacement

% Check inputs
if nargin < 9
  c = 10;
  if nargin < 8
    alpha = 1;
    if nargin < 7
        ratio = 1;
        theta = pi/2;
      if nargin < 5
        error(message('f_3_excitations:not enough input arguments'));
      end  
    end
  end
end

% the two energy gap can be evaluated with the function for only one displacement.
% We need only to take into account the correction to this configuration
[gap1, gap2] = f_two_excitations_energy_gap(delta, delta2, site1, ratio, theta, alpha, c);
energy3 = f_one_excitation_energy_gap(delta3, ratio, theta, alpha);
%% evaluation of the energy gap of the two coupled displacements
% the idea is to subtract two "fake" interactions  
% + add the two real "interactions"

% intialize the variables
diff1_f = 0; diff1_r = 0; diff2_f = 0; diff2_r = 0;
diff3_f = 0; diff3_r = 0; diff4_f = 0; diff4_r = 0;

% at the end we will need to add the quantity 
% diff1_r + diff2_r - diff1_f - diff2_f +
% diff3_r + diff4_r - diff3_f - diff4_f

a = [1, 0];
rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
b = (ratio*rot*a')';

dist1 = site1(1)*a + site1(2)*b;
dist2 = site2(1)*a + site2(2)*b;

% the displacements must be converted in the lattice coordinate
delta = delta(1)*a+delta(2)*b;
delta2 = delta2(1)*a+delta2(2)*b;
delta3 = delta3(1)*a+delta3(2)*b;

% evaluate directly the energy difference
diff1_f = 1/norm(dist2-delta)^alpha;
diff1_r = 1/norm(dist2-delta+delta3)^alpha;
diff2_f = 1/norm(-delta2-dist1+dist2)^alpha;
diff2_r = 1/norm(-delta2-dist1+dist2+delta3)^alpha;
diff3_f = 1/norm(-dist2-delta3)^alpha;
diff3_r = 1/norm(-dist2-delta3+delta)^alpha;
diff4_f = 1/norm(-delta3-dist2+dist1)^alpha;
diff4_r = 1/norm(-delta3-dist2+dist1+delta2)^alpha;

gap3 = gap2+energy3-diff1_f+diff1_r-diff2_f+diff2_r-diff3_f+diff3_r-diff4_f+diff4_r;
end