function [gap1, energy_gap] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha, c) 
% 
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2) 
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta) 
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha)
% [energy_gap] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha, c) 
%
% function for the implementation of the ewald summation following 
% notes on pages 3 and following
% suggested reference "Ewald-Rastelli"
%
% --- INPUT  ---
% delta : displacement (2D) 
% delta2: second displacement (2D)
% site  : cell (i,j) in which we apply the second displacement 
% ratio : ratio between the two basis vectors
% theta : angle between the two basis vectors
% alpha : exponent of convergence 
% c     : convergence coefficient; in general the algorithm should not depend
%         on this coefficient
%
% --- OUTPUT ---
% energy_gap : difference in term of potential energy due to the
%              displacement delta

% Check inputs
if nargin < 7
  c = 10;
  if nargin < 6
    alpha = 1;
    if nargin < 5
        ratio = 1;
        theta = pi/2
      if nargin < 3
        error(message('f_one_excitation_energy_gap:not enough input arguments'));
      end  
    end
  end
end

% the two energy gap can be evaluated with the function for only one displacement.
% We need only to take into account the correction to this configuration
gap1 = f_one_excitation_energy_gap(delta, ratio, theta, alpha);
gap2 = f_one_excitation_energy_gap(delta2, ratio, theta, alpha);

%% evaluation of the energy gap of the two coupled displacements
% the idea is to subtract two "fake" interactions  
% + add the two real "interactions"

% intialize the variables
diff1_f = 0; diff1_r = 0; diff2_f = 0; diff2_r = 0;
% at the end we will need to add the quantity 
% diff1_r + diff2_r - diff1_f  - diff2_f 

a = [1, 0];
rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
b = (ratio*rot*a')';

dist = site(1)*a + site(2)*b;

% the displacements must be converted in the lattice coordinate
delta = delta(1)*a+delta(2)*b;
delta2 = delta2(1)*a+delta2(2)*b;

% evaluate directly the energy difference
diff1_f = 1/norm(dist-delta)^alpha;
diff1_r = 1/norm(dist-delta+delta2)^alpha;
diff2_f = 1/norm(-delta2-dist)^alpha;
diff2_r = 1/norm(-delta2-dist+delta)^alpha;

energy_gap = gap1+gap2-diff1_f+diff1_r-diff2_f+diff2_r;
end