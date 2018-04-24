function Energy = f_lr_plus_nn(delta, ratio, theta, x, alpha, c) 
% this function is suited for evaluating the potential described in the
% paper "Glassy Dynamics in Geometrically Frustrated Coulomb Liquids 
% without Disorder", that is V = (1-x)*V_nn + x*V_lr
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
% delta : displacement (2D) 
% ratio : ratio between the two basis vectors
% theta : angle between the two basis vectors
% x     : convergence ratio between long range and nearest neighbor
% alpha : exponent of convergence 
% c     : convergence coefficient; in general the algorithm should not depend
%         on this coefficient
%
% --- OUTPUT ---
% energy_gap : difference in term of potential energy due to the
%              displacement delta
%
% --- NOTE   ---
% x = 0: nearest neighbor model
% x = 1: longrange ordering

% Check inputs
if nargin < 6
  c = 10;
  if nargin < 5
    alpha = 1;
    if nargin < 4
        ratio = 1;
        theta = pi/2;
      if nargin < 1
        error(message('f_one_excitation_energy_gap:not enough input arguments'));
      end  
    end
  end
end

% long range contribution
lr = f_one_excitation_energy_gap(delta, ratio, theta, alpha, c);

% nearest neighbor contribution
nn = f_one_excitation_energy_gap(delta, ratio, theta, 40, c);

Energy = x*lr+(1-x)*nn;






