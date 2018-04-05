% main reference for the work: moleskine 4 avril '18
%
% algorithm for the evaluation of a displacement on a tirangular half
% filled lattice
%
% the function f_two_excitations_energy_gap will be used

clc; clear all; close all;

delta = [1,1]/2;
delta2 = [-1,-1]/2;
site = [1,0];
ratio = sqrt(3);
theta = pi/2;
alpha = 1
% alpha = 0.5;
disp('displacement in')
disp(site)
disp('of delta')
disp(delta2)
[gap1, gap2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha)


delta2 = [-1,-1]/2;
site1 = [0,1];
delta3 = [1,1]/2;
site2 = [1,0];
disp('displacement in')
disp(site1)
disp('of delta')
disp(delta2)
disp('and displacement in')
disp(site2)
disp('of delta')
disp(delta3)
[gap1, gap2, gap3] = f_3_excitations(delta, delta2, site1, delta3, site2, ratio, theta, alpha) 

% %% other possible configurations
% delta2 = [-1,-1]/2;
% site = [0,1];
% disp('displacement in')
% disp(site)
% disp('of delta')
% disp(delta2)
% [gap1, gap2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha)
% 
% delta2 = [1,-1]/2;
% site = [1,0];
% disp('displacement in')
% disp(site)
% disp('of delta')
% disp(delta2)
% [gap1, gap2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha)
% 
% delta2 = [-1,-1]/2;
% site = [1,0];
% disp('displacement in')
% disp(site)
% disp('of delta')
% disp(delta2)
% [gap1, gap2] = f_two_excitations_energy_gap(delta, delta2, site, ratio, theta, alpha)
% 
% %% three prossible excitations
% delta2 = [-1,-1]/2;
% site1 = [1,0];
% delta3 = [-1,-1]/2;
% site2 = [-1,1];
% disp('displacement in')
% disp(site1)
% disp('of delta')
% disp(delta2)
% disp('and displacement in')
% disp(site2)
% disp('of delta')
% disp(delta3)
% [gap1, gap2, gap3] = f_3_excitations(delta, delta2, site1, delta3, site2, ratio, theta, alpha) 
% 
% delta2 = [-1,-1]/2;
% site1 = [0,1];
% delta3 = [1,1]/2;
% site2 = [1,0];
% disp('displacement in')
% disp(site1)
% disp('of delta')
% disp(delta2)
% disp('and displacement in')
% disp(site2)
% disp('of delta')
% disp(delta3)
% [gap1, gap2, gap3] = f_3_excitations(delta, delta2, site1, delta3, site2, ratio, theta, alpha) 