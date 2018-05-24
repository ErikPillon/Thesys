function [eff_H] = f_disp_diag(s,site,d,States,v,eff_H,nSt,fds,nc)
% This function is intended to modify the effective Hamiltonian.
% 
% A significative part of the code is dedicated to the implementation of 
% the Slater determinant
%
% Idea:
% [ 1 0 1 0 0 1 0 1 1 0 1 0 0 1 0 1]
% is described by 
% c1 c3 c6 c8 c9 c11 c14 c16
% and then if we want annihilate 6 and create in 10 we have to skip 
% initially 3, ie, (-1)^3 and then 4 sites, ie, (-1)^4 
% for obtaining
% [ 1 0 1 0 0 0 0 1 1 1 1 0 0 1 0 1]
% 
%
% first part of the Slater determinant

%% Slater determinant
num = sum(s(1:site-1));
dims = sqrt(length(s));
x = 1+mod(site-1,dims);
y = ceil(site/dims);
xn = 1+mod(x+d(2)-1,dims);
yn = 1+mod(y+d(1)-1,dims);
new_site = (yn-1)*dims+xn;

% create the new state
r = f_swap(s,site,d);   

num = num+sum(r(1:new_site-1)); % second part of Slater determinant

%% add v ot -v
for pos = [1:nSt]
    if r == States(pos,1:end-2)
        eff_H(fds,pos) = v*(-1)^(num);
    end
end