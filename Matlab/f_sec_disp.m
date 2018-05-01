function [States2, eff_H] = f_sec_disp(s,site,d,nc,States2,t,eff_H,nSt,fds)
r = f_swap(s,site,d);
pos = 1;
flag = 0;
for l=1:length(States2(:,1))
    if r == States2(l,1:end-2)
        % SLATER DETERMINANT TO BE ADDED
        eff_H(fds,nSt+l) = t;
        flag = 1;
        break
    end
    pos = pos+1;
end
pos = pos-1;

if flag == 0
    States2(end+1,:) = [r,site,nc];
    eff_H(fds,nSt+pos) = t;
end






