x = 1:50;
X = reshape(x,[5,10])';

%r = 1:5;
r = 2:6;
flag = 0;
pos = 1;
for l=1:length(X(:,1))
    if r == X(l,:)
        flag = 1;
        break
    end
    pos = pos+1;
end

if flag == 0
    X(end+1,:) = r;
end











