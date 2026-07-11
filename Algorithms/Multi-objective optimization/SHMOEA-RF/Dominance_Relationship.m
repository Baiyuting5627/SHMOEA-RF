function [x] = Dominance_Relationship(a,b,m)
%基于目标的支配关系
%Input:
%   a,b-两个目标向量，如1*2
%   c-code length
%   m-objective number
%Output
%   x-dominance relationship
%       1-a dominates b ,a比b小,最小化问题
%       2-b dominates a
%       3-a = b
%       4-a and b are non-dominated by each other

t=0;
q=0;
p=0;
for i=1:m
    if a(1,i)<=b(1,i)
        t=t+1;
    end
    if a(1,i)>=b(1,i)
        q=q+1;
    end
    if a(1,i)==b(1,i)
        p=p+1;
    end
end
if t==m&p~=m
    x=1;
elseif q==m&p~=m
    x=2;
elseif p==m
    x=3;
else
    x=4;
end
end

