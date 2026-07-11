function [ NPOP ] = Nondominated_Cons( POP,c,m,nc)
%基于目标和约束的支配排序,正常解254列
%Input:
%   m-objective number2
%   c-code length 250
%   POP-population,100*254
%   nc-the number of constrains2
%Output
%   NPOP-selected non-dominated individuals
n=size(POP,1);%1000
Con =POP(:,c+m+1:c+m+nc);%约束列 253-254
np=zeros(n,1);
NPOP=[];
for i=1:n
    for j=1:n
        if j~=i
            %x=Dominance_Relationship_C(POP(i,:),POP(j,:),Con(i),Con(j),m,c);
            x=Dominance_Relationship_Cons(POP(i,:),POP(j,:),Con(i),Con(j),m,c);
            if x==2 %j支配i，j比i好，i被j支配，一个解被认为是非支配的，如果不存在另一个解在所有目标上都更优（即没有任何目标被牺牲）
                np(i)=np(i)+1;
            end
        end
    end
    if np(i)==0 % i是非支配的
        NPOP=[NPOP;POP(i,:)];
    end
end
end
% x-dominance relationship
%  1-a dominates b
%  2-b dominates a
%  3-a = b
%  4-a and b are non-dominated to each other
