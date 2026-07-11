function [ NPOP ] = Nondominated( POP,c,m)
%基于目标的非支配排序
%Input:
%   m-objective number
%   c-code length 问题实例个数250
%   POP- 10*254
%Output
%   NPOP-selected non-dominated individuals
pop = [POP.decs,POP.objs,POP.cons];
n=size(pop,1);%可行解个数,5
Objs=pop(:,c+1:c+m);%提取目标函数列251:252列
%Objs=POP.objs;
np=zeros(n,1);
NPOP=[];
for i=1:n
    for j=1:n
        if j~=i
            x=Dominance_Relationship(Objs(i,:),Objs(j,:),m);
            if x==2 %j支配i，j比i好，i被j支配，一个解被认为是非支配的，如果不存在另一个解在所有目标上都更优（即没有任何目标被牺牲）
                np(i)=np(i)+1;
            end
        end
    end
    if np(i)==0 % i是非支配的，也可能x=1,3,4,
       %NPOP=[NPOP;pop(i,:)];
       NPOP=[NPOP,POP(1,i)];
    end
end
end
% x-dominance relationship
%  1-a dominates b，a比b小
%  2-b dominates a,a比b大
%  3-a = b
%  4-a and b are non-dominated to each other