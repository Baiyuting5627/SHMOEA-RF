function [ POP ] = Select_Solution_NSGA(Selected_Pop,predict_objs,num_selet)
% Input:
% Selected_Pop - 初步筛选出的 SOLUTION 好解 
% predict_objs  - 由RF预测后的鲁棒多次平均评估采样的目标函数值
% num_selet - 进一步选择的解的个数，eg.,25
% Output
% POP - 进一步选择出的 SOLUTION 解

Objectives = predict_objs;                 % 由RF预测后的多次平均评估采样的目标函数值
[~, isPareto] = paretofront(Objectives);   % 对预测值进行非支配排序，得到帕累托最优值
non_Dominated_Indices = find(isPareto);    % 选择帕累托最优值所在的行索引

num_NonDominated = length(non_Dominated_Indices);  % 找到的帕累托最优值的数量
if num_NonDominated >= num_selet                   % 如果找到的帕累托最优值的数量大于或等于num_selet   
   selected_Indices = non_Dominated_Indices(1:num_selet);  % 只选择前num_selet个索引
else                                                       % 如果找到的帕累托最优值的数量小于num_selet，
    selected_Indices = non_Dominated_Indices;      % 选择所有已找到的帕累托最优值        
    dominated_Indices = ~isPareto;                 % 提取所有帕累托支配解的索引，eg.,1*42
    dominated_Solutions = Objectives(dominated_Indices, :);  % 提取所有非帕累托最优解的目标值
    sum_Objectives = sum(dominated_Solutions, 2);            % 计算所有非帕累托最优解的目标值和
    [~, sort_Index] = sort(sum_Objectives);                  % 对所有非帕累托最优解的目标值之和排序，得到排序后的索引sortIndex，eg.,42*1
    sort_Index = sort_Index'; 
    num_To_Select = num_selet - num_NonDominated;                       % 目前找到的好解与期望选择好解的数量num_selet的差
    selected_Indices = [selected_Indices, sort_Index(1:num_To_Select)]; % 按照顺序逐个选择非帕累托最优解目标值之和最小的解，直到总数达到25个
end
for i=1:size(Selected_Pop,2)            
    Selected_Pop(i) = SOLUTION(Selected_Pop(i).dec,predict_objs(i,:),Selected_Pop(i).con);  % 预测值放入种群构成SOLUTION解
end
POP = Selected_Pop(:,selected_Indices);                                                     % 进一步选择出的好解
end

function [front, isPareto] = paretofront(M)
    % PARETOFRONT returns the logical Pareto Front of a set of points.
    % 
    % synopsis: [front, isPareto] = paretofront(M)
    % 
    % INPUT ARGUMENT
    % - M n x m array, of which (i,j) element is the j-th objective value of the i-th point;
    % 目标值
    % OUTPUT ARGUMENT
    % - front: indices of the Pareto front points
    % - isPareto: logical array indicating whether each point is on the Pareto front

    [n, m] = size(M);
    isPareto = true(1, n);
    front = [];

    for i = 1:n
        for j = 1:n
            if i ~= j
                if all(M(j,:) <= M(i,:)) && any(M(j,:) < M(i,:))
                    isPareto(i) = false;
                    break;
                end
            end
        end
        if isPareto(i)
            front = [front, i];
        end
    end
end