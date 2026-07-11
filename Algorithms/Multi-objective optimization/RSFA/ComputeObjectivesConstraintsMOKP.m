% 
% function [obj, Constraints] = ComputeObjectivesConstraintsMOKP(obj, POP, m, c, nc)
%     con = sum(obj.W, 2)' / 2;            % 1×2
% 
%     % 计算目标值，使用临时变量 temp_obj
%     temp_obj = [];
%     for i = 1 : m
%         T = POP * obj.P(i, :)';          % obj.P 正常访问
%         temp_obj = [temp_obj, T];
%     end
%     temp_obj = 0 - temp_obj;              % 取负
% 
%     % 计算约束
%     Constraints = [];
%     for i = 1 : nc
%         T = POP * obj.W(i, :)' - con(i);
%         Constraints = [Constraints, T];
%     end
%     Constraints(Constraints < 0) = 0;     % 更简洁的写法
% 
%     obj = temp_obj;                       % 输出目标值
% end

function [obj, Constraints] = ComputeObjectivesConstraintsMOKP(Problem, POP, m, c, nc)
% ComputeObjectivesConstraintsMOKP - PlatEMO兼容的评估函数
% 
% 输入:
%   Problem    - PlatEMO问题对象（含有Evaluation方法）
%   POP        - 决策变量矩阵，尺寸 N x c
%   m          - 目标个数（可选，若不提供则从Problem.M获取）
%   c          - 决策变量维数（可选，若不提供则从Problem.D获取）
%   nc         - 约束个数（可选，若不提供则从Problem.cons获取）
%
% 输出:
%   obj        - 目标值矩阵，尺寸 N x m，已取负（将原最大化问题转为最小化）
%   Constraints- 约束值矩阵，尺寸 N x nc，已执行 max(0, constraints) 操作

    % 参数自动推断：如果未提供或为空，则从Problem中获取
    if nargin < 3 || isempty(m)
        m = Problem.M;
    end
    if nargin < 4 || isempty(c)
        c = Problem.D;
    end
    if nargin < 5 || isempty(nc)
        nc = Problem.cons;
    end

    % 调用PlatEMO的标准评估接口（自动累加评估次数）
    Population = Problem.Evaluation(POP);   % 返回SOLUTION对象数组

    % 提取目标值和约束值
    obj = Population.objs;          % N × m
    Constraints = Population.cons;  % N × nc

    % 原算法中目标值取负（将最大化问题转为最小化，适应进化框架）
    obj = -obj;

    % 原算法中约束违反值被修剪为非负数：负值置0
    % 注意：PlatEMO中约束通常以 ≤0 表示满足，此处转换为非负违反度
    Constraints(Constraints < 0) = 0;
end