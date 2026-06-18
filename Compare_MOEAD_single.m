classdef Compare_MOEAD_single < ALGORITHM
% <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II 
% 单独 MOEAD + 全部单次评估 + 贪婪修正，最后一次多次评估

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)                                   
            %% Parameter setting
            num_Samples = 10;                                      % 真实评估采样次数   
            %% Generate the weight vectors
            [WV,Problem.N] = UniformPoint(Problem.N,Problem.M);   % 生成100*2的均匀分布的权重向量W
            T =  ceil(Problem.N/10);                              % 邻居个数，返回大于等于的最小整数，eg.,10           
            %% Detect the neighbours of each solution
            B = pdist2(WV,WV);     % 计算两组权重之间两两成对的距离
            [~,B] = sort(B,2);     % eg.,排序是升序
            B = B(:,1:T);          % 提取每个解和其他解间权重距离最小的前10个解所在的索引，100*10          
            %% Generate random population
            Population = Problem.Initialization();                               % 初始化种群
            [Population] = Greedy_Repair_Scalar(WV,Problem, Population);         % 贪婪修正
            Z = min(Population.objs,[],1);               % MOKP问题目标是最小值，理想点
            [FrontNo,~] = NDSort(Population.objs,inf);   % 对种群进行非支配排序，inf表示对整个总体进行排序。
            A = Population(FrontNo == 1);                % 选出初始种群第一层的解作为非支配解集          
            %% Optimization
            while Algorithm.NotTerminated(A)
                for i = 1 : Problem.N               % 遍历种群中的每个解
                    P = B(i,randperm(size(B,2)));   % 将每个解的邻居顺序打乱并提取                           
                    Offspring = OperatorGAhalf(Problem,Population(P(1:2)));           % 交叉变异
                    [Offspring] = Greedy_Repair_Scalar(WV,Problem, Offspring);        % 启发修正  
                    A = [A, Offspring];                 % 将后代解并入非支配解集
                    [FrontNo,~] = NDSort(A.objs,inf);   % 非支配排序
                    A = A(FrontNo == 1);                % 更新非支配解集                 
                    Z = min(Z,Offspring.obj);           % 更新全局最优解                    
                    g_old = sum(Population(P).objs.*WV(P,:),2);   % 加权聚合函数，eg.,g=lamda*f,10*1
                    g_new = sum(Offspring.obj.*WV(P,:),2);        % 计算后代解的加权聚合函数值
                    Population(P(g_old>=g_new)) = Offspring;      % 好解替换种群差解
                end
            end
        end
    end
end
