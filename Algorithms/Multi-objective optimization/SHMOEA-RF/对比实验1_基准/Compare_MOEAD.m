classdef Compare_MOEAD < ALGORITHM
% <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
% 单独 MOEAD + 全部多次评估 + 贪婪修正

    methods
        function main(Algorithm,Problem)                                   
            %% Parameter setting
            num_Samples = 30;                                                              % 真实评估采样次数   
            %% Generate the weight vectors
            [WV,Problem.N] = UniformPoint(Problem.N,Problem.M);                            % 生成100*2的均匀分布的权重向量W
            T =  ceil(Problem.N/10);                                                       % 邻居个数，返回大于等于的最小整数，eg.,10           
            %% Detect the neighbours of each solution
            B = pdist2(WV,WV);                                                             % 计算两组权重之间两两成对的距离
            [~,B] = sort(B,2);                                                             % eg.,排序是升序
            B = B(:,1:T);                                                                  % 提取每个解和其他解间权重距离最小的前10个解所在的索引，100*10          
            %% Generate random population
            Population = Problem.Initialization();                                         % 初始化种群                  
            Population= Greedy_Repair_Nonscalar(Problem, Population);  
            [Population] = Multiple_Sampling(Problem,Population,num_Samples);              % 多次评估种群       
            Z = min(Population.objs,[],1);                                                 % MOKP问题目标是最小值，理想点
            [FrontNo,~] = NDSort(Population.objs,inf);                                     % 对种群进行非支配排序，inf表示对整个总体进行排序。
            NPOP = Population(FrontNo == 1);                                               % 选出初始种群第一层的解作为非支配解集
            gen = 1;
            %% Optimization
            while Algorithm.NotTerminated(NPOP)
                for i = 1 : Problem.N                                                     % 遍历种群中的每个解
                    P = B(i,randperm(size(B,2)));                                         % 将每个解的邻居顺序打乱并提取                           
                    Offspring = OperatorGAhalf(Problem,Population(P(1:2)));               % 交叉变异
                    
                    Offspring = Greedy_Repair_Nonscalar(Problem, Offspring);       % 使用当前子问题的权重修复子代
                    [Offspring] = Multiple_Sampling(Problem,Offspring,num_Samples);       % 多次评估     
                    NPOP = [NPOP, Offspring];                                             % 将后代解并入非支配解集
                    [FrontNo,~] = NDSort(NPOP.objs,inf);                                  % 非支配排序
                    NPOP = NPOP(FrontNo == 1);
                    % 维持外部集大小为N
                    if length(NPOP) > Problem.N
                       CrowdDis = CrowdingDistance(NPOP.objs);
                       [~, idx] = sort(CrowdDis, 'descend');
                       NPOP = NPOP(idx(1:Problem.N));
                    end
                    % 更新外部非支配解集                 
                    Z = min(Z,Offspring.objs);                                            % 更新全局最优解                    
                    g_old = sum(Population(P).objs.*WV(P,:),2);                           % 计算后代加权聚合函数，eg.,g=lamda*f,10*1
                    g_new = sum(Offspring.objs.*WV(P,:),2);                                % 计算父代加权聚合函数
                    Population(P(g_old>=g_new)) = Offspring;                              % 好解替换种群差解
                end
                gen = gen + 1;
            end
        end
    end
end
