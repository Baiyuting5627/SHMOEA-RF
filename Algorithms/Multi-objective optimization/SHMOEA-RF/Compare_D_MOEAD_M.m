classdef Compare_D_MOEAD_M < ALGORITHM
    % <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
    % 一半标准MOEAD,一半动态权重MOEAD，全部多次评估下
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            num_Samples = 30;                                                  % 真实评估采样次数
            num_weight_set = 100;    % 加权分解子问题候选(均匀生成)权重向量总量
            num_weights = 5;         % MOEA/D 每一代演化所需权重数量（包括每一个单目标子问题对应权重向量)
            num_pop_elite = 20;      % 单个加权子问题关联最好解的数量     
            [WV,Problem.N] = UniformPoint(Problem.N,Problem.M);                            % 生成100*2的均匀分布的权重向量W
            T =  ceil(Problem.N/10);                                                       % 邻居个数，返回大于等于的最小整数，eg.,10           
            %% Detect the neighbours of each solution
            B = pdist2(WV,WV);                                                             % 计算两组权重之间两两成对的距离
            [~,B] = sort(B,2);                                                             % eg.,排序是升序
            B = B(:,1:T);                                                                  % 提取每个解和其他解间权重距离最小的前10个解所在的索引，100*10          

            %% Generate random population
            Population = Problem.Initialization();                                            % 初始化种
            [Population] = Greedy_Repair_Nonscalar(Problem, Population);
            Z = min(Population.objs,[],1);                                                    % MOKP问题目标是最小值，理想点
            [FrontNo,~] = NDSort(Population.objs,inf);                                        % 对种群进行非支配排序，inf表示对整个总体进行排序。
            POP = Population(FrontNo == 1);                                                   % 选出初始种群第一层的解作为非支配解集
            POP = Multiple_Sampling(Problem,POP,num_Samples);                                 % 外部集多次评估;

            
            while Algorithm.NotTerminated(POP)
                % ---- 根据评估次数是否过半切换搜索方式 ----
                if Problem.FE < Problem.maxFE / 2          % 前半段：标准 MOEA/D
                    for i = 1 : Problem.N                                                     % 遍历种群中的每个解
                        P = B(i,randperm(size(B,2)));                                         % 将每个解的邻居顺序打乱并提取                           
                        Offspring = OperatorGAhalf(Problem,Population(P(1:2)));               % 交叉变异
                        
                        Offspring = Greedy_Repair_Nonscalar(Problem, Offspring);       % 使用当前子问题的权重修复子代
                        [Offspring] = Multiple_Sampling(Problem,Offspring,num_Samples);       % 多次评估     
                        POP = [POP, Offspring];                                               % 将后代解并入外部集
                        [FrontNo,~] = NDSort(POP.objs,inf);                                   % 非支配排序
                        POP = POP(FrontNo == 1);
                        % 维持外部集大小为N
                        if length(POP) > Problem.N
                           CrowdDis = CrowdingDistance(POP.objs);
                           [~, idx] = sort(CrowdDis, 'descend');
                           POP = POP(idx(1:Problem.N));
                        end
                        % 更新全局最优解                 
                        Z = min(Z,Offspring.objs);                                            % 更新全局最优解                    
                        g_old = sum(Population(P).objs.*WV(P,:),2);                           % 计算后代加权聚合函数，eg.,g=lamda*f,10*1
                        g_new = sum(Offspring.objs.*WV(P,:),2);                                % 计算父代加权聚合函数
                        Population(P(g_old>=g_new)) = Offspring;                              % 好解替换种群差解
                    end
                else                                           % 后半段：动态权重搜索
                    Weights = Generate_Dynamic_Weights(num_weight_set, num_weights, Problem.M);   % 每一个单目标子问题根据最远距离生成5个权重向量                                                                                 
                    for i = 1:size(Weights, 1)                                                    % 遍历子问题下的每一个权重    
                        weight = Weights(i, :);                                                   % 获取当前子问题对应的第i个权重值 
                        aggregate_function = sum(Population.objs .* weight, 2);     % 根据加权和计算当前种群在该权重下的聚合函数(eg.,100*1 的向量）
                        [~, indices] = sort(aggregate_function);                    % 将种群加权聚合函数按照从小到大的顺序进行排列
                        selected_indices = indices(1:num_pop_elite);                % 根据排序，提取前num_pop_elite个最小聚合函数值对应的索引
                        candidate_pop = Population(selected_indices);               % 提取种群中在当前权重下，具有最小聚合函数值的num_pop_elite个候选解 
                        Offspring = OperatorGA(Problem, candidate_pop);             % 候选解进行交叉变异,生成后代个体
                        [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);
                        Offspring = Multiple_Sampling(Problem, Offspring,num_Samples);    
                        g_old = aggregate_function(selected_indices);                     % 提取num_pop_elite个候选解的聚合函数值
                        g_new = sum(Offspring.objs.*weight,2);                            % 计算后代个体在当前权重下的加权聚合函数值
                        for k = 1:size(g_new,1)                  % 遍历每一个后代解k，以便对种群进行更新
                            value = g_new(k);                    % 提取当前第k个后代解的加权聚合函数值                             
                            idx = find(g_old >= value);          % 找到候选解的聚合函数值中不小于第k个后代解聚合函数值的索引集
                            if ~isempty(idx)                     % 如果满足条件的索引集非空 （认为第k个后代解优于至少一个候选解）                                               
                                idx = idx(randperm(length(idx)));                          % 随机打乱索引集的排列顺序                                
                                idx = idx(1:min(5, length(idx)));                          % 令索引集中的索引数量不超过5                             
                                for n = 1:length(idx)                                      % 遍历位于索引集中的每一个候选解
                                    Population(selected_indices(idx(n))) = Offspring(k);   % 用第k个后代解更新当前进化种群的候选解
                                end  
                            end
                        end
                        POP = [POP, Offspring];                   % 将修正后单次评估的20个后代个体放入POP中
                    end                                           % 遍历所有权重，每进化1代存储100个解到POP                  
                    POP = Select_Individuals(POP, Problem.N);     % POP按照分层排序保留解，保证始终存储种群大小数量的好解，用于后续初步选择好解       
                end
            end                                          
        end
    end
end