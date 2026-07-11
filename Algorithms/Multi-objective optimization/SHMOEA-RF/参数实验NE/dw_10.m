classdef dw_10 < ALGORITHM
% <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II 
% 改进NSGAII + 动态MOEAD 间隔20代协同进化+ RF
% weights 修正

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
%--------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            %% 参数设置
            num_x = Problem.D;   % 问题决策空间维数(变量个数)(eg.,50)       
            num_f = Problem.M;   % 问题目标空间维数(函数个数)(eg.,2)
            num_cons = num_f;    % 约束条件的数量(eg., 2, 背包问题中目标函数个数等于约束个数)
            num_beta = 20;       % NSGA-II 与 MOEA/D 切换优化代数
            num_rf = 100;        % 随机森林模型中节点的数量
            num_sample = 30;     % 鲁棒多次平均采样评估的样本个数
            num_select_archive = 20;       % 从外部集挑选鲁棒平均多次采样评估解用于更新当前种群数量(每20代执行1次)
            num_selet_pop_first = 50;      % 从当前种群中初次选择用于鲁棒平均多次采样评估并更新外部集的(单次评估)解的数量
            num_selet_pop_second = 25;     % 从初次选择的解中，根据RF预测后的函数值进一步选择用于鲁棒平均多次采样评估并更新外部集的(单次评估)解的数量          
            num_weight_set = 100;          % 加权分解子问题候选(均匀生成)权重向量总量
            num_weights = 10;               % MOEA/D 每一代演化所需权重数量（包括每一个单目标子问题对应权重向量)
            num_pop_elite = 20;            % 单个加权子问题关联最好解的数量              
            
            %% 初始化
            Population = Problem.Initialization();                                    % 在Pij是随机的背包问题下，所有个体均为单次采样评估           

            [Population] = Greedy_Repair_Nonscalar(Problem, Population);             % 随机初始种群进行启发式修正(单次评估的函数值
            Z = min(Population.objs, [], 1);                                          % platemo下MOKP问题的目标是最小值，Z是理想点  
            [~, FrontNo, CrowdDis] = EnvironmentalSelection(Population, Problem.N);   % 基于目标值进行排序的环境选择，拥挤距离越大表示解更分散，有助于维持解的多样性  
            NPOP = Population(FrontNo == 1);      % 将位于第一层的非支配解放入外部集中，外部集存储的始终是非支配解
            gen = 1;                              % 初始化种群的迭代进化次数为1    
            Trees = [];  S = [];  NCtree = [];    % 初始化随机森林参数    
            POP = [];                             % 初始化存储权重相关联好解的外部集
            
            %% 优化阶段
            while Algorithm.NotTerminated(NPOP)            
               %% 训练随机森林模型               
                if size(NPOP, 1) >= 100                                          % 当外部集NPOP达到一定数量后，利用RF预测出初次选择解的鲁棒多次平均采样评估值 
                    num_data_RF = min(num_selet_pop_first, size(NPOP, 1));       % 防止NPOP中元素个数不足50时越界
                    NPOP_Data = NPOP(end - num_data_RF + 1:end);                 % 从NPOP中提取最新放入的50个鲁棒多次平均采样评估解
                    TrainData = [NPOP_Data.decs, NPOP_Data.objs, NPOP_Data.cons];                           % 构成用作训练RF模型的训练集
                    [Trees, S, NCtree] = Surrogate_ObjCon_RF(TrainData, num_x, num_f, num_cons, num_rf);    % 训练集训练RF模型（鲁棒多次评估评估后的解才可以构成训练集）
                end                 
               %% 两种方式生成后代新解
                if mod(gen, 2 * num_beta) < num_beta
                   %% NSGA-II 分层排序方法生成新解， eg.,gen = 1-19,40-59,80-99...   
                    MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);                             % 锦标赛选择形成交配池，如果层数相同，则优先选最大拥挤距离。eg.,返回1*100的索引
                    Offspring = OperatorGA(Problem, Population(MatingPool));                                        %  交配池中的候选种群进行交叉变异，生成后代解
                    %[Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);                                     % 对后代解进行启发式修正（单次评估的函数值)
                    [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);
                    [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);   % 父本和后代环境选择出下一代进化种群                               
                else
                   %%  MOEA/D 子问题分解方法生成新解，eg.,gen = 20-39,60-79,100-119,...            
                    Weights = Generate_Dynamic_Weights(num_weight_set, num_weights, Problem.M);   % 每一个单目标子问题根据最远距离生成5个权重向量                                                                                 
                    for i = 1:size(Weights, 1)                                                    % 遍历子问题下的每一个权重    
                        weight = Weights(i, :);                                                   % 获取当前子问题对应的第i个权重值 
                        aggregate_function = sum(Population.objs .* weight, 2);     % 根据加权和计算当前种群在该权重下的聚合函数(eg.,100*1 的向量）
                        [~, indices] = sort(aggregate_function);                    % 将种群加权聚合函数按照从小到大的顺序进行排列
                        selected_indices = indices(1:num_pop_elite);                % 根据排序，提取前num_pop_elite个最小聚合函数值对应的索引
                        candidate_pop = Population(selected_indices);               % 提取种群中在当前权重下，具有最小聚合函数值的num_pop_elite个候选解 
                        Offspring = OperatorGA(Problem, candidate_pop);             % 候选解进行交叉变异,生成后代个体
                        Offspring = Greedy_Repair_Nonscalar(Problem, Offspring);    % 后代个体基于当前权重进行启发式修正(单次评估的函数值)                                                     
                        %Offspring = Greedy_Repair_Scalar(weight, Problem, Offspring); 
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
                end                                               % 两种生成新解的方式执行完毕                
                %% 筛选好解鲁棒评估
                if mod(gen, num_beta) == 0                                                                          % 当迭代进化次数gen为20的倍数时进行初步选解             
                   %% 初步筛选阶段
                    if mod(gen, 2 * num_beta) < num_beta                                                            % NSGA-II 初步选解                                   
                       current_pop = Population;                                                                    % 当前种群为环境选择后的解
                       Selected_Pop_first = Select_Solution(current_pop, FrontNo, CrowdDis, num_selet_pop_first);   % 按照层级和拥挤距离优先选择，已计算好层级和拥挤距离          
                    else                                                                              % MOEA/D 初步选解
                       current_pop =  POP;                                                            % 当前种群为权重向量保存的解
                       Selected_Pop_first = Select_Individuals(current_pop , num_selet_pop_first);    % 按照层级和拥挤距离优先选择，计算分层排序保留第一层的解，若数量不足再计算拥挤距离选择剩余解
                    end                                                                               % 两种初步筛选解的方式执行完毕                    
                   %% RF二次筛选阶段          
                    if ~isempty(Trees)                                                                               % 如果随机森林模型已训练，使用RF预测目标值并进一步筛选解
                        Predict_Pop = [Selected_Pop_first.decs, Selected_Pop_first.objs, Selected_Pop_first.cons];   % 用初次选择的好解构建预测集，eg.,50*254数组
                        Predict_objs = Predict_Obj_RF(Predict_Pop, Trees, S, NCtree, num_x, num_f);                  % 用RF预测初次选择的好解的多次平均采样目标值，eg.,50*2数组  
                        Selected_Pop_second = Select_Solution_NSGA(Selected_Pop_first, Predict_objs, num_selet_pop_second);    % 从初次选择的解中根据RF预测后的函数值进一步选择出的好解    
                        Selected_Pop_second = Multiple_Sampling(Problem, Selected_Pop_second, num_sample);                     % 对第二层选择出的特定好解进行鲁棒多次平均采样评估     
                        TPOP = [NPOP, Selected_Pop_second];                                               % 将鲁棒多次评估后的解加入初次外部集TPOP
                    else                                                                                  % 外部集数量不足，模型没有训练好的情况下
                        Selected_Pop_first = Multiple_Sampling(Problem,Selected_Pop_first,num_sample);    % 对从当前进化种群选择出的特定好解进行鲁棒多次平均采样评估 
                        TPOP = [NPOP, Selected_Pop_first];                                                % 将鲁棒多次评估后的解加入初次外部集TPOP
                    end

                    Population = Replace_Pop(TPOP, Population, num_select_archive);   % 从外部集任选num_select_archive个解，更新当前进化种群中的任意num_select_archive个解 
                    [FrontNo_moead,~] = NDSort(TPOP.objs,inf);                        % 对外部集进行非支配排序，确保外部集中始终存储的是鲁棒多次平均采样评估的非支配解   
                    NPOP = TPOP(FrontNo_moead == 1);           % 将第一层的非支配解放入外部集中
                    Z = min([Z;NPOP.objs],[],1);               % 更新全局最优解 
                    if length(NPOP) > Problem.N
                       CrowdDis_np = CrowdingDistance(NPOP.objs);            % 计算拥挤距离
                       [~, idx_np] = sort(CrowdDis_np, 'descend');              % 从大到小排序
                       NPOP = NPOP(idx_np(1:Problem.N));                     % 保留拥挤距离最大的 Problem.N 个解
                    end
                end                                            
                gen = gen + 1;                                 % 累加种群的进化迭代次数
            end
        end
    end
end