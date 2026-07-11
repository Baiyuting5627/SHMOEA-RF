classdef MOEA_RF < ALGORITHM
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
            num_selet_pop_first = 50;     % 从当前种群中初次选择用于鲁棒平均多次采样评估并更新外部集的(单次评估)解的数量
            num_selet_pop_second = 25;     % 从初次选择的解中，根据RF预测后的函数值进一步选择用于鲁棒平均多次采样评估并更新外部集的(单次评估)解的数量          

            %% Generate the weight vectors
            [WV,Problem.N] = UniformPoint(Problem.N,Problem.M);                % 生成100*2的均匀分布的权重向量W
            T =  ceil(Problem.N/10);                                           % 邻居个数，返回大于等于的最小整数，eg.,10
            %% Detect the neighbours of each solution
            B = pdist2(WV,WV);                                                 % 计算两组权重之间两两成对的距离
            [~,B] = sort(B,2);                                                 % eg.,排序是升序
            B = B(:,1:T);                                                      % 提取每个解和其他解间权重距离最小的前10个解所在的索引，100*10

            %% 初始化
            Population = Problem.Initialization();                                    % 在Pij是随机的背包问题下，所有个体均为单次采样评估           

            [Population] = Greedy_Repair_Nonscalar(Problem, Population);             % 随机初始种群进行启发式修正(单次评估的函数值
            Z = min(Population.objs, [], 1);                                          % platemo下MOKP问题的目标是最小值，Z是理想点  
            [~, FrontNo, CrowdDis] = EnvironmentalSelection(Population, Problem.N);   % 基于目标值进行排序的环境选择，拥挤距离越大表示解更分散，有助于维持解的多样性  
            NPOP = Population(FrontNo == 1);      % 将位于第一层的非支配解放入外部集中，外部集存储的始终是非支配解
            gen = 1;                              % 初始化种群的迭代进化次数为1    
            Trees = [];  S = [];  NCtree = [];    % 初始化随机森林参数    
            A = NPOP;                             % MOEAD的外部集；     
             
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
                    [~, FrontNo, CrowdDis] = EnvironmentalSelection(Population, Problem.N);
                    MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);                             % 锦标赛选择形成交配池，如果层数相同，则优先选最大拥挤距离。eg.,返回1*100的索引
                    Offspring = OperatorGA(Problem, Population(MatingPool));                                        %  交配池中的候选种群进行交叉变异，生成后代解
                    [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);
                    [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);   % 父本和后代环境选择出下一代进化种群                               
                else
                    %%  MOEA/D 子问题分解方法生成新解，eg.,gen = 20-39,60-79,100-119,...            
                    for i = 1 : Problem.N                                                      % 遍历种群中的每个解
                        P = B(i,randperm(size(B,2)));                                          % 将每个解的邻居顺序打乱并提取
                        Offspring = OperatorGAhalf(Problem,Population(P(1:2)));                % 交叉变异
                        [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);        % 启发修正
                        A = [A, Offspring];                                               % 将后代解并入非支配解集
                        [FrontNo_M,~] = NDSort(A.objs,inf);                                 % 非支配排序
                        A = A(FrontNo_M == 1);                                              % 更新非支配解集
                        % 维持外部集大小为N
                        if length(A) > Problem.N
                           CrowdDis = CrowdingDistance(A.objs);
                           [~, idx] = sort(CrowdDis, 'descend');
                           A = A(idx(1:Problem.N));
                        end
                        Z = min(Z,Offspring.objs);                                        % 更新全局最优解 
                        g_old = sum(Population(P).objs.*WV(P,:),2);                       % 加权聚合函数，eg.,g=lamda*f,10*1
                        g_new = sum(Offspring.objs.*WV(P,:),2); 
                        Population(P(g_old>=g_new)) = Offspring;                          % 好解替换种群差解
                    end
                    
                end                                               % 两种生成新解的方式执行完毕                
                %% 筛选好解鲁棒评估
                if mod(gen, num_beta) == 0                                                                          % 当迭代进化次数gen为20的倍数时进行初步选解             
                   %% 初步筛选阶段
                    if mod(gen, 2 * num_beta) < num_beta                                                            % NSGA-II 初步选解                                   
                       current_pop = Population;                                                                    % 当前种群为环境选择后的解
                       Selected_Pop_first = Select_Solution(current_pop, FrontNo, CrowdDis, num_selet_pop_first);   % 按照层级和拥挤距离优先选择，已计算好层级和拥挤距离          
                    else                                                                              % MOEA/D 初步选解
                       current_pop =  A;                                                            % 当前种群为权重向量保存的解
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