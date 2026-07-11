classdef Compare_NSGAII_RF< ALGORITHM
% <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
% 单独 NSGAII + 单次+多次评估 + RF + 启发式修正
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
        function main(Algorithm, Problem)
            %% 参数设置
            num_x = Problem.D;   % 问题决策空间维数(变量个数)(eg.,50)       
            num_f = Problem.M;   % 问题目标空间维数(函数个数)(eg.,2)
            num_cons = num_f;    % 约束条件的数量(eg., 2, 背包问题中目标函数个数等于约束个数)
            num_rf = 100;        % 随机森林模型中节点的数量
            num_sample = 30;     % 鲁棒多次平均采样评估的样本个数
            num_select_archive = 20;     % 从外部集挑选鲁棒平均多次采样评估解用于更新当前种群数量(每10代执行1次)
            num_selet_pop_first = 50;    % 从当前种群中初次选择用于鲁棒平均多次采样评估并更新外部集的(单次评估)解的数量
            num_selet_pop_second = 25;   % 从初次选择的解中，根据RF预测后的函数值进一步选择用于鲁棒平均多次采样评估并更新外部集的(单次评估)解的数量
            num_beta = 20;               % NSGA-II的算法参数

            %% 初始化
            Population = Problem.Initialization();                                    % 在Pij是随机的背包问题下，所有个体均为单次采样评估           
            [Population] = Greedy_Repair_Nonscalar(Problem, Population);              % 随机初始种群进行启发式修正(单次评估的函数值)
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);% 基于目标计算的obj值和con分层.拥挤距离CrowdDis越大表示解更分散，有助于维持解的多样性  
            NPOP = Population(FrontNo == 1);      % 将位于第一层的非支配解放入外部集中，外部集存储的始终是非支配解
            gen = 1;                              % 初始化种群的迭代进化次数为1          
          %% Optimization
            while Algorithm.NotTerminated(NPOP)
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);                             % 锦标赛选择形成交配池，如果层数相同，则优先选最大拥挤距离。eg.,返回1*100的索引
                Offspring = OperatorGA(Problem, Population(MatingPool));                                        %  交配池中的候选种群进行交叉变异，生成后代解
                [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);                                      % 对后代解进行启发式修正（单次评估的函数值)
                [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);   % 父本和后代环境选择出下一代进化种群            
               %%  满足条件进行多次迭代进化  
                   if mod(gen,num_beta) == 0                                                                  % 当迭代进化次数gen为10的倍数时执行                                                                     
                      Selected_Pop_first = Select_Solution(Population,FrontNo,CrowdDis,num_selet_pop_first);  % 从当前进化种群初步选择出特定数量的解，层级小的优先选择，同一前沿层下拥挤距离大的优先选择                                       
                      if size(NPOP,1) >=100                                      % 当外部集NPOP达到一定数量后，利用RF预测出初次选择解的鲁棒多次平均采样评估值                         
                         num_data_RF = min(num_selet_pop_first, size(NPOP, 1));  % 防止NPOP中元素个数不足50时越界
                         NPOP_Data = NPOP(end - num_data_RF + 1 : end);          % 从NPOP中提取最新放入的50个鲁棒多次平均采样评估解
                         TrainData = [NPOP_Data.decs,NPOP_Data.objs,NPOP_Data.cons];                     % 构成用作训练RF模型的训练集
                         [Trees,S,NCtree] = Surrogate_ObjCon_RF(TrainData,num_x,num_f,num_cons,num_rf);  % 训练集训练RF模型（鲁棒多次评估评估后的解才可以构成训练集）
                         Predict_Pop = [Selected_Pop_first.decs,Selected_Pop_first.objs,Selected_Pop_first.cons];              % 用初次选择的好解构建预测集，eg.,50*254数组
                         [Predict_objs] = Predict_Obj_RF( Predict_Pop,Trees,S,NCtree,num_x,num_f);                             % 用RF预测初次选择的好解的多次平均采样目标值，eg.,50*2数组           
                         [Selected_Pop_second] = Select_Solution_NSGA(Selected_Pop_first,Predict_objs,num_selet_pop_second );  % 从初次选择的解中根据RF预测后的函数值进一步选择出的好解    
                         Selected_Pop_second = Multiple_Sampling(Problem,Selected_Pop_second,num_sample);  % 对第二层选择出的特定好解进行鲁棒多次平均采样评估           
                         TPOP = [NPOP,Selected_Pop_second];                                                % 将鲁棒多次评估后的解加入初次外部集
                      else
                         Selected_Pop_first = Multiple_Sampling(Problem,Selected_Pop_first,num_sample);    % 对从当前进化种群选择出的特定好解进行鲁棒多次平均采样评估 
                         TPOP = [NPOP,Selected_Pop_first];                                                 % 将鲁棒多次评估后的解加入初次外部集
                      end                  
                      Population = Replace_Pop(TPOP,Population,num_select_archive);   % 从外部集任选num_select_archive个解，更新当前进化种群中的任意num_select_archive个解 
                      NPOP = Nondominated( TPOP,num_x,num_f);                         % 对外部集进行非支配排序，确保外部集中始终存储的是鲁棒多次平均采样评估的非支配解           
                   end          % gen不满足条件，则继续进行当前种群的单次迭代进化
                   gen=gen+1;   % 累加种群的进化迭代次数  
            end
        end
   end
end