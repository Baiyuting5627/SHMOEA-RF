classdef Compare_NSGAII_single < ALGORITHM
% <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
% 单独 NSGAII + 全部单次评估 + 启发修正,最后一次多次评估


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
            num_sample = 30;  % 最终评估的采样次数
            flag = 0;
            %% 种群初始化（单次评估）
            Population = Problem.Initialization();
            [Population] = Greedy_Repair_Nonscalar(Problem, Population);
            [~, FrontNo, CrowdDis] = EnvironmentalSelection(Population, Problem.N);
            
            %% 优化过程（全部使用单次评估）
            gen = 1;
            while Algorithm.NotTerminated(Population)
                % 交配池
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                
                % 生成后代（单次评估）
                Offspring = OperatorGA(Problem, Population(MatingPool));
                [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);
                
                % 环境选择生成新种群
                [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);            
                gen = gen + 1;
                 %% 输出最终种群（使用多次评估结果）

                % if flag == 0 && Problem.FE >= Problem.maxFE - num_sample*length(Population)
                %    Population = Multiple_Sampling(Problem, Population, num_sample);                                     
                %    flag = 1;
                % end
            end
        end
    end
end