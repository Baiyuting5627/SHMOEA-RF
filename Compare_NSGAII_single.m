classdef Compare_NSGAII_single < ALGORITHM
% <2002> <multi> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II 
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
        function main(Algorithm,Problem)         
            num_Samples = 10;                                                      % 真实评估采样次数             
            %% Generate random population
            Population = Problem.Initialization();                                 % 初始化生成种群           
            [Population] = Greedy_Repair_Nonscalar(Problem, Population);           % 启发修正
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);   % 环境选择
            %% Optimization
            gen=1;
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);   % 锦标赛选择              
                Offspring  = OperatorGA(Problem,Population(MatingPool));           % 交叉变异
                [Offspring] = Greedy_Repair_Nonscalar(Problem, Offspring);         % 启发修正           
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);   % 环境选择生成新种群
                gen=gen+1;                                                                                  % 累加进化次数
            end
        end
    end
end