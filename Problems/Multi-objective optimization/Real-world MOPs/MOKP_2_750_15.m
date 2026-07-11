classdef MOKP_2_750_15< PROBLEM
% <1999> <multi/many> <binary> <large/none> <constrained>
% The multi-objective knapsack problem
% 无修正，修正在repair中
% 10%的变量进行扰动
%------------------------------- Reference --------------------------------
% E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
% comparative case study and the strength Pareto approach, IEEE
% Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        P;	% Profit of each item according to each knapsack
        W;  % Weight of each item according to each knapsack
    end
    methods
        %% Default settings of the problem
        function Setting(obj)   % obj: PROBLEM class
            % Parameter setting
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 750; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = 4 + zeros(1,obj.D);%离散
            % Randomly generate profits and weights
            file = sprintf('MOKP-M%d-D%d_15.mat',obj.M,obj.D);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'P','W');
            else
                P = randi([10,100],obj.M,obj.D);  %Profit
                W = randi([10,100],obj.M,obj.D);  %Weight
                save(file,'P','W');
            end
            obj.P = P;
            obj.W = W;
        end            
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)                           
             normal_mu = 2;   normal_sigma = 0.5;                          % Profit参数P_ij正态分布均值与方差
             ratio_randprofit = 0.15;                                       % 对应Item的Profit随机扰动比例（前10%的变量进行扰动）
             P_bak = obj.P;                                                % 原始价值向量
             [numRows, numCols] = size(obj.P);                             % 获取矩阵维度            
             num_Perturb = round(ratio_randprofit * numCols);                          % 前百分之a的变量数,确保是整数
             Idx = 1:num_Perturb;                                                      % 选择前百分之a的变量索引,1*25,1...25
             profit_perturb = lognrnd(normal_mu, normal_sigma, numRows, num_Perturb);  %对数正态分布扰动,1*25
             P_bak(:,Idx) = P_bak(:,Idx) + profit_perturb;                 %生成扰动后的价值向量  
             P_bak = min(max(P_bak, 10), 100);                             % 限定扰动后的值还在 [10, 100] 区间                         
             total_value = PopDec* P_bak';                                 % 计算总价值 (100×2)
             PopObj = repmat(100*obj.D, size(PopDec, 1), 2)-total_value;     %100*2 将目标值从最大值基准中扣除，转化为一个极小化问题
        end
        %% Calculate constraint violations 输入的是已经修正且离散的变量
        function PopCon = CalCon(obj,PopDec)
            PopCon = PopDec*obj.W' - repmat(sum(obj.W,2)'/2,size(PopDec,1),1);
        end

       %% Generate points on the Pareto front
       function R = GetOptimum(obj,~)
           %data=load('D:\白雨婷\PlatEMO-master-bai\PlatEMO\Real_PF\PF_2_750_15.mat');% 2-250
           data=load('D:\白雨婷\PlatEMO_latest\PlatEMO-master\PlatEMO\Real_PF\Re_PF_2_750_15.mat');
           R=data.PS.objs;          
           %R=data.result{1,2}.objs;%提取最小化后的PF
           %PF=100*obj.D-data.result{1,2}.objs;%计算真实的最大化PF
           %scatter(PF(:,1),PF(:,2));%画出真实帕累托前沿
       end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(Population.decs*obj.P',{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
        end
    end
end