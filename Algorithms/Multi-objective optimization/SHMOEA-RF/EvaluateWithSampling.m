function Obj_mean = EvaluateWithSampling(Population, Problem, num_Samples)
    % 原有的单次评估值算作一次采样，故真实采样9次
    % 输入：种群，问题 Problem，多次采样数量 numSamples
    % 输出：种群的目标均值 Obj_mean
    
    Pop_decs = Population.decs;           % 提取种群的所有决策变量
    Obj_first = Population.objs;          % 提取随机采样的单次目标函数值作为第一次采样值   
    num_Individuals = size(Pop_decs, 1);  % 种群的大小100
    num_Objectives = Problem.M;           % 背包问题的目标数量（eg.,2）

    Objs_sampled = zeros(num_Individuals, num_Objectives, num_Samples); % 初始化矩阵，用于存储多次采样的目标值, 100*2*10
    Objs_sampled(:,:,1) = Obj_first;                                    % 单次目标函数值作为第一次采样值  
    
    for i = 2:num_Samples                                  % 开始第2次到第num_Samples次采样
        Objs_sampled(:, :, i) = Problem.CalObj(Pop_decs);  % 使用背包问题中定义的目标函数计算方法
        Problem.FE = Problem.FE + size(Pop_decs,1);        % 累加算法的总评估次数
    end
    
    Obj_mean = mean(Objs_sampled, 3);    % 在采样维度上对所有已采样的目标函数取均值,100*2
end
