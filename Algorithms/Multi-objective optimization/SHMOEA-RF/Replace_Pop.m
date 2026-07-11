function Pop = Replace_Pop(PopulationA,PopulationB,num_Select)

% num_Select   当前进化种群即将被替换的解的数量
% PopulationA  用于替换解的种群 A
% PopulationB  被替换的解的种群 B,含有差解

% 从种群 A 中随机选择 num_Select 个解，用于替换种群 B 中的 num_Select 个解
selected_Indices_A = randperm(size(PopulationA, 2), num_Select);
selected_Solutions_A = PopulationA(:,selected_Indices_A);

% 从种群 B 中随机选择 num_Select 个解
selected_Indices_B = randperm(size(PopulationB, 2), num_Select);

% 将种群 B 中选择的 num_Select 个解替换为种群 A 中已选择的解
PopulationB(:,selected_Indices_B) = selected_Solutions_A;
% 输出替换后的种群 B 
Pop = PopulationB

end