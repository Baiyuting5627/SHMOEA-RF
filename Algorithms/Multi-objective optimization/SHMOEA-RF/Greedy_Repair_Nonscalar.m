function [Population] = Greedy_Repair_Nonscalar(obj,Population) 
  %% 提取决策变量及函数值
             PopDec = Population.decs;
             PopObj = Population.objs;  
             PopCon = Population.cons;
             
  %% 启发式约束修正
             % 通过逐步移除性价比最低的物品（利润除以重量的值较小），优先保留高性价比（单位重量利润较大）物品，修正种群中不满足重量约束的解，尽量减少对解质量的影响。
             C = sum(obj.W,2)/2;                         % 得到2*1的矩阵
             
             [~,rank] = sort(max(obj.P./obj.W));         
             % 从两行矩阵中找到单位重量利润的最大值，并对物品进行升序排序。
             % sort函数返回排序后的值及排序后的索引rank，这些索引表示排序后物品的位置。
             % 排序后的所在的索引1*250，如106，138，68，表示单位重量利润最小的物品是第106，然后是第138
             for i = 1 : size(PopDec,1)                  % 遍历所有决策变量
                 current_P =  obj.N*obj.D - PopObj(i,:); % 当前利润向量,1*2
                 current_W = PopCon(i,:) + C';           % 当前重量向量,1*2               
                 while any(current_W > C')               % 返回1*2的二维逻辑数组，1 1 或0 1,检查给定的逻辑矩阵中是否至少有一个元素为 true(1)
                       k = find(PopDec(i,rank),1);  
                       %在排序后的物品列表 rank 中，寻找当前解 PopDec(i,:) 中第一个值为 1的物品索引
                       % PopDec(i,rank)，返回单位重量利润最小的物品对应的决策变量，如第106个决策为0，第138为1...
                       % find(X,1),返回输入数组X中,第一个非零元素的索引。                      
                       % 检查 k是否为空数组。若为空，跳出循环
                       if isempty(k)                
                          break;                                    
                       end           
                       item_idx = rank(k);                           % 被选择移除的物品编号            
                       PopDec(i, item_idx) = 0;                      % 移除物品后的决策向量                      
                       current_P = current_P - obj.P(:, item_idx)';  % 更新移除该物品后的总价值
                       current_W = current_W - obj.W(:, item_idx)';  % 更新移除该物品后的总重量
                 end                      
                 PopObj(i,:) =  obj.N*obj.D - current_P;             % 更新最终修正后的目标值                          
                 PopCon(i,:) = current_W-C';                         % 更新最终修正后的约束值 
             end
             
  %% 结构解
       for i=1:size(Population,2)           
           Population(i) = SOLUTION(PopDec(i,:),PopObj(i,:), PopCon(i,:));
       end
end
