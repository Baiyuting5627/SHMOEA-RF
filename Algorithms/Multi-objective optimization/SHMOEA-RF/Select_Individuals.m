function Selected_Pop = Select_Individuals(Population, N)
    % 适用于MOEAD中,需重新计算分层和拥挤距离时，针对一个种群如何选择出好解
   
     Selected_Pop = [];                                                 % 初始化选择的个体集合    
    [FrontNo, ~] = NDSort(Population.objs, size(Population.objs, 1));   % 对输入种群进行非支配排序，得到1*100的数组   
        
    for i = 1:max(FrontNo)                                              % 遍历位于每一层的个体，eg.,i=1:17         
        current_front_indices = find(FrontNo == i);                     % 获取位于当前前沿的所有个体索引             
        if length(current_front_indices) >= N - length(Selected_Pop)    % 如果当前前沿的个体数量超出剩余需选择的数量，根据拥挤距离进一步选择           
            CrowdDis = CrowdingDistance(Population(current_front_indices).objs, FrontNo(current_front_indices));              % 计算拥挤距离                    
            [~, sorted_indices] = sort(CrowdDis, 'descend');                                                                  % 按拥挤距离降序排序            
            Selected_Pop = [Selected_Pop, Population(current_front_indices(sorted_indices(1:(N - length(Selected_Pop)))))];   % 选择剩余需要的个体            
            break;        
        else           
            Selected_Pop = [Selected_Pop,Population(current_front_indices)];                                                  % 否则，将当前前沿的所有个体加入选择集合        
        end        
    end   
end
