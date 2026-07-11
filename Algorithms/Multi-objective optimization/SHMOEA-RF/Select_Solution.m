function Selected_Pop = Select_Solution(Population,FrontNo,CrowdDis,num_selet)
% 适用于已经得到分层的排序和拥挤距离后，如何选择个体
% 前沿层数小的优先选择，同一前沿层数下拥挤距离大的优先选择
% 根据前沿层数进行升序排序，得到对应的索引
[~, sort_Index1] = sort(FrontNo); 

% 首先将 CrowdDis 数组按照 sortIndex1 重新排列，确保 CrowdDis 的顺序与 FrontNo 的排序顺序一致
% 在此基础上，对重新排列后的CrowdDis进行降序排序，并返回排序后对应的索引 sortIndex2。
[~, sort_Index2] = sort(CrowdDis(sort_Index1), 'descend'); 

% 将 sortIndex2 用于 sortIndex1，确保先按前沿编号升序排序，同一前沿编号内按拥挤距离降序排序。
final_SortIndex = sort_Index1(sort_Index2);

% 选择要选择的前 num_selet 个解的索引
selected_Indices = final_SortIndex(1:num_selet);

% 提取要选择的前 num_selet 个解
Selected_Pop = Population(:,selected_Indices);
end
