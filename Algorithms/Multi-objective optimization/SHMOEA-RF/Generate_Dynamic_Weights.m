function weights = Generate_Dynamic_Weights(samples_weight,num_weights,num_dimensions)
    % input 
    % num_dimensions - 权重向量的维度
    % num_weights - 需要生成的权重数量（包括端点)
    % samples_weight - 权重向量库的权重个数 100
    % output
    % weights - 生成的权重矩阵，大小为 num_weights x 2

    %% 初始化权重集合
    C = UniformPoint(samples_weight, num_dimensions);    % 生成100*2个均匀分布的权重向量
    weights = zeros(num_weights,  num_dimensions);       % 存储最终选择的权重

    %% 选择端点权重
    weights(1, :) = [1,0];     % 提取出右端点对应的权重
    weights(2, :) = [0, 1];    % 提取出左端点对应的权重
    selected_count = 2;        % 记录已选择的权重数量
    C(1,:)=[];                 % 去掉左端点值，避免随机选中
    C(end,:) =[];              % 去掉右端点值，避免随机选中
    
    %% 第3个权重从C中随机选择
    if selected_count + 1 <= num_weights       % 如果已选择的权重数量没达到需要生成的权重数量
        rand_index = randi(size(C, 1));        % 从C中随机选择一个索引
        weights(selected_count + 1, :) = C(rand_index, :);    % 提取出对应的权重
        C(rand_index, :) = [];                                % 从C中移除已选择的权重
        selected_count = selected_count + 1;                  % 更新已选择的权重数量
    end

    %% 从第4个权重开始，按最远距离选择
    while selected_count < num_weights           % 如果已选择的权重数量没达到需要生成的权重数量    
        min_distances = zeros(size(C, 1), 1);    % 初始化最小距离
        for i = 1:size(C, 1)                                                         % 对于其他所有权重
            min_distances(i) = min(pdist2(C(i, :), weights(1:selected_count, :)));   % 计算C中每个权重与已选择权重的最小距离
        end
        
        [~, max_index] = max(min_distances);                  % 找到距离最远的权重索引
        weights(selected_count + 1, :) = C(max_index, :);     % 提取出距离最远的权重
        C(max_index, :) = [];                                 % 从C中移除已选择的权重
        selected_count = selected_count + 1;                  % 更新已选择的权重数量
    end
end

