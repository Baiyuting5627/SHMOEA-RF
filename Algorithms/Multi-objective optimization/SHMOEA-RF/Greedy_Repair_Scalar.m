function [Population] = Greedy_Repair_Scalar(WV, obj, Population)
    % 参数设置与扰动生成（与CalObj保持一致）
    PopDec = Population.decs;%
    normal_mu = 2; normal_sigma = 0.5;
    ratio_randprofit = 0.1;
    P_bak = obj.P;
    [numRows, numCols] = size(obj.P);               % 2*250
    num_Perturb = round(ratio_randprofit * numCols);
    Idx = 1:num_Perturb;
    profit_perturb = lognrnd(normal_mu, normal_sigma, numRows, num_Perturb);
    P_bak(:,Idx) = P_bak(:,Idx) + profit_perturb;  % 扰动后价值矩阵,2*250
    min_threshold = 1e-10;

    C = sum(obj.W,2)/2;  % 容量向量
    
    for i = 1:size(PopDec,1)
        x = PopDec(i,:);         %1*250
        % 计算当前总价值（按目标）
        total_val = x * P_bak';  % 1*2 向量
        current_U = log(max(total_val, min_threshold));  % 1x2 对数效用
        current_W = x * obj.W';  % 1x2 当前重量
        % 当任一约束违反时循环
        while any(current_W > C')
            % 当前聚合值
            aggVal = WV .* current_U';  % WV若只有一个向量，即为1*2的标量，(1*2)*(2*1）=1*1,一个数值
            % 找到选中物品
            indicesOne = find(x == 1);
            exceeded = find(current_W > C');  % 超重目标索引，1或2
            minImpact = inf;
            bestIdx = -1;
            
            for j = indicesOne
                % 临时移除物品j后的总价值
                %tempDec=x;
                %tempDec(j)=0;
                new_total = total_val - P_bak(:,j)';  % 1*2
                new_U = log(max(new_total, min_threshold));% 1*2
                new_agg = WV .* new_U';% (1*2)*(2*1)=1*1
                % 影响计算（分母为物品j在超重背包上的重量和）
                impact = (aggVal - new_agg) / sum(obj.W(exceeded, j));
                if impact < minImpact
                    minImpact = impact;
                    bestIdx = j;
                end
            end
            if bestIdx ~= -1
                x(bestIdx) = 0;
                total_val = total_val - P_bak(:,bestIdx)';
                current_U = log(max(total_val, min_threshold));% 1*2;
                current_W = current_W - obj.W(:,bestIdx)';%1*2
            else
                break;
            end
        end
        % 更新种群
        PopDec(i,:) = x;
        PopObj(i,:) = obj.N*obj.D - current_U;  % 转为最小化形式
        PopCon(i,:) = current_W - C';
    end
    % 重构SOLUTION对象
    for m = 1:size(Population,2)
        Population(m) = SOLUTION(PopDec(m,:), PopObj(m,:), PopCon(m,:));
    end
end