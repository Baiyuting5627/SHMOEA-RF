function [Population] = Multiple_Sampling(obj,Population,numSamples) 
            
            % 提取待采样评估种群的决策变量      
            PopDec = Population.decs; 
                                
            %计算离散变量多次采样评估的平均目标函数值obj,100*2
            mean_pop_obj = EvaluateWithSampling(Population, obj, numSamples); 
            
            % 决策变量结构化
            for i = 1:size(Population,2)       
                
                Population(i) = SOLUTION(PopDec(i,:),mean_pop_obj(i,:),Population(i).con);
            
            end  
            
end
