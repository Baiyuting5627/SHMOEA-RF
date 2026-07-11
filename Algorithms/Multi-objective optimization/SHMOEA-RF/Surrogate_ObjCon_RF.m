function [ Trees,S,NCtree ] = Surrogate_ObjCon_RF(TrainData,c,m,nc,ntree )
Trees = [];                   % 初始化
S = zeros(ntree,c,m+nc);      % 生成m+nc 3个100*50的零矩阵,eg.,（：，：，1） （：，：，2）（：，：，3）
NCtree = zeros(ntree,m+nc);   % 生成100*3的零矩阵，和树的数量相关
for i=1:m                     % 遍历每一个目标
     [trees,S(:,:,i),NCtree(:,i)]=GenerateRFa( TrainData(:,1:c),TrainData(:,c+i),ntree); % 等式左边是输出，将每个目标选取的特征数量都输出
     Trees=[Trees,trees];
end

%for i = m+1:m+nc % 对于约束函数
%    [trees,S(:,:,i),NCtree(:,i)] = GenerateRFa( TrainData(:,1:c),TrainData(:,c+i),ntree);
%    Trees=[Trees,trees];
%end

end

