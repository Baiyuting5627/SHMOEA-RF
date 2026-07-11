function [ trees,S,nc] = GenerateRFa( X,Y,N)
c=size(X,2);
n=size(X,1);
nc=round(c.^0.5)*ones(N,1);
trees=cell(N,1);
S=[];
for i=1:N
    t=randperm(c);
    S=[S;t];
end
for i=1:N
    I=randi([1,n],[n,1]);
    tree=RegressionTree.fit(X(I,S(i,1:nc(i))),Y(I));
    trees{i,1}=tree;
end


end

