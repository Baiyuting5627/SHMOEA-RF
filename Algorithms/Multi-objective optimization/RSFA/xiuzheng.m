function offspring=xiuzheng(offspring)
[N,D]=size(offspring);

for i=1:N
    if find(offspring(i,:)==0)
        ind=find(offspring(i,:)==0);
        k=0;
        for j=1:D
            if find(offspring(i,:)==j)
                continue
            else
                k=k+1;
                offspring(i,ind(k))=j;
            end
        end
    end
end

a=1
        




end