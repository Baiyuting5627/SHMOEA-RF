function [obj,s ] = Predict_Obj_RF(POP,Trees,S,NCtree,c,m)

n=size(POP,1);

obj=zeros(n,m);

s=zeros(n,m);

for i=1:m 
    
    [ obj(:,i),s(:,i) ] = EnsemblePredicta( Trees(:,i),S(:,:,i),NCtree(:,i),POP(:,1:c) );

end

end

