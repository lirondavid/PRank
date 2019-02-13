function [prob,k] = merge(prob1,N1,prob2,N2,index1,index2)
m=zeros(N1*N2,3);
y=1;
for i=1:N1
    for j=1:N2
        m(y,1)=prob1(i)*prob2(j);
        m(y,2)=i;
        m(y,3)=j;
        y=y+1;
    end
end
sortedM=sortrows(m,1);
sortedM=flip(sortedM);
k=1;
for i=1:size(sortedM,1)
    if sortedM(i,2)==index1 && sortedM(i,3)==index2
        k=i;
        break
    end
end
prob=sortedM(1:size(sortedM,1),1);
end
