function sum = upper_rank(alpha,a,d,p,n)

i=1;
sum=0;
j=1;
pro=1;
while(j<=d)
    pro=pro*a(j)/(n(j)^(alpha(j)));
    j=j+1;
end
pro=pro/p;
    
N=1;
j=1;
while(j<=d)
    N=N*n(j);
    j=j+1;
end
    
while(i<=d)
    A=pro^(1/alpha(i));
    A=A*N;
    j=1;
    r=1;
    while(j<=d)
        if j~=i
           r=r*(alpha(i))/(alpha(i)-alpha(j));
        end
        j=j+1;
    end   
    sum=sum+r*A;
    i=i+1;
end
end


    