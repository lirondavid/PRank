function [as,alphas] = upper_estimation_fast(prob,N,key,p)
D=10^10000;
as=0;
alphas=0;
i=1;
foundi=0;
counter=0;

while i<=N && foundi==0
   while(i+1<=N && prob(i)==prob(i+1))
        i=i+1;
   end 
   j=i+1;
   foundj=0;
   while j<=N  && foundj==0
        alpha=(log2(prob(i))-log2(prob(j))) / (log2(j)-log2(i));
        a=prob(j)*j^alpha; 
        k=j;
        if k==i || k==j
            diff=0;
        else                
            diff=-prob(k)+a/(k^alpha);
        end
            
        while( k<N )
            k=k+1;
            if k~=i && k<=N && k~=j
                diff=-prob(k)+a/(k^alpha);
            else
                diff=0;
            end
            if diff>=0
                k=min(N,floor((a/prob(k))^(1/alpha)));
            end
            if diff<0 || k==N
                break
            end
        end
        if diff<0 
            if j<N
                j=k;
            else
                foundi=1;
            end
        else
            counter=counter+1;
            foundj=1;
            as=a;
            alphas=alpha;
            if  a/(key^alpha)-p > D || counter==10
                foundi=1;
            else
                D= a/(key^alpha)-p;
            end
            if j<N
                i=j;
            else
                foundi=1;
            end
        end
   end
end
end
