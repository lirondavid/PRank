   function [upper,pareto_time] = pareto_estimation(probMatrix,index,d,n,a,alpha)

    stime=tic;
        
    %Computing prob for pareto
    p=1;
    for i=1:d
        p=p*probMatrix(i,index(i));
    end 
  
  
    %Computing shorter bounds   
    n2=zeros(1,d);    
    for i=1:d        
        probMul=1;
        k=n(i);
        for j=1:d
            if j~=i
                probMul=probMul*probMatrix(j,1);
            else
                probMul=probMul*probMatrix(i,k); 
            end
        end
        
        if probMul>=p
            n2(i)=n(i);
        else    
    
        first=1;
        last=n(i);        
        while(first+1<last)
            middle=floor((first+last)/2);
            probMul=1;
                for j=1:d
                    if j~=i
                        probMul=probMul*probMatrix(j,1);
                    else
                        probMul=probMul*probMatrix(i,middle); 
                    end     
                end
             if probMul>=p
                 first=middle;   
             end
             if probMul<p
                 last=middle;
             end    
        end
        n2(i)=first;
        end 
    end   
    
   
    B=0;
    j=1;
    while(j<=d)
        B=B+log2(a(j))-alpha(j)*log2(n2(j));
        j=j+1;
    end
    B=B-log2(p);
    B=2^B;
    
    if B>1
        upper=prod(n2)-prod(n2-index);
    else
        upper=upper_rank(alpha,a,d,p,n2);
    end
    pareto_time=toc(stime);
  
 
end
    
    