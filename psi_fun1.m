function psi1 = psi_fun1(z) %Eq. 3.18
%from 6.3.18 (p. 80 - Handbook of Mathematical Functions) -> Polygamma
%reduces to gamma (6.4.1 p. 81) followed by asymptotic expansion.
sum=0;
for j=1:100
    sum_old=sum;
    
    if j==1
            sum=sum + bernoulli(2*j)/(2*j*z^(2*j));
    else
            sum=sum + bernoulli(2*j)/(2*j*z^(2*j));
    end
    %Convergence criteria - change if needed
    if abs((sum-sum_old)/sum_old)<0.00000000000001 
       break 
    end 
end

if j>=100
    'Convergence is not reached'
end 

psi1=log(z)-(1/(2*z))-sum; %eq. 3.18 (log in matlab = ln)