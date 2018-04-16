function psi = psi_fun(z,i) %Eq. 3.20
%Asymptotic expansion of the polygamma function
sum=0;
for j=1:100
    sum_old=sum; %Calculate sum in eq. 3.20
    sum=sum + bernoulli(2*j)*gamma(2*j+i)/(gamma(2*j+1)*z^(2*j+i));
    %Convergence criteria - change if needed
    if abs((sum-sum_old)/sum_old)<0.00000000000001 
       break 
    end 
end

if j>=100
    'Convergence is not reached'
end 
%Compute final expression when the sum has converged
psi=(-1)^(i-1)*(gamma(i)/(z^i)+gamma(i+1)/(2*z^(i+1))+sum);