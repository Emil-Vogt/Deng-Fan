function alpha = alpha_fun(n,k,b,d) %Eq. 3.4
alpha=vpa(((k*b*(b+2))/(n+d)-n-d)/2);