function C0 = C0_fun(u,v,d,m,n,a,k,b) %Summand of eq. 3.11 for lambda = 0
%The numerator and denominator are first computed without the betafunction
Numerator1=vpa(Normalization_fun(v,k,b,d,a)*Normalization_fun(u,k,b,d,a)*(-1)^(m+n)*pochhammer(-v,m)*pochhammer(-u,n)*pochhammer(1-v-2*d,m)*pochhammer(1-u-2*d,n));
Denominator1=vpa(a*pochhammer((1+2*alpha_fun(v,k,b,d)),m)*pochhammer((1+2*alpha_fun(u,k,b,d)),n)*gamma(n+1)*gamma(m+1));
%The beta function is then included for the final calculation
C0=vpa(Numerator1/Denominator1*beta(vpa(2*d+u+v-n-m+1),vpa(alpha_fun(u,k,b,d)+alpha_fun(v,k,b,d)+m+n)));