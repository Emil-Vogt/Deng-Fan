function N = Normalization_fun(n,k,b,d,a) %Eq. 3.5
%The numerator and denominator are first computed without the betafunction
taeller=a*(alpha_fun(n,k,b,d)+n+d)*gamma(2*alpha_fun(n,k,b,d)+n+1);
naevner=(n+d)*gamma(n+1)*gamma(2*alpha_fun(n,k,b,d)+1);
%The beta function is then included for the final calculation
N=(taeller/(naevner*beta(2*alpha_fun(n,k,b,d),n+2*d)))^(0.5);