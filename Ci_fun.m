function Ci = Ci_fun(i,u,v,d,m,n,a,k,b) %eq. 3.16 and 3.19
if i==1 %If i=1 => use eq. 3.16
    Ci=vpa((psi_fun1(alpha_fun(u,k,b,d)+alpha_fun(v,k,b,d)+m+n)-psi_fun1(alpha_fun(u,k,b,d)+alpha_fun(v,k,b,d)+u+v+2*d+1))/(a));
else %If i>1 (C0 is computed with C0_fun) => use eq. 3.19 
    Ci=vpa((psi_fun(alpha_fun(u,k,b,d)+alpha_fun(v,k,b,d)+m+n,i-1)-psi_fun(alpha_fun(u,k,b,d)+alpha_fun(v,k,b,d)+u+v+2*d+1,i-1))/(a^(i)));
end