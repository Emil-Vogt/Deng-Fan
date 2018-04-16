%Remember to have all Deng-Fan scripts in the same folder.
%This script is used to calculated oscillator strengths for the Deng-Fan
%potential according to the article "Accuracy of XH-stretching Intensities 
%with the Deng-Fan Potential" - by Emil Vogt, Daniel Sage and 
%Henrik G. Kjaergaard. All references are to equations in the original
%manuscript. To use this script dipole moment functions need to
%be expressed as sixth-order polynomials in the displacement coordinate
%and parameters of the Deng-Fan potential should be specified.
%These can be obtained by fitting eq. 3.7 to transition frequencies
%calculate by other methods. This script is written by Emil Vogt and
%Aleksandrs Smilgins.

% Load constants
format compact
%Use 'format long g' if needed.
digits(64)
%Physical constants
N_A = 6.022140857*10^23;
hhat=1.054571726*10^(-34);
h=6.62607004*10^(-34);
c=vpa(physconst('LightSpeed')*10^2); %Conversion from m/s to cm/s

%% All parameters in this section need to be put in manually. - START
%Parameters - Should be obtained by fitting eq. 3.7 to the available
%transition frequencies (calculated from eq. 2.5 + 2.4 in the article)
a=1.585082*10^10;
re=0.9701*10^(-10);
D=vpa(48625.05407*h*c);

%Specify masses in atomic units
%Hydrogen
mass_1 = 1.00782503223;
%Other atom Cl(34.96885268),F(18.99840322),O(15.99491461956)
mass_2 = 15.99491461956;

%Dipolemoment coefficients (units of D/Å, D/Å^2,...,D/Å^6)
x=[0.0; -0.00001 ; 0.00005; -0.00006; -0.00025; 0.00028];
y=[1.03327; -0.63141; -0.75318; -0.24755; 0.20279; 0.12008];
z=[0.40192; 0.27323; -0.42830; -0.16269; -0.16443; 0.16097];
%Do not remove a component but set zero if not used!
%The section for defining parameters ends here.

%% Calculate new parameters based on the input parameters
%Conversion from atomic units to g/mol
ma_1=mass_1*10.^(-3.0)/N_A;
ma_2=mass_2*10.^(-3.0)/N_A;
%Reduced mass
my=ma_1*ma_2/(ma_1+ma_2);

%Intermediate constants
b=vpa(exp(a*re)-1);
k=vpa(2*my*D/(a^2*hhat^2));
d=vpa(0.5*(1+(1+4*k*b^2)^(0.5)));


%% Transition frequencies - verify your transition frequencies from fit!
'Transition Frequencies'
for n=1:5
    (D*(n+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(n+d)))^2-1))*n
end

%% matrix elements
%This part of the script deals with the matrix elements of integral powers
%of the coordinate (eq. 3.17 + 3.21-3.25 for u=0 and v=1...5)
u=0;

res=zeros(5,6);

for v=1:5;
    v;
    for n=0:u;
        n;
        for m=0:v;
            res(v,1)=res(v,1)-C0_fun(u,v,d,m,n,a,k,b)*Ci_fun(1,u,v,d,m,n,a,k,b); %eq. 3.17
            res(v,2)=res(v,2)+C0_fun(u,v,d,m,n,a,k,b)*(Ci_fun(1,u,v,d,m,n,a,k,b)^2+Ci_fun(2,u,v,d,m,n,a,k,b)); %eq. 3.21
            res(v,3)=res(v,3)-C0_fun(u,v,d,m,n,a,k,b)*(Ci_fun(1,u,v,d,m,n,a,k,b)^3+3*Ci_fun(1,u,v,d,m,n,a,k,b)*Ci_fun(2,u,v,d,m,n,a,k,b)+Ci_fun(3,u,v,d,m,n,a,k,b)); %eq. 3.22
            
            res(v,4)=res(v,4)+C0_fun(u,v,d,m,n,a,k,b)*(Ci_fun(1,u,v,d,m,n,a,k,b)^4+6*Ci_fun(1,u,v,d,m,n,a,k,b)^2*Ci_fun(2,u,v,d,m,n,a,k,b)+4*Ci_fun(1,u,v,d,m,n,a,k,b)*Ci_fun(3,u,v,d,m,n,a,k,b)+...
               3*Ci_fun(2,u,v,d,m,n,a,k,b)^2+Ci_fun(4,u,v,d,m,n,a,k,b)); %eq. 3.23
            
            res(v,5)=res(v,5)-C0_fun(u,v,d,m,n,a,k,b)*(Ci_fun(1,u,v,d,m,n,a,k,b)^5+10*Ci_fun(1,u,v,d,m,n,a,k,b)^3*Ci_fun(2,u,v,d,m,n,a,k,b)+...
                10*Ci_fun(1,u,v,d,m,n,a,k,b)^2*Ci_fun(3,u,v,d,m,n,a,k,b)+15*Ci_fun(1,u,v,d,m,n,a,k,b)*Ci_fun(2,u,v,d,m,n,a,k,b)^2+...
                10*Ci_fun(2,u,v,d,m,n,a,k,b)*Ci_fun(3,u,v,d,m,n,a,k,b)+5*Ci_fun(1,u,v,d,m,n,a,k,b)*Ci_fun(4,u,v,d,m,n,a,k,b)+Ci_fun(5,u,v,d,m,n,a,k,b)); %eq. 3.24
            
            res(v,6)=res(v,6)+C0_fun(u,v,d,m,n,a,k,b)*(Ci_fun(1,u,v,d,m,n,a,k,b)^6+15*Ci_fun(1,u,v,d,m,n,a,k,b)^4*Ci_fun(2,u,v,d,m,n,a,k,b)+20*Ci_fun(1,u,v,d,m,n,a,k,b)^3*Ci_fun(3,u,v,d,m,n,a,k,b)+...
                45*Ci_fun(1,u,v,d,m,n,a,k,b)^2*Ci_fun(2,u,v,d,m,n,a,k,b)^2+15*Ci_fun(1,u,v,d,m,n,a,k,b)^2*Ci_fun(4,u,v,d,m,n,a,k,b)+...
                60*Ci_fun(1,u,v,d,m,n,a,k,b)*Ci_fun(2,u,v,d,m,n,a,k,b)*Ci_fun(3,u,v,d,m,n,a,k,b)+6*Ci_fun(1,u,v,d,m,n,a,k,b)*Ci_fun(5,u,v,d,m,n,a,k,b)+...
                15*Ci_fun(2,u,v,d,m,n,a,k,b)*Ci_fun(4,u,v,d,m,n,a,k,b)+10*Ci_fun(3,u,v,d,m,n,a,k,b)^2+Ci_fun(6,u,v,d,m,n,a,k,b)); %eq. 3.25
        end
    end
end

'The result is:'
res

%% Last part of the binomial section
%In matrix res we have v as label for different rows, and r^1, r^2, r^3 as
%columns. A new matrix is made that has v labels as rows and j as columns. 
%Nevertheless, we leave out j=0 and v=0.

'START'
res2=nan(5,6);
for v1=1:5
    r_vector=res(v1,:); %Here we have r1, r2, ..., r6 for a v1 value
    
    for j=1:6 %columns
        sum=0;
        sum=r_vector(j); %first term
        
        if j>=2 %binomial
            for q=1:j-1 %The binomial coefficient
                sum=sum+nchoosek(j,q)*r_vector(j-q)*re^(q)*(-1)^q;
            end
        end
        
        %Last term is removed due to orthogonality of the Deng-Fan
        %wavefunctions (sum=sum+re^j*(-1)^j)
        
        res2(v1,j)=sum*(10^10)^j;
    end
end

'Displacement matrix elements'
res2

%% Last part (calculate oscillator strengths)
res3=nan(5,6);

%Remember the transition frequencies (eq. 3.7)
for rk=1:5
    %y-component
    res3(rk,1)=4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(y(1)*res2(rk,1))^2;
    res3(rk,2)=4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(y(1)*res2(rk,1)+y(2)*res2(rk,2))^2;
    res3(rk,3)=4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(y(1)*res2(rk,1)+y(2)*res2(rk,2)+y(3)*res2(rk,3))^2;
    res3(rk,4)=4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(y(1)*res2(rk,1)+y(2)*res2(rk,2)+y(3)*res2(rk,3)+y(4)*res2(rk,4))^2;
    res3(rk,5)=4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(y(1)*res2(rk,1)+y(2)*res2(rk,2)+y(3)*res2(rk,3)+y(4)*res2(rk,4)+y(5)*res2(rk,5))^2;
    res3(rk,6)=4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(y(1)*res2(rk,1)+y(2)*res2(rk,2)+y(3)*res2(rk,3)+y(4)*res2(rk,4)+y(5)*res2(rk,5)+y(6)*res2(rk,6))^2;

    %x-component
    res3(rk,1)=res3(rk,1)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(x(1)*res2(rk,1))^2;
    res3(rk,2)=res3(rk,2)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(x(1)*res2(rk,1)+x(2)*res2(rk,2))^2;
    res3(rk,3)=res3(rk,3)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(x(1)*res2(rk,1)+x(2)*res2(rk,2)+x(3)*res2(rk,3))^2;
    res3(rk,4)=res3(rk,4)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(x(1)*res2(rk,1)+x(2)*res2(rk,2)+x(3)*res2(rk,3)+x(4)*res2(rk,4))^2;
    res3(rk,5)=res3(rk,5)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(x(1)*res2(rk,1)+x(2)*res2(rk,2)+x(3)*res2(rk,3)+x(4)*res2(rk,4)+x(5)*res2(rk,5))^2;
    res3(rk,6)=res3(rk,6)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(x(1)*res2(rk,1)+x(2)*res2(rk,2)+x(3)*res2(rk,3)+x(4)*res2(rk,4)+x(5)*res2(rk,5)+x(6)*res2(rk,6))^2;

    %z-component
    res3(rk,1)=res3(rk,1)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(z(1)*res2(rk,1))^2;
    res3(rk,2)=res3(rk,2)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(z(1)*res2(rk,1)+z(2)*res2(rk,2))^2;
    res3(rk,3)=res3(rk,3)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(z(1)*res2(rk,1)+z(2)*res2(rk,2)+z(3)*res2(rk,3))^2;
    res3(rk,4)=res3(rk,4)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(z(1)*res2(rk,1)+z(2)*res2(rk,2)+z(3)*res2(rk,3)+z(4)*res2(rk,4))^2;
    res3(rk,5)=res3(rk,5)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(z(1)*res2(rk,1)+z(2)*res2(rk,2)+z(3)*res2(rk,3)+z(4)*res2(rk,4)+z(5)*res2(rk,5))^2;
    res3(rk,6)=res3(rk,6)+4.702*10^(-7)*(D*(rk+2*d)/(4*k*h*c)*(((b*(b+2)*k)/(d*(rk+d)))^2-1))*rk*(z(1)*res2(rk,1)+z(2)*res2(rk,2)+z(3)*res2(rk,3)+z(4)*res2(rk,4)+z(5)*res2(rk,5)+z(6)*res2(rk,6))^2;
end

'row = Oscillator strength for the transition 0->1,2...5, column = truncation of dipole moment'
res3

