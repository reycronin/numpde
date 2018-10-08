close all 
clear

xx = linspace(0,1,100);
a=2;
b=2;
n=50;

%% solving the ODE

[t, L, rhs] = laplacian_cheb(n,a,b);
uu = inv(L)*rhs'; 


figure(1)
plot(xx, rhsfunc(xx,0),'k',t,uu ,'r.');
title('ODE Solution')
legend('exact solution','numerically solved solution', 'Location', 'SouthEast')
xlabel('t')
ylabel('y(t)')

rhs_t=rhsfunc(t,0);

error = norm(uu-rhs_t',inf);

%% convergence study

num = 10; %number of trials
n = 10:num:100;
errors = zeros(1,length(n));
count=0;

for nn = n
count=count+1;
[x, L, rhs] = laplacian_cheb(nn,a,b);
uu = L\rhs'; % u double prime
errors(1,count) = norm(uu-rhsfunc(x,0)',inf);
end

figure(2)

semilogy(n,errors,'r.',n,1./n.^2,n,1./n.^3,n,1./n.^4)
xlabel('number of points')
ylabel('error')
title('convergence study')
legend('errors', '1/n^2','1/n^3', '1/n^4')

% this plot shows there is spectral convergence, the errros converge faster
% than any of the polynomial functions


%% f(t) = 0

n=50;
count = 0;

uu=zeros(n,3);
for a=[1,20,400]
    count = count+1;
b=100;
[t,L,rhs] = laplacian_cheb_f0(n,a,b);
uu(:,count) = inv(L)*rhs'; 
end

figure(3)
plot(t,uu(:,1), t,uu(:,2), t,uu(:,3))
title('b=100')
legend('a=1','a=20', 'a=400' )
xlabel('t')





