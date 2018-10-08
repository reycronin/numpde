function [t,L,rhs] = laplacian_cheb_f0(n,a,b)

x = cosspace(-1,1,n);

% transformations to [0,1] grid t
t = (x+1)./2;
t1 = 2; % inverse first der of t for first derivative equal to dx/dt

T0 = chebeval(x, length(x) ,0);
T1 = chebeval(x, length(x) ,1).*t1;
T2 = chebeval(x, length(x) ,2).*(t1^2);

L = T2 + a.*T1 + b.*T0;

% first row is x = 0
% for u(0) = 1
% for u'(0) = 0

L(end,:)=T1(1,:);
L(1,:) = T0(1,:);



rhs = rhsfunc(t,2) + a.*rhsfunc(t,1) + b.*rhsfunc(t,0);
rhs = rhs*0;
% rhs(ix)= 0;
rhs(1) = 1;
T = T0;
invT = inv(T);

L=L*invT;

    
