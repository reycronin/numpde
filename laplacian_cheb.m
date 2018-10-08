function [t,L,rhs] = laplacian_cheb(n,a,b, bound)


% n - number of points to evaluate
% bounds_matrix each row is a boundary
% [left/right , Dirichlet/Neumann]
% [1,1] -> u(0)'' = u'(0) [left, Neumann]
% [n,0] -> u(n)'' = u(n) [right, Dirichlet]
% laplacian_cheb(50, [1,0;50,0]) is the example from FDHighOrder.ipynb

x = cosspace(-1,1,n);

% transformations to [0,1] grid t
t = (x+1)./2;
t1 = 2; % inverse first der of t for first derivative equal to dx/dt


T0 = chebeval(x, length(x) ,0);
T1 = chebeval(x, length(x) ,1).*t1;
T2 = chebeval(x, length(x) ,2).*(t1^2);

% first row is x=0
% for u(0) = 1
% for u'(0) = 0

L = T2 + a.*T1 + b.*T0;
%L = T2 + a.*T1 + b.*T0;

L(end,:)=T1(1,:);
L(1,:) = T0(1,:);

rhs = rhsfunc(t,2) + a.*rhsfunc(t,1) + b.*rhsfunc(t,0);
rhs(end)= rhsfunc(t(1),1);
rhs(1) = rhsfunc(t(1),0);

T = T0;
invT = inv(T);

if exist('bound','var')
    [row,~]=size(bound);
    for i = 1:row
        T = chebeval(x,length(x),bound(i,2));
        L(bound(i,1),:) = T(bound(i,1),:);
        rhs(bound(i,1)) = rhsfunc(x(bound(i,1)),bound(i,2));
    end
end


L=L*invT;

    
