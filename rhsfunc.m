function rhs = rhsfunc(x , der)
% x - points to be evaluated
% der - derivative inputs are 0 for function, 1 for first derivative and 2
% for second derivative

k = 5;
x0 = 0.3; %centering

if der == 0
    rhs = tanh(k.*(x-x0));
end

if der == 1
    rhs = k.*cosh(k.*(x-x0)).^(-2);
end

if der == 2
    rhs = -2.*k^2*tanh(k.*(x-x0)).*cosh(k.*(x-x0)).^-2;
end

    
