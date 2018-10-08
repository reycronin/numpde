function T = vander_chebyshev(x,n)
% This function creates a matrix with the columns as the chebyshev
% polynomials for each point of the input vector

% x - input points (rows)
% n - number of chebyshev polynomials (columns)

lx = length(x);

if ~exist('n','var')
    % if no input n then make V a square matrix by using length of x
    n = lx;
end

T = ones(lx,n);
if n > 1
    T(:,2) = x;
end

x=x'; % make x a column vector
for i = 3:n
    T(:,i) = 2.*x.*T(:,i-1) - T(:,i-2);
end


