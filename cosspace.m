function c = cosspace(a, b, n)
% a - first point
% b - end point
% n - number of points
c = ((a + b)/2 + (b - a)/2 )* cos(linspace(-pi,0,n));
