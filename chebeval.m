function cheb = chebeval(z,n,der)

% z - point you want to evaluate
% n - number of chubby chef polynomials
% der - derivative input

if (der ~= 0 && der ~= 1 && der ~= 2)
    error('derivative input must be equal to 0 for interpolation, 1 for first derivative or 2 for second derivative')
end


    
if ~exist('n','var')
    % if no input n then make V a square matrix by using length of x
    n = length(z);
end

Tz = vander_chebyshev(z,n);
dTz = zeros(size(Tz)); % first column is der of 1 -> all zeros
% initialize:
dTz(:,2) = 1;
dTz(:,3) = 4*z;
ddTz = zeros(size(Tz)); % 1st and 2nd columns der of 0 and 1 = 0
% initialize:
ddTz(:,3) = 4;
for i = 4:n
    num = i-1; % num is the n in Tn
    %first der
    dTz(:,i) = num*(2*Tz(:,i-1) + dTz(:,i-2)/(num-2));
    %second der
    ddTz(:,i) = num*(2*dTz(:,i-1) + ddTz(:,i-2)/(num-2));
end

% isolate the output
if der == 0
    cheb = Tz;
end
if der == 1
    cheb = dTz;
end
if der == 2
    cheb = ddTz;
end



