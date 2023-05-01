% a_i*x(i) = b_i*x(i+1) + c_i*x(i-1) + d_i


% a, b, c, d are column vectors
function [x] = TDMA(a, b, c, d)

x = zeros(size(a)) ;
n = size(a, 1);


P = zeros(size(a));
Q = zeros(size(a));

% forward substitution
P(1) = b(1)/a(1);
Q(1) = d(1)/a(1);

for i = 2:n
    m = (a(i)-c(i)*P(i-1));
    P(i) = b(i)/m ;
    Q(i) = (d(i)+c(i)*Q(i-1))/m ;
end


% backward substitution
x(n) = Q(n);

for i = (n-1):-1:1
    x(i) = P(i)*x(i+1) + Q(i) ;
end

end