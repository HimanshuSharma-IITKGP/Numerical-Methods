clear;
clc;


hold all;

% legend;
for i = 0:10
    disp(expand(leg(i)));
end
grid on

function [y] = leg(n)
    N = floor(n/2) ;


    syms x;
    syms y;
    y = 0;
    for m=0:N
        y = y + ((-1)^m*factorial(2*(n-m))*x^(n-2*m))/((2^n)*factorial(m)*factorial(n-m)*factorial(n-2*m));
    end
  
    fplot(y, [-1, 1], linewidth=2);
end