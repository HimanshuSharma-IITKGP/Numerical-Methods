% solution of d^2u/dx^2 = k (constant) using FDM
% with boundary conditions u(x1) = u1 and u(x2) = u2

%%
clear;
clc;

x1 = 0;
x2 = 1;
u1 = 0;
u2 = 1;
k = 5;
n = 30;

x_domain = x2 - x1;
num_points = n;
h = x_domain/(num_points-1) ;
x = x1:h:x2; 


u(1) = u1;
u(num_points) = u2;

u_new(1) = u1;
u_new(num_points) = u2;

error_now = 1 ;
error_req = 1e-7 ;

iterations = 0 ;

%%
while error_now > error_req
    for i=2:num_points-1
        u_new(i) = 0.5*(u(i-1) + u(i+1)) - 0.5*k*h*h  ;
    end

    error_now = 0;
    for i=2:num_points-1
        error_now = error_now + abs(u(i) - u_new(i)) ;
    end
    iterations = iterations + 1;
    u = u_new ;
end

%%
plot(x, u, linewidth=2 );
xlim([0, 1]);
ylim([-0.5, 1]) ;
disp(u) ;
