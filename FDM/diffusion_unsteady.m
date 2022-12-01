% solution of the 2-D unsteady diffusion Equation
% Boundary conditions : T = 1 at the upper boundary and T = 0 at other 3
% boundaries
% Initial  condition  : T = 0 initially everywhere
clear ;
close all;
clc;

%% initialization

dom_size = 1;
n_points = 51;
h = dom_size/(n_points - 1) ;
dt = 0.0001 ;
alpha = dt/(h*h) ;


T = zeros(n_points, n_points) ; 
T(1, :) = 1; 

T_new = zeros(n_points, n_points) ;
T_new(1, :) = 1;


error_mag = 1;
error_req = 10e-6;
iterations = 0;
error_track = 0;


%% Iterations


while error_mag > error_req
    for i=2:(n_points - 1)
        for j=2:(n_points - 1)
            a = T(i, j) + alpha*(T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1) - 4*T(i, j));
            T_new(i, j) = a ;
        end
    end


    iterations = iterations + 1;

    error_mag = 0;

    for i=2:(n_points - 1)
        for j=2:(n_points - 1)
            error_mag = error_mag + abs(T(i, j)-T_new(i, j));
        end
    end

    if(rem(iterations, 1000) == 0)
        iterations 
        error_mag
    end
    error_track(iterations) = error_mag;
    T = T_new ;
end

%% Plotting

x_dom = 0:h:dom_size ;
y_dom = x_dom ;

[X,Y] = meshgrid(x_dom,y_dom);
contourf(X, Y, T_new) ;
colorbar