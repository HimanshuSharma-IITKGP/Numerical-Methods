clear ;
close all;
clc;

%% initialization

dom_size = 1;
n_points = 6;
h = dom_size/(n_points - 1) ;
 dt = 0.01 ;
 alpha = dt/(h*h) ;


T = zeros(n_points, n_points) ;
T(1, :) = 1;

T_new = zeros(n_points, n_points) ;
T_new(1, :) = 1;

T_transient = 0;

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
        iterations %#ok<*NOPTS> 
        error_mag
    end
    error_track(iterations) = error_mag;
    T = T_new ;
    T_transient(1:n_points, 1:n_points, iterations) = T_new ;
end

%% Plotting

x = 0:h:dom_size;
y = x;
[X, Y] = meshgrid(x, y) ;

for iter=1:iterations
    contourf(X, Y, T_transient(:, :, iter)) ;
    colorbar;
    disp(iter);
    pause(0.2) ;
end
