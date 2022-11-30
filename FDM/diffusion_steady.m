clear ;
close all;
clc;

%% initialization

x_dom_size = 1;
y_dom_size = 1;

x_num_points = 5;
y_num_points = 5;

deltaX = x_dom_size/(x_num_points - 1) ;
deltaY = y_dom_size/(y_num_points - 1) ;



T = zeros(x_num_points, y_num_points) ;
T(1, :) = 1;

T_new = zeros(x_num_points, y_num_points) ;
T_new(1, :) = 1;


error_now = 1;
error_req = 10e-4;

%% Iterations
iterations = 1;

while error_now > error_req
    for i=2:(x_num_points - 1)
        for j=2:(y_num_points - 1)
            a = 0.25*(T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1));
            T_new(i, j) = a ;
        end
    end


    iterations = iterations + 1;
    error_now = 0;
    for i=2:(x_num_points - 1)
        for j=2:(y_num_points - 1)
            error_now = error_now + abs(T(i, j)-T_new(i, j));
        end
    end
    T = T_new ;
end

%% Plotting
disp(T_new);
x_dom = 0:deltaX:x_dom_size ;
y_dom = 0:deltaY:y_dom_size ;
[X,Y] = meshgrid(x_dom,y_dom);
contourf(X, Y, T_new) ;
colorbar