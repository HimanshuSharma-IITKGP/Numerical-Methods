% solution of the equation (dC/dt) - D*(d^2C/dx^2) = 0, t>0 and 0<=x<=1 ;
% with the boundary conditions dT/dx = 0 at x = 0 and dT/dx = k at x = 1;


function [C_transient] = diffusion(D, a, b, n, L_c, D_c, k, dt, tFinal)
    h = (b-a)/(n-1) ;
    x_arr = a:h:b ;

    T_c = (L_c^2)/D_c ; 
    
    
    total_iterations = tFinal/dt ;
    disp(total_iterations);
    
    alpha = D*dt/(h*h) ;
    disp(alpha);
    
    C = zeros(n, 1) ;
    
    C_new = zeros(n, 1) ;

    C_transient = 0;
    
    for iter = 1:total_iterations
    
 
        for i = 1:n
            if i == 1
                C_new(i) = C(i) + 2*alpha*(C(i+1) - C(i)) ;
            elseif i == n
                C_new(i) = C(i) + 2*alpha*(h*k - C(i) + C(i-1)) ;
            else
                C_new(i) = C(i)  + alpha*(C(i+1) - 2*C(i) + C(i-1)) ;
            end
        end
    
    
        C = C_new ;
        C_transient(iter, 1:n) = C ;
    
%         plot(x_arr, C, linewidth=2) ;
%         ylim([0, 2]) ;
%         pause(0.1) ;



    end

    figure ;
    plot(x_arr, C, linewidth=2) ;
%     title(['D = ', num2str(D), ' At final time']);
    title(['k = ', num2str(k), ' At final time']);
    ylim([0, 2]) ;
end


