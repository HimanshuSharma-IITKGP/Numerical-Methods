% solution of the ODE -d^2u/dx^2 = -f for 0<x<L
% with boundary conditions du/dx (at x=0 and x=L) = 0
% using Trigonometric functions for approximation

clear;
clc;

syms x ;

L = 1;
N = 2;
f0 = 1;
fx = f0*cos(pi*x/L) ;




phi_0 = 0 ; % phi_0(0) = 0 and phi_0(1) = 0 ;
phi_array = [cos(pi*x/L), cos(2*pi*x/L)] ;
f = createF(N, phi_array, phi_0, x, L, fx) ;
B = createB(N, phi_array, x, L);
c = pinv(B)*f ;

U = createSolution(phi_0, phi_array, N, c) ;
Usol = matlabFunction(U) ;
figure;
fplot(x, U, linewidth=2) ;
% xlim([0, 1]) ;
% ylim([-0.05, 0.05]) ;

function [U] = createSolution(phi_0, phi_array, N, c)

    U = phi_0 ;
    for i=1:N
        U = U + c(i)*phi_array(i) ;
    end

    disp(U);
    U = simplify(U) ;
    disp(U) ;
end

function [phi_array] = createPhi(N, x, L)
    phi_array(1) = (x)*(x-L) ;

    for i = 2:N
        phi_array(i) = (x^(i-1))*phi_array(1) ;
    end
end

function [f] = createF(N, phi_array, phi_0, x, L, fx)
    f = zeros(N, 1) ;
    for i=1:N
        f(i) = linearFunctional(phi_array, phi_0, i, x, L, fx) ;
    end
end

function [B] = createB(N, phi_array, x, L)
    B = zeros(N, N);

    for i=1:N
        for j=1:N
            B(i, j) = bilinearFunctional(phi_array, i, j, x, L);
        end
    end
end
function [f_i] = linearFunctional(phi_array, phi_0, i, x, L, fx)
    phi_i = phi_array(i) ;

    f_i = -int(phi_i*fx, 0, L) - int(phi_i*diff(phi_0, x, 2), 0, L) ;
end


function [B_ij] = bilinearFunctional(phi_array, i, j, x, L)
    phi_i = phi_array(i) ;
    phi_j = phi_array(j) ;

    B_ij = int(phi_i*diff(phi_j, x, 2) , 0, L) ;
end