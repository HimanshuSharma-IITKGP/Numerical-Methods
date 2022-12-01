% solution of the ODE d^2(EI*(d^2w/dx^2))/dx^2  = fx for 0<x<L
% with boundary conditions W(0) = 0, w(L) = 0, EI*d^2w/dx^2 (at x=0) = 0
% and EI*d^2w/dx^2 (at x=L) = 0 using Trigonometric functions for
% approximation

clear;
clc;

syms x ;

L = 1;
N = 2;
f0 = 1;
fx = f0*dirac(x-0.5*L) ;
E = 1;
I = 1;



phi_0 = 0 ; % phi_0(0) = 0 and phi_0(1) = 0 ;
phi_array = [sin(pi*x/L), sin(2*pi*x/L)] ;
f = createF(N, phi_array, phi_0, x, L, E, I, fx) ;
B = createB(N, phi_array, x, L, E, I);
c = pinv(B)*f ;

U = createSolution(phi_0, phi_array, N, c) ;
Usol = matlabFunction(U) ;
figure;
fplot(x, U, linewidth=2) ;
xlim([0, 1]) ;
ylim([0, 0.025]);

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

function [f] = createF(N, phi_array, phi_0, x, L, E, I, fx)
    f = zeros(N, 1) ;
    for i=1:N
        f(i) = linearFunctional(phi_array, phi_0, i, x, L, E, I, fx) ;
    end
end

function [B] = createB(N, phi_array, x, L, E, I)
    B = zeros(N, N);

    for i=1:N
        for j=1:N
            B(i, j) = bilinearFunctional(phi_array, i, j, x, L, E, I);
        end
    end
end
function [f_i] = linearFunctional(phi_array, phi_0, i, x, L, E, I, fx)
    phi_i = phi_array(i) ;
    didx = diff(phi_i, x) ;
    d0dx = diff(phi_0, x) ;

    f_i = int(phi_i*fx, 0, L) - int(E*I*didx*d0dx, 0, L) ;
end


function [B_ij] = bilinearFunctional(phi_array, i, j, x, L, E, I)
    phi_i = phi_array(i) ;
    phi_j = phi_array(j) ;

    didx = diff(phi_i, x, 2) ;
    djdx = diff(phi_j, x, 2) ;

    B_ij = int(E*I*didx*djdx, 0, L) ;
end