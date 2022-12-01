% solution of the ODE d((1+x)du/dx)dx = 0 for 0<x<1
% with boundary conditions u(0) = 0 and u(1) = 1


clear;
clc;


x1 = 0;
x2 = 1;
u1 = 0;
u2 = 1;
N = 2;

syms x ;

phi_0 =  ((x2-x1)*(x-x1)/(u2-u1)) + u1; % phi_0(x1) = u1 and phi_0(x2) = u2;
phi_array = createPhi(N, x, x1, x2);
f = createF(N, phi_array, phi_0, x, x1, x2) ;
B = createB(N, phi_array, x, x1, x2);
c = pinv(B)*f ;

U = createSolution(phi_0, phi_array, N, c) ;
Usol = matlabFunction(U) ;
figure;
fplot(x, U, linewidth=2) ;
xlim([0, 1]) ;
ylim([0, 1]) ;

function [U] = createSolution(phi_0, phi_array, N, c)

    U = phi_0 ;
    for i=1:N
        U = U + c(i)*phi_array(i) ;
    end

    disp(U);
    U = simplify(U) ;
    disp(U) ;
end

function [phi_array] = createPhi(N, x, x1, x2)
    phi_array(1) = (x-x2)*(x-x1) ;

    for i = 2:N
        phi_array(i) = ((x-x1)^(i-1))*phi_array(1) ;
    end
end

function [f] = createF(N, phi_array, phi_0, x, x1, x2 )
    f = zeros(N, 1) ;
    for i=1:N
        f(i) = linearFunctional(phi_array, phi_0, i, x, x1, x2 ) ;
    end
end

function [B] = createB(N, phi_array, x, x1, x2)
    B = zeros(N, N);

    for i=1:N
        for j=1:N
            B(i, j) = bilinearFunctional(phi_array, i, j, x, x1, x2);
        end
    end
end
function [f_i] = linearFunctional(phi_array, phi_0, i, x, x1, x2)
    phi_i = phi_array(i) ;
    didx = diff(phi_i, x) ;
    d0dx = diff(phi_0, x) ;

    f_i = -int(didx*d0dx*(1+x), x1, x2) ;
end


function [B_ij] = bilinearFunctional(phi_array, i, j, x, x1, x2)
    phi_i = phi_array(i) ;
    phi_j = phi_array(j) ;

    didx = diff(phi_i, x) ;
    djdx = diff(phi_j, x) ;

    B_ij = int(didx*djdx*(1+x), x1, x2) ;
end