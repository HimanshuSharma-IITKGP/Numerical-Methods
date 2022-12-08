% solution of the ODE -k*(d^2T/dx^2) + -k*(d^2T/dy^2)  = q0  for 0<x<1 and 0<y<1 
% with boundary conditions T(x, 1) = T(1, y) = 0 and (dT/dn) on x=0 side and on y=0 side = 0

clear;
clc;
close all ;

syms x y ;

N = 1;
q0 = 1;
k = 1;



phi_0 = 0 ; % phi_0(0) = 0 and phi_0(1) = 0 ;

phi_array = createPhi(N, x, y);
f = createF(N, phi_array, phi_0, x, y, k, q0) ;
B = createB(N, phi_array, x, y, k);
c = pinv(B)*f ;

U = createSolution(phi_0, phi_array, N, c) ;
Usol = matlabFunction(U) ;
figure;
fcontour(U, linewidth=2) ;
% xlim([0, 1]) ;
% ylim([-0.01, 0.01]);

function [U] = createSolution(phi_0, phi_array, N, c)

    U = phi_0 ;
    for i=1:N
        U = U + c(i)*phi_array(i) ;
    end

    disp(U);
    U = simplify(U) ;
    disp(U) ;
end

function [phi_array] = createPhi(N, x, y)
    phi_array(1) = (1-x*x)*(1-y*y) ;

    for i = 2:N
        phi_array(i) = ((x*y)^(i-1))*phi_array(1) ;
    end
end



function [f] = createF(N, phi_array, phi_0, x, y, k, q0)
    f = zeros(N, 1) ;
    for i=1:N
        disp(linearFunctional(phi_array, phi_0, i, x, y, k, q0));
        f(i) = linearFunctional(phi_array, phi_0, i, x, y, k, q0) ;
    end
end



function [B] = createB(N, phi_array, x, y, k)
    B = zeros(N, N);

    for i=1:N
        for j=1:N
            disp(bilinearFunctional(phi_array, i, j, x, y, k)) ;
            B(i, j) = bilinearFunctional(phi_array, i, j, x, y, k);
        end
    end
end


function [f_i] = linearFunctional(phi_array, phi_0, i, x, y, k, q0)
    phi_i = phi_array(i) ;

    a = matlabFunction(q0*phi_i) ;
    b = matlabFunction(k*( diff(phi_i, x)*diff(phi_0, x) + diff(phi_i, y)*diff(phi_0, y)));
    b = 0; % hard coding the value of b because integral2 gives an error for b = 0
%     disp(a) ;
%     disp(b) ;
    f_i = integral2(a, 0, 1, 0, 1) - b  ;
end


function [B_ij] = bilinearFunctional(phi_array, i, j, x, y, k)
    phi_i = phi_array(i) ;
    phi_j = phi_array(j) ;

    c = matlabFunction(k*( diff(phi_i, x, 1)*diff(phi_j, x, 1) + diff(phi_i, y, 1)*diff(phi_j, y, 1)));
%     disp(c) ;
    B_ij = integral2(c, 0, 1, 0, 1) ;
end