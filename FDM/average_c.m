omega = 1;
nu = 0.33;
a = 0;
b = 1;
E = 1;

C_av = zeros(size(C_transient)) ;

for i=1:size(C_transient, 1)
    for j=1:size(C_transient, 2)
        C_av(i, j) = integral(C_transient, i, j, h, a, n);
    end
end


U_transient = zeros(size(C_transient));
for i =1:size(U_transient, 1)
    for j=1:size(U_transient, 2)
        k = omega*(1+nu)/(3*(1-nu));
        U_transient(i, j) = -2*C_av(i, n)*omega/(9*(1-nu)) + k*(a+(j-1)*h)*C_av(i, j)/3 ;
    end
end

figure 
plot(r_arr, U_transient(total_iterations, :), linewidth=2);


function [I] = integral(C_transient, i ,j, h, a, n)
    I = 0;
    for iter = 1:j
        if(iter == 1 || iter == n)
            I = I + C_transient(i, iter)*(a + (iter-1)*h)^2;
        else
            I = I + 2*C_transient(i, iter)*(a + (iter-1)*h)^2;
        end
    end

    I = I*h/2;
    I = 3*I/(C_transient(i, j)^3) ;
end
