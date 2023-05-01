% A = [80 -20 -20; -20 40 -20; -20 -20 130];
% [L, U] = Doolittle_LU(A);
% L
% U

function [L, U] = Doolittle_LU_Decomposition(C)
   
    n = size(C, 1);
    o = (1:n)';L8

    L = eye(n) ;
    
    for i = 1:n-1
        max = -1;
        s = -1;
        for g=i:n
            temp = abs(C(g, i))/getMax(abs(C(g, :)'));
            if (temp>max)
                max = temp;
                s = g;
            end
        end
    
        temp = o(i);
        o(i) = o(s);
        o(s) = temp;
    end

    for i = 1:n-1
        pivot_row_ind = o(i) ;

        for j = i+1:n
                em =  C(o(j), pivot_row_ind)/C(pivot_row_ind, pivot_row_ind) ;
                C(o(j), :) = C(o(j), :) - em*C(pivot_row_ind, :);

                L(o(j), o(i)) = em ;
        end
    end

    U = C(1:n, 1:n);
end

