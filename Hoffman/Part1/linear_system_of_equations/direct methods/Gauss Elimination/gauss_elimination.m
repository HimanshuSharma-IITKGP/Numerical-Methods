


function [C, o] = gauss_elimination(A, b)
    C = [A, b];
    n = size(A, 1);
    o = (1:n)';
    
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
        end
    end
end

