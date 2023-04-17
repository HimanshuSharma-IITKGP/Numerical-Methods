
% inverse of a matrix using row operations %
function [invA] = matrix_inverse(A)
    n = size(A, 1);

    I = eye(n);
    C = [A, I];
    
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

    disp(o)
    for i = 1:n
        pivot_row_ind = o(i) ;

        for j = 1:n
            if(i~=j)
                em =  C(o(j), pivot_row_ind)/C(pivot_row_ind, pivot_row_ind) ;
                C(o(j), :) = C(o(j), :) - em*C(pivot_row_ind, :);
            end
        end

        C(o(i), :) = C(o(i), :)/C(o(i), o(i));
    end

    invA = C(1:n, n+1:2*n);
   
end



