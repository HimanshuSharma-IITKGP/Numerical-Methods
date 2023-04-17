
%  this algorithm computes the determinant by performing row operations on the matrix and  %
%  to convert it to a diagonal matrix   %
%  and then multiplying the diagonal terms %
function [det] = determinant(C)
    
    n = size(C, 1);
    o = (1:n)';
    
    det = 1;

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

   
    for i = 1:n
        pivot_row_ind = o(i) ;

        for j = 1:n
            if(i~=j)
                em =  C(o(j), pivot_row_ind)/C(pivot_row_ind, pivot_row_ind) ;
                C(o(j), :) = C(o(j), :) - em*C(pivot_row_ind, :);
            end
        end
        det = det*C(i, i);
    end

    disp(C)
end


