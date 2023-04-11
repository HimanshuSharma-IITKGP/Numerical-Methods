clear;
clc;

A = [3 2 105;2 -3 103; 1 1 3];
b = [104; 98; 3];

n = size(A, 1);
[C, o, X] = eliminate(A, b);
X

function [C, o, X] = eliminate(A, b)
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

    X = C(:,n+1);
end



function [m] = getMax(X)
    m = -1;
    for i = 1:size(X)
        if(m<X(i))
            m = X(i);
        end
    end
end