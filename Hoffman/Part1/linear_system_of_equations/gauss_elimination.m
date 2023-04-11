clear;
clc;

A = [3 2 105; 2 -3 103; 1 1 3];
b = [104; 98; 3];

[C, o] = eliminate(A, b);
X = backSubstitution(C, o)

function [C, o] = eliminate(A, b)
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

function [X] = backSubstitution(C, o)
    n = size(o, 1);
    X = zeros(size(o)) ;
    
    X(o(n)) = C(o(n), n+1)/C(o(n), o(n)) ;
    
    for i = (n-1):-1:1
        sum = 0;
        for j=i+1:n
            sum = sum + C(o(i), o(j))*X(o(j)) ;
        end
        X(o(i)) = (C(o(i), n+1) - sum)/C(o(i), o(i)) ; 
    end
end

function [m] = getMax(X)
    m = -1;
    for i = 1:size(X)
        if(m<X(i))
            m = X(i);
        end
    end
end