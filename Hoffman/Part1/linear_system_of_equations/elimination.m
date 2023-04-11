clear;
clc;

A = [3 2 105; 2 -3 103; 1 1 3];
b = [104; 98; 3];

eliminate(A, b)

function [C] = eliminate(A, b)
    C = [A, b];
    n = size(A, 1);

    for i = 1:n-1

                                % Pivoting %
%                 l = max(abs(C(i:n, i)));
%                 s = find(C(:, i) == l);
%                 [C(i, :), C(s, :)] = deal(C(s, :), C(i, :));
                                %  ***  %


                                % Scaled Pivoting
        max = -1;
        s = -1;
        for g=i:n
            temp = abs(C(g, i))/getMax(abs(C(g, :)'));
            if (temp>max)
                max = temp;
                s = g;
            end
        end

        [C(i, :), C(s, :)] = deal(C(s, :), C(i, :));

                                % *** %
        for j = i+1:n
                em =  C(j, i)/C(i, i) ;
                C(j, :) = C(j, :) - em*C(i, :);
        end
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