function [p, q] = getMaxIndex(A)
    
    m = -1 ;
    for i=1:size(A, 1)
        for j =(i+1):size(A, 1)
            if (m<abs(A(i, j)))
                m = abs(A(i, j));
                p = i;
                q = j;
            end
        end
    end
end