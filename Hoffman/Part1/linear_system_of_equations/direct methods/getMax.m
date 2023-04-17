
function [m] = getMax(X)
    m = -1;
    for i = 1:size(X)
        if(m<X(i))
            m = X(i);
        end
    end
end