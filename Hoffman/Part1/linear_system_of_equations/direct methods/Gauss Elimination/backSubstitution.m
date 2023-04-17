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
