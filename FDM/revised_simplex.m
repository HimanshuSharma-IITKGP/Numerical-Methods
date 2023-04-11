c = [2, 1];
A = [6 4; 6 1];
b = [6; 3];



num_constraints = size(A, 1);
num_dv = size(A, 2);
num_slacks = num_constraints;


augumented_A = [A, eye(num_slacks)];

c_b = zeros(1, num_slacks);
B = eye(num_slacks);

X_b_index = [] ;

while(true)


    X_b = B\b;
    Z = c_b*X_b ;

    A_1 = c_b*inv(B)*A - c;
    A_2 = c_b*inv(B);

    a = inf;
    ind = -1;

    for i = 1:size(A_1, 2)
        if(A_1(i)<inf)
            a = A_1(i);
            ind = i;
        end
    end

    for i = 1:size(A_2, 2)
        if(A_2(i)<inf)
            a = A_2(i);
            ind = i + num_dv;
        end
    end

    leaving_ind = -1;
    col = zeros(num_constraints , 1);

    if(ind <= num_dv)
        A_3 = inv(B)*A;
        col = A_3(:, ind);
    else
        A_4 = inv(B);
        col = A_4(:, ind-num_dv);
    end



end
