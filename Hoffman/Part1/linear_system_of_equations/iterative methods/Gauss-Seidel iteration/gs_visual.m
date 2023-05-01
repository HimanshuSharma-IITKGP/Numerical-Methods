
function [] = gs(A, b, conv)
    
    figure;
    hold on
    X = zeros(2, 1);
    plot(X(1), X(2), '.');

    X_new = zeros(2, 1);

    err = 1;

    while (err>conv)
        hold on
        X_new(1) = (b(1) - A(1, 2)*X(2))/A(1, 1);
        X_new(2) = (b(2) - A(2, 1)*X_new(1))/A(2, 2);

        plot([X_new(1), X(1)], [X_new(2), X(2)]);

        err = norm(X-X_new)/norm(X) ;
    end
end