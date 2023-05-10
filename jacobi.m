% Jacobi implementation
% png - if empty = "" then without graph else name.type
% A - matrix
% b - vector
% N - matrix size
function time = jacobi(A, b, N, png)
    op = matrixOperations;
    x = op.ones_i(N);

    L = op.tril_i(A,-1);
    U = op.triu_i(A,1);
    D = op.diag_i(A);

    calc1 = op.multiply_i(-op.reverse_diagonal(D), op.sum_i(L,U)); % D\(L + U)
    calc2 = op.multiplyTable_i(op.reverse_diagonal(D),b); % D\b

    iterations = 0;

    tic
    while true
        iterations = iterations + 1;
        x = op.sum_i(op.multiplyTable_i(calc1,x), calc2); % calc1*x + calc2
        res = op.decTables_i(op.multiplyTable_i(A,x),b); % Ax - b
        resXY(iterations) = op.norm_i(res);
        if resXY(iterations) <= 10^(-9) || isnan(resXY(iterations)) || resXY(iterations) >= 10^(100)
            break;
        end
    end
    time = toc;

    if png ~= ""     
        semilogy(resXY);
        title('Rozwiązanie iteracyjne metodą Jacobiego');
        xlabel('Nr iteracji');
        ylabel('Residuum');
        saveas(gcf, png);
    end
end