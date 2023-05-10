% Gauus-Seidl implementation
% png - if empty = "" then without graph else name.type
% A - matrix
% b - vector
% N - matrix size
function time = gauss_seidl(A, b, N, png)
    op = matrixOperations;
    x = op.ones_i(N);

    L = op.tril_i(A,-1);
    U = op.triu_i(A, 1);
    D = op.diag_i(A);
    
    calc1 = -(op.sum_i(D,L)); % -(D + L)
    calc2 = op.forward_i(op.sum_i(D,L), b); % (D + L)\b
    iterations = 0;

    tic
    while true
        iterations = iterations + 1;
        calc3 = op.multiplyTable_i(U,x); % (U*x)
        x = op.sum_i(op.forward_i(calc1, calc3), calc2); % calc1\calc3 + calc2
        res = op.decTables_i(op.multiplyTable_i(A,x),b); % Ax - b
        resXY(iterations) = op.norm_i(res);
        if resXY(iterations) <= 10^(-9) || isnan(resXY(iterations)) || resXY(iterations) >= 10^(100)
            break;
        end
    end
    time = toc;
    
    if png ~= ""
        semilogy(resXY);
        title('Rozwiązanie iteracyjne metodą Gauusa-Seidla');
        xlabel('Nr iteracji');
        ylabel('Residuum');
        saveas(gcf, png);
    end
end