% Factorization LU implementation
function time2 = faktoryzacjaLU(A, b)

    op = matrixOperations;
    % 1 method
    %tic
    %[L, U] = op.lu_i(A);
    % x = op.backward_i(U, op.forward_i(L, b));
    % time = toc;
    % x = op.LU_frac_i(A, b);
    % res = op.decTables_i(op.multiplyTable_i(A,x),b);
    % op.norm_i(res)
    % disp(time);
    % 2 method
    tic
    x2 = op.LU_frac_i(A, b);
    time2 = toc;
    disp(time2)
end