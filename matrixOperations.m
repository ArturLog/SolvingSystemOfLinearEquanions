classdef matrixOperations
    methods
    
    % Creator for vector 
    function b = init_b(~, f, N)
        for i=1:N
            b(i) = sin(i*(f+1));
        end
    end

    % Creator for matrix
    function A = initMatrix(~, a1, a2, a3, N)
        for i=1:N
            A(i,i) = a1;
            if i < N
                A(i+1,i) = a2;
                A(i,i+1) = a2;
            end
            if i < N-1
                A(i+2,i) = a3;
                A(i,i+2) = a3;
            end
        end
    end
    
    % Creator for diagonal matrix
    function D = diag_i(obj, A)
        [numRows, ~] = size(A);
        obj.checkSquareMatrix(A);
        for i = 1:numRows
            D(i,i) = A(i,i);
        end
    end

    % Creator for lower triangle
    function U = triu_i(obj, A, n)
        [numRows, numCols] = size(A);
        obj.checkSquareMatrix(A);
        for i = 1:numRows
            for j = i+n:numCols
                U(i,j) = A(i,j);
            end
        end
        if n ~= 0
            U(numRows, numCols) = 0;
        end
    end

    % Creator for upper triangle
    function L = tril_i(obj, A, n)
        [numRows, numCols] = size(A);
        obj.checkSquareMatrix(A);
        for i = 1:numRows
            for j = 1:(i+n)
                L(i, j) = A(i, j);
            end
        end
        if n ~= 0
            L(numRows, numCols) = 0;
        end
    end

    % Creator for unit vector
    function x = ones_i(~, n)
        for i = 1:n
            x(i) = 1;
        end
    end

    % Creator for unit diagonal matrix
    function A = diag_ones_i(~, n)
        for i = 1:n
            A(i,i) = 1;
        end
    end

    % Function for summarising two matrix
    function C = sum_i(~, A, B)
        [numRows, numCols] = size(A);
        for i = 1:numRows
            for j = 1:numCols
                C(i,j) = A(i,j) + B(i,j);
            end
        end
    end

    % Function for substraction two matrix
    function C = dec_i(~, A, B)
        [numRows, numCols] = size(A);
        for i = 1:numRows
            for j = 1:numCols
                C(i,j) = A(i,j) - B(i,j);
            end
        end
    end

    % Function for substraction two vectors
    function c = decTables_i(~, a, b)
        for i = 1:length(a)
            c(i) = a(i) - b(i);
        end
    end

    % Function for multiply two matrix
    function C = multiply_i(~, A, B)
        [numRowsA, numColsA] = size(A);
        [numRowsB, ~] = size(B);
        % axb * cxd = bxc
        for i = 1:numColsA % rows
            for j = 1:numRowsB % cols
                tmp = 0;
                for k = 1:numRowsA % numRowsA == numColsB
                    tmp = tmp + (A(k,j) * B(i,k));
                end
                C(i,j) = tmp;
            end
        end
    end

    % Function for multiply matrix by vector
    function c = multiplyTable_i(~, A, b)
        [numRowsA, numColsA] = size(A);
        for i = 1:numRowsA
            tmp = 0;
            for j = 1:numColsA
                tmp = tmp + (A(i,j) * b(j));
            end
            c(i) = tmp;
        end
    end

    % Function for multiply matrix by const
    function C = multiplyConst_i(~, A, b)
        [numRowsA, numColsA] = size(A);
        for i = 1:numRowsA
            for j = 1:numColsA
                C(i,j) = A(i,j)*b;
            end
        end
    end

    % Function for reverse diagonal matrix
    function D = reverse_diagonal(obj, A)
        [numRowsA, ~] = size(A);
        obj.checkSquareMatrix(A);
        D = A;
        for i = 1:numRowsA
            D(i,i) = 1/A(i,i);
        end
    end
    
    % Function for forward substitution
    % A - square matrix
    % b - vector
    % return vector
    function y = forward_i(obj, A, b)
        [numRowsA, ~] = size(A);
        obj.checkSquareMatrix(A);
        L = A; % obj.tril_i(A, 0);
        y(1) = b(1)/L(1,1);
        for i = 2:numRowsA
            obj.checkSingular(L,i);
            sum = b(i);
            for j = 1:(i-1)
                sum = sum - L(i,j)*y(j);
            end
            y(i) = sum/L(i,i);
        end
    end

    % Function for backward substitution
    % A - square matrix
    % b - vector
    % return vector
    function x = backward_i(obj, A, b) 
       [numRowsA, ~] = size(A);
       n = numRowsA;
       U = obj.triu_i(A, 0);
       obj.checkSquareMatrix(A);
       for i = n:-1:1
            obj.checkSingular(U,i);
            x(i) = b(i) / U(i,i);
            for j = i-1:-1:1
                b(j) = b(j) - U(j,i) * x(i);
            end
        end
    end

    % Function for factorization
    % A - square matrix
    % return
    % L - lower triangle matrix
    % U - upper triangle matrix
    function [L, U] = lu_i(obj, A)
        [numRowsA, numColsA] = size(A);
        obj.checkSquareMatrix(A);
        L = obj.diag_ones_i(numRowsA);
        U = A;
        for k = 1:numRowsA-1
            for j = k+1:numRowsA
                L(j,k) = U(j,k) / U(k,k);
                for i = k:numColsA
                    U(j,i) = U(j,i) - (L(j,k)*U(k,i));
                end
            end
            obj.checkSingular(A, k);
        end
    end

    % Function for factorization ( another implementation )
    % A - square matrix
    % return vector
    function x = LU_frac_i(obj, A, b)
        [numRowsA, numColsA] = size(A);
        [L, U] = obj.lu_i(A);
        % A*x = b
        % L*y = b
        % U*x = y
        y(1) = b(1);
        for i = 2:numRowsA
            sum = b(i);
            for j = 1:(i-1)
                sum = sum - L(i,j)*y(j);
            end
            y(i) = sum;
        end
        
        x(numRowsA) = y(numRowsA)/U(numRowsA,numColsA); 
        for i = (numRowsA-1):-1:1
            sum = y(i);
            for j = numColsA:-1:(i+1)
                sum = sum - U(i,j)*x(j);
            end
            x(i) = sum/U(i,i);
        end
    end

    % Function for normalize
    % x - vector
    % return const
    function c = norm_i(~, x)
        sum = 0;
        for i = 1:length(x)
            sum = sum + (x(i)*x(i));
        end
        c = sqrt(sum);
    end

    % Error checking
    function checkSquareMatrix(~, A)
        [n, m] = size(A);
        if n ~= m
            error("Matrix is not Square")
        end
    end

    % Error checking
    function checkSingular(~, A, i)
        if A(i,i) == 0
            error("Matrix is singular")
        end
    end

    end
end