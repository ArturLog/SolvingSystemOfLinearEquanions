clc
clear all
close all

op = matrixOperations;
% Ex. A
c = 6; 
d = 3; 
N = 9*100 + c*10 + d;
e = 6; 
a1 = 5 + e;
a2 = -1;
a3 = -1;
f = 8;
b = op.init_b(f, N);
A = op.initMatrix( a1,  a2, a3, N);

% Ex. B 

jacobi(A, b, N, 'exB_jacobi.png');
gauss_seidl(A, b, N, 'exB_gauss-seidl.png');

% Ex. C

a1 = 3;
B = op.initMatrix(a1, a2, a3, N);
jacobi(B, b, N, 'exC_jacobi.png');
gauss_seidl(B, b, N, 'exC_gauss-seidl.png');

% Ex. D

faktoryzacjaLU(B, b);

% Ex. E

N2 = [ 100, 500, 1000, 2000, 3000 ];
c = 6;
d = 3;
e = 6;
a1 = 5 + e;
a2 = -1;
a3 = -1;
f = 8;
for i = 1:length(N2)
    b2 = op.init_b(f, N2(i));
    A2 = op.initMatrix( a1,  a2, a3, N2(i));
    
    time_LU(i) = faktoryzacjaLU(A2, b2);
    time_J(i) = jacobi(A2, b2, N2(i), "");
    time_GS(i) = gauss_seidl(A2, b2, N2(i), "");
end

plot(N2, time_LU);
title('Czas wyznaczania metodą faktoryzacji LU');
xlabel('N');
ylabel('Czas [s]');
saveas(gcf, 'exE_LU.png');

plot(N2, time_J);
title('Czas wyznaczania metodą Jacobiego');
xlabel('N');
ylabel('Czas [s]');
saveas(gcf, 'exE_jacobi.png');

plot(N2, time_GS);
title('Czas wyznaczania metodą Gaussa-Seidla');
xlabel('N');
ylabel('Czas [s]');
saveas(gcf, 'exE_gauss-seidl.png');





