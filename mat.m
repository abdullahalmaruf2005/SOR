% SOR Method for solving linear system Ax = b
A = [5, 2, -1; 1, 6, -3; 2, 1, 4];
b = [6; 4; 7];
x = [0; 0; 0]; % initial guess (প্রাথমিক অনুমান)
maxiter = 50; % সর্বোচ্চ iteration সংখ্যা
maxerr = input('Enter the error allowed: '); % error tolerance ইনপুট
w = input('Enter the relaxation parameter: '); % relaxation parameter ইনপুট

% Check for matrix A
if size(A,1) ~= size(A,2)
    disp('ERROR: Matrix A should be a square matrix');
    return;
end

% Check for matrix b
if size(b,1) ~= size(A,1) || size(b,2) ~= 1
    disp('ERROR: Input error..pls re-enter data');
    return;
end

% Check for initial guess
if size(x,1) ~= size(A,1) || size(x,2) ~= 1
    disp('ERROR: Pls check input');
    return;
end

% Decompose A = D + L + U
D=diag(diag(A)); % diagonal matrix(D)
L=tril(A)-D;% lower triangular part(L)
U=triu(A)-D;% upper triangular part(U)

% SOR iteration matrix and constant vector
H=inv(D + w*L) * ((1-w)*D - w*U); % iteration matrix তৈরি
C=w*inv(D+w*L)*b; % constant vector

% Check spectral radius
spectral_radius=max(abs(eig(H)));% spectral radius বের করা
fprintf('Spectral radius: %.6f\n',spectral_radius);

% Convergence condition
if spectral_radius >= 1
    disp('ERROR: Method will not converge (spectral radius >= 1)');
    disp('Try a relaxation parameter w between 0 and 2');
    return;
end

% Iterative solution
iter = 0;
err = inf;
X = x; %initial solution vector

fprintf('\nIteration\t x1\t\t x2\t\t x3\t\t Error\n');
fprintf('-----------------------------------------------------------\n');

while err > maxerr && iter < maxiter
    X_old = X;
    X = H * X_old + C; % নতুন iteration update
    err = norm(X - X_old, inf); % error calculation
    iter = iter + 1;
    
    fprintf('%d\t\t %.6f\t %.6f\t %.6f\t %.6f\n', iter, X(1), X(2), X(3), err);
end

% Display results
fprintf('-----------------------------------------------------------\n');
if err <= maxerr
    fprintf('CONVERGED! Solution after %d iterations:\n', iter);
    fprintf('x1 = %.6f\n', X(1));
    fprintf('x2 = %.6f\n', X(2));
    fprintf('x3 = %.6f\n', X(3));
else
    fprintf('NOT CONVERGED! Maximum iterations (%d) reached\n', maxiter);
    fprintf('Final approximation:\n');
    fprintf('x1 = %.6f\n', X(1));
    fprintf('x2 = %.6f\n', X(2));
    fprintf('x3 = %.6f\n', X(3));
end

% Compare with exact solution
exact = A\b; % direct solution
fprintf('\nComparison with exact solution:\n');
fprintf('Exact: x1 = %.6f, x2 = %.6f, x3 = %.6f\n', exact(1), exact(2), exact(3));
fprintf('Error: %.6f\n', norm(X - exact, inf)); % final error