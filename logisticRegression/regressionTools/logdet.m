function x = logdet(A)
% LOGDET - computes the log-determinant of a matrix A
%
% x = logdet(A);
%
% This is faster and more stable than using log(det(A))
%
% Input:
%     A NxN - A must be sqaure, positive semi-definite

%[R, p] = chol(A);
%x = 2*sum(log(diag(R)));
%x = 2*sum(log(diag(chol(A))));
 x = 2*sum(log(det(A)));
