function [X,A] = ut_sigmas(M,P,c,sqrtfx)


% [A,err] = cholcov(P);
% 
% if isempty(A)
%     [n,m] = size(P);
%     disp('Semi')
%         % Can get factors of the form Sigma==T'*T using the eigenvalue
%         % decomposition of a symmetric matrix, so long as the matrix
%         % is positive semi-definite.
%         [U,D] = eig(full((P+P')/2));
% 
%         % Pick eigenvector direction so max abs coordinate is positive
%         [ignore,maxind] = max(abs(U),[],1);
%         negloc = (U(maxind + (0:n:(m-1)*n)) < 0);
%         U(:,negloc) = -U(:,negloc);
% 
%         D = diag(D);
%         tol = eps(max(D)) * length(D);
%         t = (abs(D) > tol);
%         D = D(t);
%         err = sum(D<0); % number of negative eigenvalues
% 
%         if (err ==0)
%             A = diag(sqrt(D)) * U(:,t)';
%         else
%             error('Super Failure')
%             A = [];
%         end
% end

% [A err] = cholcov(P, 0);
% if (err ~= 0)
%     % the covariance matrix is not positive definite!
%     [v d] = eig(P);
% 
%     % set any of the eigenvalues that are <= 0 to some small positive value
%     for n = 1:size(d,1)
%         if (d(n, n) <= ZERO)
%             d(n, n) = EPS;
%         end
%     end
%     % recompose the covariance matrix, now it should be positive definite.
%     Pn = v*d*v';
% 
%     [A err] = cholcov(Pn, 0);
%     if (err ~= 0)
%         error('ut_sigma failure');
%     end
% end
% 
% assert(isreal(A),'ut_sigmas should give real');

if nargin == 3
    sqrtfx = [];
end
if isempty(sqrtfx) == 0
    A = sqrtfx(P);
else
    A = chol(P);
end

A = A';

X = [zeros(size(M)) A -A];
X = sqrt(c)*X + repmat(M,1,size(X,2));

