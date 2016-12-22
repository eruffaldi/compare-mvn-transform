function Y = svdsqrt(X)

[U,S,V] = svd(X);
R = V;
S = sqrt(S);
Y = R*S;
disp('svd')
Y*Y'-X
disp('cholcov')
C=cholcov(X)';
if isempty(C) == 0
    C*C'-X
end

