%% The following imports: http://becs.aalto.fi/en/research/bayes/ekfukf/install.html
addpath external
%OR use: mysetup from https://github.com/eruffaldi/matlabaddons
% mysetup('ekfukf');

%% Define the function and compute the gradient, then an evaluable matlab function
rng(1)
xx = symmtx('x',[2,1],'real');
f = @(x) [x(2); 2*(x(1)).*(x(2)+x(1))];
fs = f(xx);

Js = jacobian(fs); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function
x0 = [0.0054,1.01];

a = pi/3;
P0 = diag([0.01,0.2]);
R = [cos(a),-sin(a); sin(a) cos(a)];
P = R*P0*R';
usesqrtm = @sqrtm; %cjo
usesqrtm = [];
[r,e]= compareapprox(f,J,x0,P,500,[],{'bhattacharyya_r','jeffreys'},usesqrtm);
close all
displayapprox(r,1,0);
Hm = hessianmax(fs,xx,x0)
e.x.bhattacharyya_r.ds
e.y.bhattacharyya_r.ds
e.y.jeffreys.ds

%% Another
xx = symmtx('x',[2,1],'real');
f = @(x) [exp(x(2).^2); 200*x(1).^2.*(x(2)+x(1))];
fs = f(xx);
Js = jacobian(fs); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function

% function, jacobian, point, number of samples for MC method, parameters
% for the Unscented Transformation
x0 = [0.0054,1.01];
usesqrtm = [];
[r,e] = compareapprox(f,J,x0,0.01^2*eye(2),500,[],{'bhattacharyya_r','jeffreys'},usesqrtm);
close all
displayapprox(r,1,1);
Hm = hessianmax(fs,xx,x0)
e.x.bhattacharyya_r.ds
e.y.bhattacharyya_r.ds
e.y.jeffreys.ds

%% for fx: n -> m, jacobian is: m by n, hessian: m by n by n
