%% The following imports: http://becs.aalto.fi/en/research/bayes/ekfukf/install.html
mysetup('ekfukf');
% Externals:
% ut_mweights
% ut_sigmas
% ut_weights
% mvnrnd2

%% Define the function
xx = sym(sym('x',[2,1]),'real');
f = @(x) [x(2); 2*(x(1))*(x(2)+x(1))];
Js = jacobian(f(xx)); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function

a = pi/3;
P0 = diag([0.01,0.2]);
R = [cos(a),-sin(a); sin(a) cos(a)];
P = R*P0*R';
[r,e]= compareapprox(f,J,[0.0054,1.01],P,500,[],'kl');
close all
displayapprox(r,1,0);

e.y.ds


%% Define the function and compute the graident, then an evaluable matlab function
xx = sym(sym('x',[2,1]),'real');
f = @(x) [exp(x(2)^2); 200*x(1)^2*(x(2)+x(1))];

Js = jacobian(f(xx)); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function

% function, jacobian, point, number of samples for MC method, parameters
% for the Unscented Transformation
[r,e] = compareapprox(f,J,[0.0054,1.01],0.01^2*eye(2),500,'kl');
close all
displayapprox(r,1,1);

e.y.dclosclse