%% The following imports: http://becs.aalto.fi/en/research/bayes/ekfukf/install.html
addpath external
%OR use: mysetup from https://github.com/eruffaldi/matlabaddons
% mysetup('ekfukf');

%% Define the function and compute the gradient, then an evaluable matlab function
xx = sym(sym('x',[2,1]),'real');
f = @(x) [x(2); 2*(x(1))*(x(2)+x(1))];
fs = f(xx);
fs1 = fs;

Js = jacobian(fs); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function
x0 = [0.0054,1.01];

P0 = diag([0.01,0.2]);
[r,e]= compareapprox(f,J,x0,P0,500,[],'kl');
Hm = hessianmax(fs,xx,x0);

e.Hm = Hm;

rr{1} = e;


%% Define the function and compute the gradient, then an evaluable matlab function
xx = sym(sym('x',[2,1]),'real');
f = @(x) [x(2)^3; 2*(x(1))*(x(2)+x(1))];
fs = f(xx);
fs2 = fs;
Js = jacobian(fs); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function
x0 = [0.0054,1.01];

P0 = diag([0.01,0.2]);
[r,e]= compareapprox(f,J,x0,P0,500,[],'kl');
Hm = hessianmax(fs,xx,x0);

e.Hm = Hm;

rr{2} = e;

%% Define the function and compute the gradient, then an evaluable matlab function
xx = sym(sym('x',[2,1]),'real');
f = @(x) [exp(x(2)^2); 200*(x(1)^2)*(x(2)+x(1))];
fs = f(xx);
fs3 = fs;
Js = jacobian(fs); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function
x0 = [0.0054,1.01];

P0 = diag([0.01,0.2]);
[r,e]= compareapprox(f,J,x0,P0,500,[],'kl');
Hm = hessianmax(fs,xx,x0);

e.Hm = Hm;

rr{3} = e;
displayapprox(r,1,1);


%%
%
D = zeros(length(rr),3);
for I=1:length(rr)
    D(I,:) = [rr{I}.Hm,rr{I}.y.ds.linear(1),rr{I}.y.ds.ut(1)];
end
ds = dataset([D,{'H','EKF','UT'}],'Obs',{'f1','f2','f3'});
ds


