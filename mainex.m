%% The following imports: http://becs.aalto.fi/en/research/bayes/ekfukf/install.html
addpath external
%OR use: mysetup from https://github.com/eruffaldi/matlabaddons
% mysetup('ekfukf');
rng(2);

%%
f1 = @(x) [x(2); 2*(x(1))*(x(2)+x(1))];
f2 = @(x) [x(2)^3; 2*(x(1))*(x(2)+x(1))];
f3 = @(x) [exp(x(2)^2); 200*(x(1)^2)*(x(2)+x(1))];

ff = {f1,f2,f3};

%%
xx = symmtx('x',[2,1],'real');

x0 = [0.0054,1.01];
P0 = diag([0.01,0.2]);
rr = {};
er = {};
for I=1:length(ff)
    f = ff{I};
    fs = f(xx);

    Js = jacobian(fs); % Jacobian
    J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function
    [r,e]= compareapprox(f,J,x0,P0,500,[],'kl');
    Hm = hessianmax(fs,xx,x0);

    e.Hm = Hm;
    rr{I} = r;
    er{I} = e;
end

%%
displayapprox(rr{3},1,1);


%%
%
D = zeros(length(er),3);
for I=1:length(er)
    D(I,:) = [er{I}.Hm,er{I}.y.ds.linear(1),er{I}.y.ds.ut(1)];
end
ds = dataset([D,{'H','EKF','UT'}],'Obs',{'f1','f2','f3'});
ds


