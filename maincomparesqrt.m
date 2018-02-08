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
P0 = diag([0.01,0.3]);
R = [cos(a),-sin(a); sin(a) cos(a)];
P = R*P0*R';

usesqrtm = [];
[r0,e0]= compareapprox(f,J,x0,P,500,[],{'sampling','bhattacharyya_r','jeffreys','bhattacharyya_mean'},usesqrtm);
%
usesqrtm = @(x) svdsqrt(x); %(x) chol(eigv(x)); %cjo
[r,e]= compareapprox(f,J,x0,P,500,[],{'sampling','bhattacharyya_r','jeffreys','bhattacharyya_mean'},usesqrtm);
close all
subplot(1,2,1);
% Unscented
dims=1:2;
s = r0.ut;
sx = r0.x;
sg = r0.sampling;

draw_ellipse(sx.mu(dims),sx.cov(dims,dims),'k');
hold on
%h1 = scatter(sx.mu(dims(1)),sx.mu(dims(2)),'k');
h2 = scatter(s.x.sigma(:,dims(1)),s.x.sigma(:,dims(2)),'m');
h3 = scatter(s.x.mu(dims(1)),s.x.mu(dims(2)),'g');
draw_ellipse(s.x.mu(dims),s.x.cov(dims,dims),'g');
hold off
legend([h2,h3],{'Sigma Points','Est Sigma Point'});
axis equal
title('Input with Sigma points CHOL');
xlim([-0.9,0.9]);

subplot(1,2,2);
% Unscented
s = r.ut;
sx = r.x;
sg = r.sampling;

draw_ellipse(sx.mu(dims),sx.cov(dims,dims),'k');
hold on
%h1 = scatter(sx.mu(dims(1)),sx.mu(dims(2)),'k');
h2 = scatter(s.x.sigma(:,dims(1)),s.x.sigma(:,dims(2)),'m');
h3 = scatter(s.x.mu(dims(1)),s.x.mu(dims(2)),'g');
draw_ellipse(s.x.mu(dims),s.x.cov(dims,dims),'g');
hold off
legend([h2,h3],{'Sigma Points','Est Sigma Point'});
axis equal
title('Input with Sigma points SQRT');
xlim([-0.9,0.9]);

%%
set(gcf,'Position',[   440   550   414   248]);
mysetup('exportfig');
export_fig('comparesqrt.pdf','-transparent')


