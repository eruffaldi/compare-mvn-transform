function [r,e]= compareapprox(fx,Jx,x0,P0,nsamples,uc,mode)
%
% Compare Approximation
%
% fx = function
% J0  = Jacobian
% x0 = mean
% P0 = covariance
% nsamples = approximation for sampling
% sigmac = sigma point options (see ekfukf)
% mode = function for comparison: kl db bc fichet, default kl
%
% TODO: sigma point options
if isempty(mode)
    mode = 'kl';
end

if isempty(uc)
    uc.alpha = 1;
    uc.beta = 0;
    uc.kappa = 3-length(x0);
end

x0 = x0(:);

r = [];
s = [];
s.mu = x0;
s.cov = P0;
r.x = s;

% Sampling
R = mvnrnd2(x0,P0,nsamples);
Y = zeros(nsamples,length(x0));
for i=1:size(Y,1)
    Y(i,:) = fx(R(i,:));
end

s = [];
s.mu = mean(Y);
s.cov = cov(Y);
s.pts = Y;

s.x = [];
s.x.mu = mean(R);
s.x.cov = cov(R);
s.x.pts = R;

r.sampling = s;

% Gaussian Linearization
s = [];
J0 = Jx(x0);
s.mu = fx(x0);
s.cov = J0*P0*J0';

r.lin = s;

% Unscented
[WM,W,c] = ut_mweights(length(x0),uc.alpha,uc.beta,uc.kappa);

X = ut_sigmas(x0,P0,c);
s =[];
Y = zeros(length(x0),size(X,2));
for i=1:size(Y,2)
    Y(:,i) = fx(X(:,i));
end
s.sigma = Y';
s.mu = Y*WM;
s.cov = Y*W*Y';
s.covxy = X*W*Y';
s.uc = uc;
s.x = [];
s.x.sigma = X';
s.x.mu = X*WM;
s.x.cov = X*W*X';

r.ut = s;

if nargout > 1
    % result as table using 
    e = [];
    e.y = [];
    d = zeros(3);
    e.y.names = {'sampling','linear','ut'};
    values = {r.sampling,r.lin,r.ut};
    for I=1:3
        for J=I+1:3
            w = comparemvn(values{I}.mu,values{I}.cov,values{J}.mu,values{J}.cov,mode);
            d(I,J) = w;
            d(J,I) = w;
        end
    end
    
    e.y.ds = dataset([d,e.y.names],'ObsNames',e.y.names);
    e.y.d = d;
    
    e.x = [];
    e.x.names = {'ref','sampling','sigma'};
    d = zeros(3);
    values = {r.x,r.sampling.x,r.ut.x};
    for I=1:3
        for J=I+1:3
            w = comparemvn(values{I}.mu,values{I}.cov,values{J}.mu,values{J}.cov,mode);
            d(I,J) = w;
            d(J,I) = w;
        end
    end
    e.x.ds = dataset([d,e.x.names],'ObsNames',e.x.names);
    e.x.d = d;
end




