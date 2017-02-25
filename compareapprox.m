function [r,e]= compareapprox(fx,Jx,x0,P0,nsamples,uc,mode,usesqrtm)
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
% TODO: add other quadrature examples
%
% Emanuele Ruffaldi Scuola Superiore Sant'Anna 2016
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
Ys = zeros(nsamples,length(x0));
for i=1:size(Ys,1)
    Ys(i,:) = fx(R(i,:));
end

s = [];
s.mu = mean(Ys);
s.cov = cov(Ys);
s.pts = Ys;

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

X = ut_sigmas(x0,P0,c,usesqrtm);
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

Yp = ut_sigmas(s.mu,s.cov,c,usesqrtm);

s.sigmapost = Yp';

r.ut = s;
r.f = fx;

if nargout > 1
    % NEW table:
    % - variable(x or y)    
    % - measure
    % - value
    ts = [];
    
    % result as table using 
    e = [];

    if iscell(mode) == 0
        modes = {mode};
    else
        modes = mode;
    end
    
    e.y = [];
    e.y.names = {'sampling','linear','ut'};
    e.x = [];
    e.x.names = {'ref','sampling','sigma'};
    
    ulinY = [];
    uutY = [];
    uutX = [];
    % for sampling modes
    for K=1:length(modes)
        if strcmp(mode,'sampling')
            ulinY = mvnrnd2(r.lin.mu,r.lin.cov,nsamples);
            uutY =  mvnrnd2(r.ut.mu,r.ut.cov,nsamples);
            uutX =  mvnrnd2(r.ut.x.mu,r.ut.x.cov,nsamples);
            break;
        end
    end
    
    for K=1:length(modes)
        mode = modes{K};
        
        d = zeros(3);
        values = {r.sampling,r.lin,r.ut};
        
        % sampling
        uss = {Ys, ulinY,uutY };
        
        for I=1:3
            for J=I+1:3
                if strcmp(mode,'sampling') == 0
                    w = comparemvn(values{I}.mu,values{I}.cov,values{J}.mu,values{J}.cov,mode);
                else
                    %Z = squareform(pdist(uss{I},uss{J}));
                    w = comparesampled(uss{I},uss{J});
                end
                d(I,J) = w;
                d(J,I) = w;
                ets = [];
                ets.mode = mode;
                ets.what = 'y';
                ets.var1 = e.y.names{I};
                ets.var2 = e.y.names{J};
                ets.value = w;
                if isempty(ts)
                    ts = ets;
                else                    
                    ts = [ts; ets];
                end
            end
        end

        e.y.(mode).ds = dataset([d,e.y.names],'ObsNames',e.y.names);
        e.y.(mode).d = d;

        d = zeros(3);
        values = {r.x,r.sampling.x,r.ut.x};
        uss = {R,R,uutX};
        
        for I=1:3
            for J=I+1:3
                if strcmp(mode,'sampling') == 0
                    w = comparemvn(values{I}.mu,values{I}.cov,values{J}.mu,values{J}.cov,mode);
                else
                    w = comparesampled(uss{I},uss{J});
                end
                d(I,J) = w;
                d(J,I) = w;
                ets = [];
                ets.mode = mode;
                ets.what = 'x';
                ets.var1 = e.x.names{I};
                ets.var2 = e.x.names{J};
                ets.value = w;
                if isempty(ts)
                    ts = ets;
                else                    
                    ts = [ts; ets];
                end
            end
        end
        e.x.(mode).ds = dataset([d,e.x.names],'ObsNames',e.x.names);
        e.x.(mode).d = d;
    end
    e.table = struct2table(ts);
end




