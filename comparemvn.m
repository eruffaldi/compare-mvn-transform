function [e,emu,eP] = comparemvn(m1,P1,m2,P2,mode,alpha)
% covariance comparison using various methods
%
% References
% Used: http://like.silk.to/studymemo/PropertiesOfMultivariateGaussianFunction.pdf
% alternative: SAS https://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_glimmix_a0000001447.htm
% http://www.subcortex.net/research/code/testing-for-differences-in-multidimensional-distributions
%
% MVN:
% - Kullback-Leibler divergence: general
% - Bhattacharyya distance: general
% - MAYBE Hellinger Distance: M. Kristan, A. Leonardis, D. Skoaj, "Multivariate online Kernel Density
% Estimation", Pattern Recognition, 2011. 
% - DONE Fr?chet distance: http://www.sciencedirect.com/science/article/pii/0047259X8290077X
%       HAS issues with chol(Sx*Sy)
% Helling Distance implementations: http://it.mathworks.com/matlabcentral/fileexchange/36164-unscented-hellinger-distance-between-gmms
% Theory for general (sampled) distributions: Kolmogorov-Smirnov test in multidimension:http://articles.adsabs.harvard.edu/full/1987MNRAS.225..155F
%
%
% Better: A Note on Metric Properties for Some Divergence Measures: The Gaussian Case
% Abou-Mustafa 
% JMLR: Workshop and Conference Proceedings 25:1?15, 2012 Asian Conference on Machine Learning
m1=m1(:);
m2=m2(:);
if nargin < 6
    alpha = 0.5;
end
u = m1-m2;
switch(mode)
    case {'jeffreys','skl'}
        % Symmetric KL aka Jeffreys divergence
        Psi = inv(P1)+inv(P2);
        e = 0.5*(u'*Psi*u)+0.5*trace(inv(P1)*P2+inv(P2)*P1-2*eye(size(P1)));
    case 'jeffreys_r'
        Psi = inv(P1)+inv(P2);
        l = eig(Psi);
        dr = sqrt(sum(log(l).^2));
        e = alpha*sqrt(u'*Psi*u)+(1-alpha)*dr;
    case 'hellinger'
        e = sqrt(2*(1-bhattacharyya_coeff(u,P1,P2)));
    case 'bhattacharyya_mean'
        Gamma = (P1+P2)/2;
        e = 1/8*(u'*inv(Gamma)*u);
    case 'bhattacharyya'
        % e = -log(bhattacharyya_coeff(u,P1,P2))
        Gamma = (P1+P2)/2;
        e = 1/8*(u'*inv(Gamma)*u)+1/2*log(isqrt(det(P1))*isqrt(det(P2))*det(Gamma));
    case 'bhattacharyya_r'
        Gamma = (P1+P2)/2;
        Psi = inv(P1)+inv(P2);      
        try
        l = eig(Psi);
        catch me
            l = Inf;
        end
        dr = sqrt(sum(log(l).^2));
        e = alpha*sqrt(u'*inv(Gamma)*u)+(1-alpha)*dr;
        
    case 'frichet'       
        emu = norm(m1-m2).^2;
        % The following is the
        try
        eP = trace(P1+P2-2*chol(P1*P2)); % SUBJECT to chol ISSUES, requires schol of ekfukf
        catch me
            eP = NaN;
        end
        e = emu + eP;
        
    %case 'kl' % not symmetric
    %    D = length(P1);
    %    nP1 = det(P1);
    %    nP2 = det(P2);
    %    e = 0.5*(trace(inv(P2)*P1)+u'*inv(P2)*u-log(nP1/nP2)-D);
%     case 'bc'
%         % 
%         Gamma = (P1+P2)/2;
%         nP1 = det(P1);
%         nP2 = det(P2);
%         nP = det(Gamma);
%         e = exp(-0.8*u'*inv(Gamma)*u)*sqrt(sqrt(nP1*nP2)/nP);
    otherwise
       error('Required: jeffreys (skl) jeffreys_r frichet bhattacharyya bhattacharyya_r hellinger');
end


function r = isqrt(x)

r = 1/sqrt(x);

function b = bhattacharyya_coeff(u,P1,P2)

Gamma = (P1+P2)/2;
nP1 = det(P1);
nP2 = det(P2);
nP = det(Gamma);

b = isqrt(nP)*pow(P1,1/4)*pow(P2,1/4)*exp(-1/8*u'*inv(Gamma)*u);



