function [e,emu,eP] = comparemvn(m1,P1,m2,P2,mode)
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

m1=m1(:);
m2=m2(:);
switch(mode)
    case 'frichet'
        
        emu = norm(m1-m2).^2;

        % The following is the
        try
        eP = trace(P1+P2-2*chol(P1*P2)); % SUBJECT to chol ISSUES, requires schol of ekfukf
        catch me
            eP = NaN;
        end
        e = emu + eP;
    case 'kl'
        D = length(P1);
        nP1 = det(P1);
        nP2 = det(P2);
        e = 0.5*(trace(inv(P2)*P1)+(m2-m1)'*inv(P2)*(m2-m1)-log(nP1/nP2)-D);
    case 'bc'
        P = (P1+P2)/2;
        nP1 = det(P1);
        nP2 = det(P2);
        nP = det(P);
        e = exp(-0.8*(m2-m1)'*inv(P)*(m2-m1))*sqrt(sqrt(nP1*nP2)/nP);
    case 'db'
        P = (P1+P2)/2;
        nP1 = det(P1);
        nP2 = det(P2);
        nP = det(P);
        e = 1/8*(m2-m1)'*inv(P)*(m2-m1)+0.5*log(nP/(sqrt(nP1*nP2)));
    otherwise
       error('Required: kl frichet bc db');
end

