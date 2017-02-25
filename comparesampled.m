function w = comparesampled(X1,X2)
%Theory
%https://www.slac.stanford.edu/econf/C030908/papers/WEJT001.pdf
%https://pdfs.semanticscholar.org/aa3c/cba3ecab69cb51e53276a1d7597461104404.pdf
%https://indico.cern.ch/event/524705/contributions/2169813/attachments/1274007/1890263/IML2_2016.pdf
%http://www.math.tau.ac.il/~ruheller/Papers/1201.3522v3.pdf
%
%R Package: https://cran.r-project.org/web/packages/crossmatch/crossmatch.pdf
%Related
%Rosenbaum, P.R. (2005), An exact distribution-free test comparing two multivariate distributions
%based on adjacency, Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67,
%4, 515-530.
%
% Works by matching and with small size. 
% in R uses nonbimatch: http://biostat.mc.vanderbilt.edu/wiki/Main/MatchedRandomization
%


% Hellinger Distance
% https://it.mathworks.com/matlabcentral/fileexchange/36164-unscented-hellinger-distance-between-gmms
% M. Kristan, A. Leonardis, D. Skoaj, "Multivariate online Kernel Density
% Estimation", Pattern Recognition, 2011. 
% (url: http://vicos.fri.uni-lj.si/data/publications/KristanPR11.pdf)

%                    Z = pdist2(X1,X2);
 %                   Zm = min(Z);
  %                  w = mean(Zm);                    
w = 0;
