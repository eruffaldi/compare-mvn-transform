function h = draw_ellipse(x, c, outline_color)
%
% One single ellipse in x and with covariance c
%
% For multiple see

n = 40;					% resolution
radians = [0:(2*pi)/(n-1):2*pi,0];
unitC = [sin(radians); cos(radians)];
r = chol(c)';

if nargin < 3
  outline_color = 'g';
end

p = r*unitC;
p(1,:) = p(1,:) + x(1);
p(2,:) = p(2,:) + x(2);

h = plot(p(1,:),p(2,:),outline_color);


