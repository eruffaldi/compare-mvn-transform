%TODO: automate the handling of X/Y distributions that are infinite
%TODO: normalize the height of the gaussians (each) providing only the
%shape
mu = [0,0.2]';

rx = [-1:0.01:1];
ry = [-1:0.01:1];
SQ = [0.1,0.2; 0,0.4];
S = SQ*SQ';

%% Linear Transformation
theta = pi/3;
A = [cos(theta),-sin(theta); sin(theta),cos(theta)];
Q = zeros(2);
S2 = A*S*A';
mu2 = A*mu;

aa = prepare3plot(0.55);

% activate subplot axis, clear and hold on
set(gcf,'CurrentAxes',aa.c); cla; hold on;
title('Example of Rotation of a 2D Gaussian');

draw_ellipse(mu(1:2),S(1:2,1:2),'k');
scatter(mu(1),mu(2),[],'r+');

draw_ellipse(mu2(1:2),S2(1:2,1:2),'m');
scatter(mu2(1),mu2(2),[],'g+');

xlim([rx(1),rx(end)]);
ylim([ry(1),ry(end)]);
aa.cpost();

%left panel is y distribution: x is on the right, y bottom
set(gcf,'CurrentAxes',aa.y); cla; hold on;

q2 = normpdf(ry,mu(2),S(2));
plot(ry,q2)
q2 = normpdf(ry,mu2(2),S2(2));
plot(ry,q2)
aa.ypost();

set(gcf,'CurrentAxes',aa.x); cla; hold on;

q1 = normpdf(rx,mu(1),S(1));
plot(rx,q1)
q1 = normpdf(rx,mu2(1),S2(1));
plot(rx,q1)
aa.xpost();

%%
xo = 0.2;
mu2 = mu;
S2 = S;
mu2(2) = mu2(2) + S2(1,2)*inv(S2(2,2))*(xo-mu2(1));
mu2(1) = xo;
S2(2,2) = S2(2,2) - S2(1,2)*inv(S2(2,2))*S2(2,1);
S2(1,1) = 1e-10;
S2(2,1) = 0;
S2(1,2) = 0;

aa = prepare3plot(0.55);

% activate subplot axis, clear and hold on
set(gcf,'CurrentAxes',aa.c); cla; hold on;
title('Example of Conditioning of a 2D Gaussian');

draw_ellipse(mu(1:2),S(1:2,1:2),'k');
scatter(mu(1),mu(2),[],'r+');

draw_ellipse(mu2(1:2),S2(1:2,1:2),'m');
scatter(mu2(1),mu2(2),[],'g+');

xlim([rx(1),rx(end)]);
ylim([ry(1),ry(end)]);
aa.cpost();

%left panel is y distribution: x is on the right, y bottom
set(gcf,'CurrentAxes',aa.y); cla; hold on;

q21 = normpdf(ry,mu(2),S(2,2));
plot(ry,q21,'k')
if(S2(2,2) > 1e-5)
    q22 = normpdf(rx,mu2(2),S2(2,2));
    plot(ry,q22,'m')
end
aa.ypost();

set(gcf,'CurrentAxes',aa.x); cla; hold on;

q11 = normpdf(rx,mu(1),S(1,1));
plot(rx,q11,'k')
if(S2(1,1) > 1e-5)
    q12 = normpdf(rx,mu2(1),S2(1,1));
    plot(rx,q12,'m')
end
aa.xpost();
