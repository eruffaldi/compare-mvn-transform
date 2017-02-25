addpath external

%%
rng(1)
xx = symmtx('x',[2,1],'real');
%f = @(x) [x(2); 2*(x(1)).*(x(2)+x(1))];
%f = @(x) [exp(x(2).^2); 200*x(1).^2.*(x(2)+x(1))];
f = @(x) [cos(x(2)^3*20); 2*(cos(x(1)*1000))*(x(2)+x(1))];
f = @(x) [x(1), cos(x(2)*1000)];
%f = @(x) [cos(x(1)*10),x(2)];
fs = f(xx);


Js = jacobian(fs); % Jacobian
J = matlabFunction(Js,'Vars',{xx}); % Vars is required to make 1 input function

a = pi/2;
P0 = diag([0.01,0.1]);
R = [cos(a),-sin(a); sin(a) cos(a)];
P = R*P0*R';

%{'sampling','bhattacharyya_r','jeffreys','bhattacharyya_mean'}
mode = 'jeffreys';
[X,Y] = meshgrid(-1:0.1:1,-1:0.1:1);

nsamples = 500;

Z = zeros([size(X),2]); 
W = zeros([size(X),3]);
for I=1:size(X,1)
    for K=1:size(X,2)
        x0 = [X(I,K),Y(I,K)];
        Z(I,K,:) = f(x0);

        %usesqrtm = [];
        %[r0,e0]= compareapprox(f,J,x0,P,500,[],mode,usesqrtm);

        usesqrtm = @(x) svdsqrt(x); %(x) chol(eigv(x)); %cjo
        
        [r,e]= compareapprox(f,J,x0,P,nsamples,[],mode,usesqrtm);
        d =  e.y.(mode).d;
        W(I,K,:) = [d(1,2),d(1,3),d(2,3)]; % s-l s-u l-u
    end
end
wm = squeeze(min(min(W,[],1),[],2));
wM = squeeze(max(max(W,[],1),[],2));

q = {'S-L','S-U','L-U'};
figure(1)
for I=1:3
subplot(3,1,I);
colormap(hot)
%imshow(W(:,:,I),[wm(I),wM(I)]);
WI= W(:,:,I);
scatter(X(:),Y(:),[],WI(:),'filled');
colorbar
% flipped due to the imshow
xlabel('x');
ylabel('y');
title(sprintf('%s %s %f-%f',q{I},mode,wm(I),wM(I)));
end

figure(2)
subplot(2,1,1);
surf(X,Y,Z(:,:,1));
title('X output');
xlabel('x');
ylabel('y');
subplot(2,1,2);
surf(X,Y,Z(:,:,2));
title('Y output');
xlabel('x');
ylabel('y');
figure(3)
XY = [X(:),Y(:)];
Q = floor(size(X)/2);
x0 = [X(Q(1),Q(2)),Y(Q(1),Q(2))];
Zg = reshape(mvnpdf(XY,repmat(x0,size(XY,1),1),P),size(X));
surf(X,Y,Zg)

