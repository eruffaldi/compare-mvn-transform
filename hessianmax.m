function [Hm,H] = hessianmax(fs,xx,x0)

cxx =  num2cell(xx)';
cx0 = num2cell(x0);

% m n n 
H = zeros(length(fs),length(xx),length(xx));
for I=1:length(fs)
    Hsi = hessian(fs(I),xx); % n n 
    Hi =double(subs(Hsi,cxx,cx0));
    H(I,:,:) = Hi;
end
Hm =max(abs(H(:)));

%Hs2 = hessian(fs(2),xx);
%H2 = double(subs(Hs2,{xx(1),xx(2)},{x0(1),x0(2)}));
%H = zeros(2,2,2);
%H(1,:,:) = H1;
%H(2,:,:) = H2;
