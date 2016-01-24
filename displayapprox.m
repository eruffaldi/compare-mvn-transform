function displayapprox(r,same,sameratio)
sx = r.x;
sg = r.sampling;
dims = 1:2;
figure;

s = r.sampling;
subplot(2,3,1);
% samples
scatter(s.x.pts(:,dims(1)),s.x.pts(:,dims(2)),[],'r');
hold on
% the truth
scatter(sx.mu(dims(1)),sx.mu(dims(2)),'k');
draw_ellipse(sx.mu(dims),sx.cov(dims,dims),'k');

% the estimated froms amples
scatter(s.x.mu(dims(1)),s.x.mu(dims(2)),'g');
draw_ellipse(s.x.mu(dims),s.x.cov(dims,dims),'g');
hold off
title('Input with sampling');


subplot(2,3,4);
% the result sampled
scatter(s.pts(:,dims(1)),s.pts(:,dims(2)),[],'r');
hold on
% the estimates from sampled
draw_ellipse(sg.mu(dims),sg.cov(dims,dims),'k');
scatter(sg.mu(dims(1)),sg.mu(dims(2)),'k');
hold off
title('Output with sampling');


s = r.lin;
subplot(2,3,2);
draw_ellipse(sx.mu(dims),sx.cov(dims,dims));
hold on
scatter(sx.mu(dims(1)),sx.mu(dims(2)),'g');
hold off
title('Input');
subplot(2,3,5);
scatter(s.mu(dims(1)),s.mu(dims(2)),'g');
hold on
draw_ellipse(s.mu(dims),s.cov(dims,dims));

% and then the good output
draw_ellipse(sg.mu(dims),sg.cov(dims,dims),'k');
scatter(sg.mu(dims(1)),sg.mu(dims(2)),'k');
hold off
title('Output Linearized');




% Unscented
s = r.ut;
subplot(2,3,3);
draw_ellipse(sx.mu(dims),sx.cov(dims,dims));
hold on
scatter(sx.mu(dims(1)),sx.mu(dims(2)),'g');
scatter(s.x.sigma(:,dims(1)),s.x.sigma(:,dims(2)),'m');
hold off
title('Input with Sigma points');

subplot(2,3,6);
draw_ellipse(s.mu(dims),s.cov(dims,dims));
hold on
scatter(s.mu(dims(1)),s.mu(dims(2)),'g');
scatter(s.sigma(:,dims(1)),s.sigma(:,dims(2)),'m');
draw_ellipse(sg.mu(dims),sg.cov(dims,dims),'k');
scatter(sg.mu(dims(1)),sg.mu(dims(2)),'k');

hold off
title('Output Unscented');

if same == 0
    if sameratio
for i=1:6
    subplot(2,3,i);
    axis equal
end
    end
    return
end

lx = zeros(6,2);
ly = zeros(6,2);
for i=1:6
    subplot(2,3,i);
    lx(i,:) = xlim;
    ly(i,:) = ylim;
end

mx = min(min(lx(1:3,:)));
Mx = max(max(lx(1:3,:)));
my = min(min(ly(1:3,:)));
My = max(max(ly(1:3,:)));
ns = 10;
sx = linspace(mx(1),Mx(1),ns);
sy = linspace(mx(1),Mx(1),ns);

s = 1.0;
for i=1:3
    subplot(2,3,i);
    xlim([mx(1) Mx(1)]*s);
    ylim([my(1) My(1)]*s);
    set(gca,'XTick',sx);
    set(gca,'XTick',sy);
    grid on
    if sameratio
        axis equal
    end

end

mx = min(min(lx(4:6,:)));
Mx = max(max(lx(4:6,:)));
my = min(min(ly(4:6,:)));
My = max(max(ly(4:6,:)));
sx = linspace(mx(1),Mx(1),ns);
sy = linspace(mx(1),Mx(1),ns);

for i=4:6
    subplot(2,3,i);
    xlim([mx(1) Mx(1)]*s);
    ylim([my(1) My(1)]*s);
    set(gca,'XTick',sx);
    set(gca,'XTick',sy);
    grid on
    if sameratio
        axis equal
    end

end


if 0==1
% all the x (1:3) should have the same axis
blocks = {1:3,4:6};
blocksaxis = [1,1;1,1];

    b = blocks{j};
    for i=1:length(b)
        subplot(2,3,b(i));
        lx(i,:) = xlim;
        ly(i,:) = ylim;
    end
    xl = [min(lx(:,1)) max(lx(:,2))];
    yl = [min(ly(:,2)) max(ly(:,2))];
    for i=1:length(b)
        subplot(2,3,b(i));
        if blocksaxis(j,1) == 1
            xlim(xl);
        end
        if blocksaxis(j,2) == 1
            ylim(yl);
        end
    end
end
end

