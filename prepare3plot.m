function aa = prepare3plot(mainsize,margins)

if nargin < 1
    mainsize = 0.55;
end
if nargin < 2
    margins = [0.05,0.05,0.05,0.05,0.1,0.1];
end

margininx = margins(5);
margininy = margins(6);
marginleft = margins(1);
marginbottom = margins(2);
marginright = margins(3);
margintop = margins(4);

X = [marginleft,0,margininx,mainsize,marginright];
X(2) = 1.0-sum(X);
Y = [marginbottom,0,margininy,mainsize,margintop];
Y(2) = 1.0-sum(Y);
cX = cumsum(X);
cY = cumsum(Y);
% vertical 0 is below
posC = [cX(3),cY(3),X(4),Y(4)];
posL = [cX(1),cY(3),X(2),Y(4)];
posB = [cX(3),cY(1),X(4),Y(2)];

clf
a = subplot('Position',posC);
b = subplot('Position',posL);


c = subplot('Position',posB);

aa.c = a;
aa.y = b;
aa.x = c;
aa.ypost = @ypost;
aa.xpost = @xpost;
aa.cpost = @cpost;

function ypost()
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.XTick = [];
view(-90,90)

function cpost()

function xpost()
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'right';
ax.XTick = [];
view(180,90)
