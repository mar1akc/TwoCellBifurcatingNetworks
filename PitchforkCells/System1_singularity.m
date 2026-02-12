function System1_singularity()
mu = 0.2;
N = 200;
x = linspace(-mu,2,N);
xmax = 2;
xmin = -2;
ymax = 2;
ymin = -2;
lam = sqrt(4*(mu+x).^3/(27*mu));
figure(1); clf; hold on; grid;
lwt = 5;
fsz = 20;
plot(x,lam,'Color','k','Linewidth',lwt)
plot(x,-lam,'Color','k','Linewidth',lwt)
set(gca,'Fontsize',fsz);

    verts_pink = [[x,xmax,xmax,x(end:-1:2)]',[lam,ymax,ymin,-lam(end:-1:2)]'];
    Nverts_pink = size(verts_pink,1);

    patch('Vertices', verts_pink, 'Faces', (1:Nverts_pink), ...
      'FaceColor', [0,0,0.5], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');

        verts_pink = [[x,xmin,xmin,x(end:-1:2)]',[lam,ymax,ymin,-lam(end:-1:2)]'];
    Nverts_pink = size(verts_pink,1);

    patch('Vertices', verts_pink, 'Faces', (1:Nverts_pink), ...
      'FaceColor', [0,0.5,0], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');
axis([xmin,xmax,ymin,ymax])

%%

e1 = 0.5;
e2 = (27*mu/4)^(1/3) - mu;
e3 = 2*e2 - e1;

plot([e1-0.3,e3 + 0.3],[1,1],'Linewidth',2,'Linestyle','--','Color','k');
plot([e1,e2,e3],[1,1,1],'.','Markersize',50,'Markeredgecolor','k');
plot([e1,e2,e3],[1,1,1],'.','Markersize',30,'Markeredgecolor','y');

col_stab = [0,0.7,0];
col_unst = [0,0,1];
ln_col = [0.8500 0.3250 0.0980];
saveas(gcf,'System1_singularity','epsc');
%%
figure(2);
fsz = 30;
clf; hold on; grid;


ep = e1;

rootfun = @(y)(-y.^3+(mu+ep)*y)/sqrt(mu);
y1 = linspace(ymin,-sqrt((mu + ep)/3),N);
y2 = linspace(-sqrt((mu + ep)/3),sqrt((mu + ep)/3),N);
y3 = linspace(sqrt((mu + ep)/3),ymax,N);
plot(rootfun(y1),y1,'Color',col_stab,'Linewidth',5);
plot(rootfun(y3),y3,'Color',col_stab,'Linewidth',5);
plot(rootfun(y2),y2,'Color',col_unst,'Linewidth',5,'LineStyle','--');
set(gca,'Fontsize',fsz);
axis([xmin,xmax,ymin,ymax])
plot([1,1],[ymin,ymax],'color',ln_col,'Linewidth',2);

fun = @(y)-y.^3+(mu+ep)*y - sqrt(mu);
y0 = fzero(fun,1);
plot(1,y0,'.','Markersize',50,'Color',ln_col);
saveas(gcf,'System1_singularity_A','epsc');

%%
figure(3);
clf; hold on; grid;

ep = e2;

rootfun = @(y)(-y.^3+(mu+ep)*y)/sqrt(mu);
y1 = linspace(ymin,-sqrt((mu + ep)/3),N);
y2 = linspace(-sqrt((mu + ep)/3),sqrt((mu + ep)/3),N);
y3 = linspace(sqrt((mu + ep)/3),ymax,N);
plot(rootfun(y1),y1,'Color',col_stab,'Linewidth',5);
plot(rootfun(y3),y3,'Color',col_stab,'Linewidth',5);
plot(rootfun(y2),y2,'Color',col_unst,'Linewidth',5,'LineStyle','--');
set(gca,'Fontsize',fsz);
axis([xmin,xmax,ymin,ymax])

plot([1,1],[ymin,ymax],'color',ln_col,'Linewidth',2);

fun = @(y)-y.^3+(mu+ep)*y - sqrt(mu);
plot(1,sqrt((mu+ep)/3),'.','Markersize',50,'Color',ln_col);

y0 = fzero(fun,-1);
plot(1,y0,'.','Markersize',50,'Color',ln_col);
saveas(gcf,'System1_singularity_B','epsc');
%%

figure(4);
clf; hold on; grid;

ep = e3;

rootfun = @(y)(-y.^3+(mu+ep)*y)/sqrt(mu);
y1 = linspace(ymin,-sqrt((mu + ep)/3),N);
y2 = linspace(-sqrt((mu + ep)/3),sqrt((mu + ep)/3),N);
y3 = linspace(sqrt((mu + ep)/3),ymax,N);
plot(rootfun(y1),y1,'Color',col_stab,'Linewidth',5);
plot(rootfun(y3),y3,'Color',col_stab,'Linewidth',5);
plot(rootfun(y2),y2,'Color',col_unst,'Linewidth',5,'LineStyle','--');
set(gca,'Fontsize',fsz);
axis([xmin,xmax,ymin,ymax])

plot([1,1],[ymin,ymax],'color',ln_col,'Linewidth',2);

fun = @(y)-y.^3+(mu+ep)*y - sqrt(mu);
y0 = fzero(fun,1);
plot(1,y0,'.','Markersize',50,'Color',ln_col);

y0 = fzero(fun,-1);
plot(1,y0,'.','Markersize',50,'Color',ln_col);

y0 = fzero(fun,0);
plot(1,y0,'.','Markersize',50,'Color',ln_col);
saveas(gcf,'System1_singularity_C','epsc');

figure(5); clf; hold on; grid;
y = linspace(ymin,ymax,N);
plot(y.^3,y,'Color',col_stab,'Linewidth',5);
set(gca,'Fontsize',fsz);
axis([xmin,xmax,ymin,ymax])
saveas(gcf,'System1_singularity_D','epsc');

figure(5); clf; hold on; grid;
y = linspace(ymin,ymax,N);
plot(y.^3,y,'Color',col_stab,'Linewidth',5);
set(gca,'Fontsize',fsz);
axis off
axis([xmin,xmax,ymin,ymax])
saveas(gcf,'System1_singularity_E','epsc');

end