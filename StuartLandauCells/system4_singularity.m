function system4_singularity()
close all
% Plots bifurcation and hysteresis surfaces in (epsilon,gamma,lambda)-space
% x-axis : epsilon
% y-axis : gamma
% z-axis : lambda
save_count = 0;

fsz  = 20;
lwt  = 1;
lwt1 = 3;
colH = [0.2 0.4 1];
colB = [0,0,0];

figure(1); clf; hold on; grid on;

gmin = -1.5;
gmax =  1.5;
emax =  1.5;
lmin =  1e-10;
lmax =  1.5;

mu = 0.2;

Ngam = 101;
Nlam = 101;
Ne   = 101;

gam = linspace(gmin,gmax,Ngam);
lam = linspace(lmin,lmax,Nlam);
B_ep = linspace(-mu,emax,Ne);
B_lam = sqrt(4*(B_ep+mu).^3/(27*mu));

ep1 = zeros(Ngam,Nlam);
ep2 = zeros(Ngam,Nlam);

%% Compute hysteresis data
for jg = 1:Ngam
    for jl = 1:Nlam
        [~,~,A] = Hfun(mu,gam(jg),lam(jl));
        ep1(jg,jl) = A(1) - mu;
        ep2(jg,jl) = A(2) - mu;
    end
end

%% Build grids
[Gam,Lam] = meshgrid(gam,lam);
Ep1 = ep1.';
Ep2 = ep2.';

%% Hysteresis surfaces
surf(Ep1,Gam,Lam,...
    'FaceColor',colH,...
    'FaceAlpha',0.6,...
    'EdgeColor','none');

surf(Ep2,Gam,Lam,...
    'FaceColor',colH,...
    'FaceAlpha',0.6,...
    'EdgeColor','none');

%% Bifurcation surface
[GamB,EpB] = meshgrid(gam,B_ep);
LamB = repmat(B_lam.',1,Ngam);

surf(EpB,GamB,LamB,...
    'FaceColor',[0 0 0],...
    'FaceAlpha',0.4,...
    'EdgeColor','none');

%% Gamma planes
eps_plane = linspace(-mu,emax,50);
lam_plane = linspace(lmin,lmax,50);
[EpP,LamP] = meshgrid(eps_plane,lam_plane);

gamma_planes = [0 1 -1];
a = 0.1;
for g0 = gamma_planes
    GamP = g0*ones(size(EpP));
    surf(EpP,GamP,LamP,...
        'FaceColor','r',...
        'FaceAlpha',a,...
        'EdgeColor','none');
end

% %% Curves on gamma planes
% for g0 = gamma_planes
% 
%     [~,jg] = min(abs(gam-g0));
% 
%     % Hysteresis curves
%     plot3(ep1(jg,:), g0*ones(1,Nlam), lam,...
%         'k','LineWidth',lwt1);
% 
%     plot3(ep2(jg,:), g0*ones(1,Nlam), lam,...
%         'k','LineWidth',lwt1);
% 
%     % Bifurcation curve
%     plot3(B_ep, g0*ones(1,Ne), B_lam,...
%         'k--','LineWidth',lwt1);
% end

%% Axes and view
set(gca,'FontSize',fsz);
xlabel('\epsilon','FontSize',fsz);
ylabel('\gamma','FontSize',fsz);
zlabel('\lambda','FontSize',fsz);

axis([-mu emax gmin gmax lmin lmax])
view(3)
camlight headlight
lighting gouraud
material dull

save_count = save_count + 1;
fig_name = sprintf("fig%d",save_count);
saveas(gcf,fig_name,'epsc');

%% Curves on gamma planes
fcount = 1;
for g0 = gamma_planes
    fcount = fcount + 1;
    figure(fcount); clf;
    hold on; grid;
    title(strcat("\gamma = ",sprintf("%.0f",g0)),'FontSize',fsz);
    [~,jg] = min(abs(gam-g0));

    % Hysteresis curves
    plot(ep1(jg,:),lam,...
        'color',colH,'LineWidth',lwt1);

    plot(ep2(jg,:), lam,...
        'color',colH,'LineWidth',lwt1);

    % Bifurcation curve
    plot(B_ep, B_lam,...
        'color',colB,'LineWidth',lwt1);
    axis([-mu emax lmin lmax])
    set(gca,'Fontsize',fsz);
    xlabel('\epsilon','FontSize',fsz);
    ylabel('\lambda','FontSize',fsz);


end
%% Plots x vs sigma
detJ = @(x,mu,ep,sig,gam) (mu+ep - x).*(mu+ep - 3*x) + (sig - gam*x).*(sig - 3*gam*x);
trJ = @(x,mu,ep,sig,gam) 2*(mu+ep - 2*x);

fc = fcount;
fcount = 1;
l = 1;
eB = 3*(l^2*mu*0.25)^(1/3) - mu;
for g0 = gamma_planes
    [xh,sh,Ah] = Hfun(mu,g0,l);
    if g0 == 0
        ep_vals = [0.5, Ah(1) - mu,0.8, eB, 1];
    else
        ep_vals = sort([0, Ah(2) - mu,0.5, Ah(1) - mu,eB,1],'ascend');
    end
    fcount = fcount + 1;
    figure(fcount);
    plot([ep_vals(1)-0.3,ep_vals(end) + 0.3],l*[1,1],'Linewidth',2,'Linestyle','--','Color','k');
    plot(ep_vals,l*ones(size(ep_vals)),'.','Markersize',50,'Markeredgecolor','k');
    plot(ep_vals,l*ones(size(ep_vals)),'.','Markersize',30,'Markeredgecolor','y');
    save_count = save_count + 1;
    fig_name = sprintf("fig%d",save_count);
    saveas(gcf,fig_name,'epsc');
    gaux = 1+g0^2;
    F = @(x,s,A,mu,g,gaux,l)x.^3*gaux -2*(A+g*s).*x.^2 +(A^2+s.^2).*x - l^2*mu;
    for j = 1 : length(ep_vals)
        [s,x] = meshgrid(linspace(-2,2,101),linspace(0,2,101));
        f = F(x,s,mu+ep_vals(j),mu,g0,gaux,l);
        fc = fc+1;
        figure(fc); clf;
        c = contour(s,x,f,[0,0],'Linewidth',0.5);
        Ncontour = size(c,2);
        N1 = c(2,1);
        sx = c(:,2:2+N1-1);
        dJ = detJ(sx(2,:),mu,ep_vals(j),sx(1,:),g0);
        tJ = trJ(sx(2,:),mu,ep_vals(j),sx(1,:),g0);
        ind_stab = find(dJ > 0 & tJ < 0);
        ind_unst = find(dJ < 0 | tJ > 0);
        hold on; grid;
        plot_sigma_x(sx,ind_stab,'-',[0,0.7,0],lwt1);
        plot_sigma_x(sx,ind_unst,'--','b',lwt1);
        [n,s_d0,x_d0] = find_zero(dJ,sx(1,:),sx(2,:));
        if n > 0 
            plot(s_d0,x_d0,'.','Color','b','MarkerSize',30);
        end
        [n,s_t0,x_t0] = find_zero(tJ,sx(1,:),sx(2,:));
        if n > 0
            plot(s_t0,x_t0,'.','Color',[0.9,0,0],'MarkerSize',30);
        end
        if Ncontour > N1 + 1
            sx = c(:,N1+3:end);
            dJ = detJ(sx(2,:),mu,ep_vals(j),sx(1,:),g0);
            tJ = trJ(sx(2,:),mu,ep_vals(j),sx(1,:),g0);
            ind_stab = find(dJ > 0 & tJ < 0);
            ind_unst = find(dJ < 0 | tJ > 0);
            plot_sigma_x(sx,ind_stab,'-',[0,0.7,0],lwt1);
            plot_sigma_x(sx,ind_unst,'--','b',lwt1);
            [n,s_d0,x_d0] = find_zero(dJ,sx(1,:),sx(2,:));
            if n > 0 
                plot(s_d0,x_d0,'.','Color','b','MarkerSize',30);
            end
            [n,s_t0,x_t0] = find_zero(tJ,sx(1,:),sx(2,:));
            if n > 0
                plot(s_t0,x_t0,'.','Color',[0.9,0,0],'MarkerSize',30);
            end
        end
        if abs(ep_vals(j)- (Ah(1) - mu)) < 1e-10
            plot(sh(1),xh,'.','Markersize',50,'Color','k');
            plot(sh(1),xh,'.','MarkerSize',30,'Color',[115/255,253/255,1]);
        end
        if abs(ep_vals(j)- (Ah(2) - mu)) < 1e-10
            plot(sh(2),xh,'.','Markersize',50,'Markeredgecolor','k');
            plot(sh(2),xh,'.','MarkerSize',30,'Color',[115/255,253/255,1]);
        end
        title(strcat("\gamma = ",sprintf("%.0f",g0),...
            " mu = 0.2, \epsilon = ",sprintf("%.2f",ep_vals(j))),...
            'FontSize',fsz);
        % c = contour(s,x,f,[0,0],'Linewidth',lwt);
        axis([-2,2,0,2])
        set(gca,'Fontsize',fsz);
        xlabel('\epsilon','FontSize',fsz);
        ylabel('|u|^2','FontSize',fsz);
        save_count = save_count + 1;
        fig_name = sprintf("fig%d",save_count);
        saveas(gcf,fig_name,'epsc');

   end
end




end
%% ------------------------------------------------------------------------
function [x,s,A] = Hfun(mu,gam,lam)

if abs(gam) > 1e-10
    C = 1.5*(lam^2*mu*(1+gam^2)^2)^(1/3);
    a = 1 + gam^2;
    b = C;
    c = C^2 - gam^2*4*C^2/(3*(1+gam^2));
    D4 = b^2 - a*c;

    A1 = (b + sqrt(D4))/a;
    A2 = (b - sqrt(D4))/a;

    s1 = (C-A1)/gam;
    s2 = (C-A2)/gam;

    x = (2*C/a)/3;
    A = [A1 A2];
    s = [s1 s2];
else
    A1 = 1.5*lam^(2/3)*mu^(1/3);
    s1 =  A1/sqrt(3);
    s2 = -s1;

    x = 2*A1/3;
    A = [A1 A1];
    s = [s1 s2];
end

end

%%
function [ncomp,ifirst,ilast] = find_conn_intervals(ind)
d = ind(2:end) - ind(1:end-1);
ijump = find(d > 1);
if isempty(ijump)
    ifirst = 1;
    ilast = length(ind);
    ncomp = 1;
else
    njump = length(ijump);
    ncomp = njump + 1;
    ifirst = zeros(1,ncomp);
    ilast = zeros(1,ncomp);
    ifirst(1) = 1;
    for i = 1 : njump
        ilast(i) = ijump(i);
        ifirst(i+1) = ijump(i)+1;
    end
    ilast(ncomp) = length(ind);
end
end

function plot_sigma_x(sx,ind,lstyle,col,lwt)
[ncomp,ifirst,ilast] = find_conn_intervals(ind);
for j = 1 : ncomp
    plot(sx(1,ind(ifirst(j):ilast(j))),sx(2,ind(ifirst(j):ilast(j))), ...
        'LineStyle',lstyle,'Color',col,'LineWidth',lwt);
end
end

%%
function [n,s0,x0] = find_zero(f,s,x)
ind = find(f(1:end-1).*f(2:end) <= 0);
if ~isempty(ind)
    n = length(ind);
    s0 = zeros(1,n);
    x0 = zeros(1,n);
    for j = 1 : length(ind)
        a = -f(ind(j))/(f(ind(j) + 1)-f(ind(j)));
        s0(j) = s(ind(j)) + a*(s(ind(j)+1)-s(ind(j)));
        x0(j) = x(ind(j)) + a*(x(ind(j)+1)-x(ind(j)));
    end
else
    n = 0;
    s0 = 0;
    x0 = 0;
end
end