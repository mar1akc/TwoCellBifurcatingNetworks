function system2_reduced_simulations()
% mu = 3, sigma = 0
% mu = 2, sigma = 0.5
% mu = 1.91, sigma = 1.05
fsz = 20;
mu = 0.5; 
sigma = 1.3;
par = [mu,sigma];
options = odeset('AbsTol',1e-9,'RelTol',1e-9);
v0 = [0.4;1e-10];
fun = @(t,v)func(v,par);
[T,Y] = ode45(fun,[0,20],v0,options);

figure(1); clf; hold on; grid;
plot(T,Y(:,1),'Linewidth',2,'DisplayName','v_R');
plot(T,Y(:,2),'Linewidth',2,'DisplayName','v_I');
set(gca,'FontSize',fsz);
legend();

figure(2); clf; hold on; grid;
a = linspace(-1.5,1.5,50);
[vr,vi] = meshgrid(a,a);
[dvr,dvi] = func1(vr,vi,par);
aux = sqrt(dvr.^2 + dvi.^2+1e-10);
quiver(vr,vi,dvr./aux,dvi./aux);
set(gca,'FontSize',fsz);

tskip = 300;

v0 = [1;1e-10];
[T,Y] = ode45(fun,[0,tskip],v0,options);

[T,Y] = ode45(fun,[0,60],Y(end,:),options);
plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);

% v0 = [1;-1.4];
% [T,Y] = ode45(fun,[0,30],v0,options);
% plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);



p = [mu^2,-2*mu^2,mu^2 + sigma^2,-1];
root = roots(p);
ind = find(abs(imag(root)) < 1e-10 & real(root) > 0);
if ~isempty(ind)
    N = length(ind);
    vstar = zeros(N,2);
    for j = 1 : N
        v2 = real(root(ind(j)));
        detJ = mu^2*(1-v2)*(1-3*v2) + sigma^2;
        trJ = 2*mu*(1-2*v2);
        aux = (mu^2*(1-v2)^2 + sigma^2);
        vstar(j,1) = mu*(1-v2)/aux;
        vstar(j,2) = -sigma/aux;

        x = vstar(j,1);
        y = vstar(j,2);
        
        J = [mu*(1-v2) - 2*mu*x^2,-sigma-2*mu*x*y;
            sigma - 2*mu*x*y,mu*(1-v2) - 2*mu*y^2];


        if detJ > 0 && trJ < 0
            plot(vstar(j,1),vstar(j,2),'.','MarkerSize',50);
        else 
            if detJ > 0 && trJ >= 0
                tmax = 200;
                    [V,E] = eig(J);
                    [evals,jsort] = sort(diag(E),'descend');
                    V = V(:,jsort);
                    % 
                    % v0 = [x;y] + 1e-3*V(:,1);
                    % [T,Y] = ode45(fun,[0,tmax],v0,options);
                    % plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);
                    % 
                    % v0 = [x;y] - 1e-3*V(:,1);
                    % [T,Y] = ode45(fun,[0,tmax],v0,options);
                    % plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);
                    % 
                    %  v0 = [x;y] + 1e-3*V(:,2);
                    % [T,Y] = ode45(fun,[0,tmax],v0,options);
                    % plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);
                    % 
                    % v0 = [x;y] - 1e-3*V(:,2);
                    % [T,Y] = ode45(fun,[0,tmax],v0,options);
                    % plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);

                plot(vstar(j,1),vstar(j,2),'*','MarkerSize',20,'LineWidth', 3,...
                    'color',[0.4940, 0.1840, 0.5560]);
            else
                if detJ <= 0
                    [V,E] = eig(J);
                    [evals,jsort] = sort(diag(E),'descend');
                    V = V(:,jsort);

                    v0 = [x;y] + 1e-2*V(:,1);
                    [T,Y] = ode45(fun,[0,100],v0,options);
                    plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);

                    v0 = [x;y] - 1e-2*V(:,1);
                    [T,Y] = ode45(fun,[0,100],v0,options);
                    plot(Y(:,1),Y(:,2),'Linewidth',3,'color',[0,0.5,0]);


                    plot(vstar(j,1),vstar(j,2),'x','Color','k','MarkerSize',20,'LineWidth', 3);
                end
            end
        end
    end
    for j = 1 : N
        v2 = real(root(ind(j)));
        detJ = mu^2*(1-v2)*(1-3*v2) + sigma^2;
        trJ = 2*mu*(1-2*v2);
        aux = (mu^2*(1-v2)^2 + sigma^2);
        vstar(j,1) = mu*(1-v2)/aux;
        vstar(j,2) = -sigma/aux;
        if detJ > 0 && trJ < 0
            plot(vstar(j,1),vstar(j,2),'.','MarkerSize',50,'color',[0.7,0,0]);
        else 
            if detJ > 0 && trJ >= 0
                plot(vstar(j,1),vstar(j,2),'*','MarkerSize',20,'LineWidth', 3,...
                    'color',[0.4940, 0.1840, 0.5560]);
            else
                if detJ <= 0
                    plot(vstar(j,1),vstar(j,2),'x','Color','k','MarkerSize',20,'LineWidth', 3);
                end
            end
        end
    end
else
    fprint("|v|^2 is not found\n");
end
daspect([1,1,1])
axis tight
end
%%
function dv = func(v,par)
dv = zeros(2,1);
aux = v(1).^2 + v(2).^2;
dv(1) = par(1)*(1-aux).*v(1) - par(2)*v(2) - 1;
dv(2) = par(1)*(1-aux).*v(2) + par(2)*v(1);
end

%%
function [dvr,dvi] = func1(vr,vi,par)
dv = zeros(2,1);
aux = vr.^2 + vi.^2;
dvr = par(1)*(1-aux).*vr - par(2)*vi - 1;
dvi = par(1)*(1-aux).*vi + par(2)*vr;
end