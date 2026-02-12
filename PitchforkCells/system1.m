function system1()
% close all
% if flag == 0, plot the phase diagram
% if flag == 1, plot bifurcation diagrams
% if flag == 2, plot basins
% if flag == 3, plot jump in y as a result of mu changing from zero to a
% positive value
% if flag == 4, plot the transition time for y to reach equilibrium after a
% sudden jump in mu from zero to a positive number

flag = 3;

lambda = 1;
fsz = 20;
lwt = 5;
efun = @(mu)(1.5*sqrt(3)*lambda*sqrt(mu)).^(2/3)-mu;

mu_star = 1.5*sqrt(3)*lambda;
mu = linspace(0,1,100000);
ee = linspace(0,1.5,100);
ep = efun(mu);
mumin = 0;
mumax = 3;
emin = -1.5; 
emax = 1.5;
%%
col(1,:) = [0, 0.4470, 0.7410]; % Light blue 
col(2,:) = [0.8500, 0.3250, 0.0980]; % Orange 
col(3,:) = [0.9290, 0.6940, 0.1250]; % Light yellow: 
col(4,:) = [0.4940, 0.1840, 0.5560]; % Purple: 
col(5,:) = [0.4660, 0.6740, 0.1880]; % Green: 
col(6,:) = [0.3010, 0.7450, 0.9330]; % Cyan: 
col(7,:) = [0.6350, 0.0780, 0.1840]; % Reddish-brown:
col(8,:) = [0,0,0];
col(9,:) = [0,0,1];
col(10,:) = [1,0,0];
col(11,:) = [0,0.5,0];
%%
if flag == 0
    figure(2);
    clf; hold on; 
    plot(ep,mu,'linewidth',lwt);
    plot(ee,4*ee.^3/(27*lambda^2),'--','LineWidth',lwt,'color',[0,0,0.5]);
    set(gca,'Xscale','log','Yscale','log','Fontsize',fsz)
    % xlabel('\epsilon','FontSize',fsz)
    % ylabel('\mu','Fontsize',fsz)
    axis([4e-2,0.8,1e-5,0.5])
    
    
    figure(1);
    mu_star = 1.5*sqrt(3)*lambda;
    mu1 = linspace(0,0.5,100);
    ep1 = efun(mu1);
    mu2 = linspace(0.5,mu_star,100);
    ep2 = efun(mu2);
    mu3 = linspace(mu_star,3,100);
    ep3 = efun(mu3);

    % test
    % roots = find_cubic_roots_mu(ee,lambda);
    
    mu = linspace(0,3,1000);
    ep = efun(mu);
    
    clf; hold on; 
    verts = [[ep;mu],[emax;mumax],[emax;mumin]];
    patch(verts(1,:),verts(2,:),[0,0,1]);
    alpha(0.2);
    
    verts = [[ep;mu],[emin;mumax],[emin;mumin],[0;0]];
    patch(verts(1,:),verts(2,:),[0,0.5,0]);
    alpha(0.2);
    
    plot(ep1,mu1,'linewidth',lwt);
    plot(ep2,mu2,'linewidth',lwt);
    plot(ep3,mu3,'linewidth',lwt);
    plot([-1.5,0],[1.5,0],'--','Linewidt',1,'color','k')
    plot([0,0],[mumin,mumax],'--','Linewidt',1,'color','k')
    
    % test
    % plot(ee,real(roots(:,1)),'.');
    % plot(ee,real(roots(:,2)),'.');
    % plot(ee,real(roots(:,3)),'.');
    
    set(gca,'Fontsize',fsz)
    % xlabel('\epsilon','FontSize',fsz)
    % ylabel('\mu','Fontsize',fsz)
    axis([-1.5,1.5,0,3])
    plot([4e-2,0.8,0.8,4e-2,4e-2],[1e-5,1e-5,0.5,0.5,1e-5],'k')
    
    %%
    figure(3);
    mu_max = 6;
    mu_min = -2;
    emin = -2;
    emax = 2;
    mu_star = 1.5*sqrt(3)*lambda;
    mu1 = linspace(0,0.5,100);
    ep1 = efun(mu1);
    mu2 = linspace(0.5,mu_star,100);
    ep2 = efun(mu2);
    mu3 = linspace(mu_star,mu_max,100);
    ep3 = efun(mu3);
    
    mu = linspace(0,mu_max,1000);
    ep = efun(mu);
    
    clf; hold on; 
    verts = [[mu;ep],[mu_max;emax],[mu_min;emax]];
    patch(verts(1,:),verts(2,:),[0,0,1]);
    alpha(0.2);
    
    verts = [[mu;ep],[0;emin],[0;emax],[0;0]];
    patch(verts(1,:),verts(2,:),[0,0.5,0]);
    alpha(0.2);
    
    verts = [[mu_min;emax],[0;emax],[0;0]];
    patch(verts(1,:),verts(2,:),[0.5,0.5,0],'EdgeColor','none');
    alpha(0.2);
    
    
    plot(mu1,ep1,'linewidth',lwt);
    plot(mu2,ep2,'linewidth',lwt);
    plot(mu3,ep3,'linewidth',lwt);
    % plot([1.5,0],[-1.5,0],'--','Linewidt',1,'color','k')
    plot([0,0],[emin,emax],'Linewidt',1,'color','k')
    plot([mu_min,mu_max],[-0.8,-0.8],'Linewidt',1,'color','k')
    plot([mu_min,mu_max],[0.8,0.8],'Linewidt',1,'color','k')
    plot([mu_min,mu_max],[0,0],'Linewidt',1,'color','k')
    plot([mu_min,mu_max],[1.2,1.2],'Linewidt',1,'color','k')
    set(gca,'Fontsize',fsz)
    axis([mu_min,mu_max,emin,emax])
    % plot([1e-5,1e-5,0.5,0.5,1e-5],[4e-2,0.8,0.8,4e-2,4e-2],'k')
    % ylabel('\epsilon','FontSize',28)
    % xlabel('\mu','Fontsize',28)
    grid on
    % daspect([1,1,1])
end
%% epsilon fixed
if flag == 1
    r0 = @(mu,ep)sqrt(ep) - 0.5*lambda*sqrt(mu)./ep;
    r1 = @(mu,ep)-sqrt(ep) - 0.5*lambda*sqrt(mu)./ep;
    r2 = @(mu,ep)lambda*sqrt(mu)./ep;
    
    ep = 0.1;
    mumax = 4*ep^3/(27*lambda^2);
    N = 100;
    mu = linspace(0,mumax,N)';
    roots = find_cubic_roots(mu,ep,lambda);
    
    figure(4)
    clf; hold on;
    plot(mu,roots(:,1),'LineWidth',lwt);
    plot(mu,roots(:,2),'LineWidth',lwt);
    plot(mu,roots(:,3),'LineWidth',lwt);
    plot(mu,r0(mu,ep),'--','LineWidth',lwt,'Color',[0,0,0.5]);
    plot(mu,r1(mu,ep),'--','LineWidth',lwt,'Color',[0.5,0,0]);
    plot(mu,r2(mu,ep),'--','LineWidth',lwt,'Color',[0.5,0.5,0]);
    % plot(mu,sqrt(ep)*ones(size(mu)),'LineWidth',2,'color',[0,0,0.5]);
    % plot(mu,-sqrt(ep)*ones(size(mu)),'LineWidth',2,'color',[0.5,0,0]);
    % plot(mu,zeros(size(mu)),'LineWidth',2,'color',[0.5,0.5,0]);
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Roots of p_{+}(y)','Fontsize',fsz)
    
    %% plot bifurcation diagram at epsilon = 0
    color1 = [0,0,1]; %[0, 0.4470, 0.7410]; % Light blue 
    color2 = [1,0,0]; %[0.8500, 0.3250, 0.0980]; % Orange 
    color3 = [0.9290, 0.6940, 0.1250]; % Light yellow: 
    color4 = [0.4940, 0.1840, 0.5560]; % Purple: 
    color5 = [0.4660, 0.6740, 0.1880]; % Green: 
    color6 = [0.3010, 0.7450, 0.9330]; % Cyan: 
    color7 = [0.6350, 0.0780, 0.1840]; % Reddish-brown:
    mu_star = 0.5*3*sqrt(3);
    mu_max = 6;
    mu_min = -2;
    Nmu = 1000;
    aa = 1e-10;
    mu = linspace(aa,mu_max,Nmu);
    lwt = 3;
    lwtu = 1;
    %
    fig_x = figure;
    hold on
    grid on
    plot([mu_min,0],[0,0],'Color','k','LineWidth',lwt);
    plot([0,mu_max],[0,0],'--','Color','k','LineWidth',lwtu);
    plot(mu,sqrt(mu),'Color',color1,'LineWidth',lwt);
    plot(mu,-sqrt(mu),'Color',color2,'LineWidth',lwt);
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Equilibrium x','Fontsize',fsz)
    %
    fig_y0 = figure;
    hold on
    grid on
    plot([mu_min,0],[0,0],'Color','k','LineWidth',lwt);
    plot([0,mu_max],[0,0],'--','Color','k','LineWidth',lwtu);
    plot(mu,sqrt(mu),'--','Color','k','LineWidth',lwtu);
    plot(mu,-sqrt(mu),'--','Color','k','LineWidth',lwtu);
    
    roots = find_cubic_roots(mu,0,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    % completion
    mu_aux = linspace(0,aa, 20);
    y = mu_aux.^(1/6);
    plot(mu_aux,-y,'Color',color1,'LineWidth',lwt);
    plot(mu_aux,y,'Color',color2,'LineWidth',lwt);
    
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Equilibrium y','Fontsize',fsz)
    
    %% plot bifurcation diagram at 0 < epsilon < lambda
    ep = 0.8;
    roots = find_cubic_roots_mu(ep,lambda);
    ind = find(roots > 0);
    mu1 = min(roots(ind));
    mu2 = max(roots(ind));
    
    fig_y1 = figure;
    hold on
    grid on
    plot([mu_min,-ep],[0,0],'Color','k','LineWidth',lwt);
    plot([-ep,mu_max],[0,0],'--','Color','k','LineWidth',lwtu);
    mu = linspace(-ep,0,Nmu);
    plot(mu,sqrt(mu + ep),'Color',color4,'LineWidth',lwt);
    plot(mu,-sqrt(mu + ep),'Color',color3,'LineWidth',lwt);
    mu = linspace(0,mu_max,Nmu);
    plot(mu,sqrt(mu + ep),'--','Color',color4,'LineWidth',lwtu);
    plot(mu,-sqrt(mu + ep),'--','Color',color3,'LineWidth',lwtu);
    
    mu = linspace(0,mu1,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    mu = linspace(mu1+aa,mu2-aa,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    mu = linspace(mu2+aa,mu_max,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Equilibrium y','Fontsize',fsz)
    
    %% plot bifurcation diagram at epsilon > lambda
    ep = 1.2;
    
    fig_y2 = figure;
    hold on
    grid on
    plot([mu_min,-ep],[0,0],'Color','k','LineWidth',lwt);
    plot([-ep,mu_max],[0,0],'--','Color','k','LineWidth',lwtu);
    mu = linspace(-ep,0,Nmu);
    plot(mu,sqrt(mu + ep),'Color',color4,'LineWidth',lwt);
    plot(mu,-sqrt(mu + ep),'Color',color3,'LineWidth',lwt);
    mu = linspace(0,mu_max,Nmu);
    plot(mu,sqrt(mu + ep),'--','Color',color4,'LineWidth',lwtu);
    plot(mu,-sqrt(mu + ep),'--','Color',color3,'LineWidth',lwtu);
    
    mu = linspace(0,mu_max,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Equilibrium y','Fontsize',fsz)
    
    %% plot bifurcation diagram at epsilon > lambda
    ep = -0.8;
    
    fig_y3 = figure;
    hold on
    grid on
    plot([mu_min,0],[0,0],'Color','k','LineWidth',lwt);
    plot([0,mu_max],[0,0],'--','Color','k','LineWidth',lwtu);
    mu = linspace(-ep,mu_max,Nmu);
    plot(mu,sqrt(mu + ep),'--','Color',color4,'LineWidth',lwtu);
    plot(mu,-sqrt(mu + ep),'--','Color',color3,'LineWidth',lwtu);
    
    mu = linspace(aa,mu_max,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    if ~isempty(ind)
        test = stability_test(sqrt(mu(ind)),real(r1(ind)),mu(ind),ep);
        Iunstab = find(test == 1);
        Istab = find(test == -1);
        if ~isempty(Istab)
            plot(mu(ind(Istab)),real(r1(ind(Istab))),'Color',color1,'LineWidth',lwt);
            plot(mu(ind(Istab)),-real(r1(ind(Istab))),'Color',color2,'LineWidth',lwt);
        end
        if ~isempty(Iunstab)
            plot(mu(ind(Iunstab)),real(r1(ind(Iunstab))),'--','Color',color1,'LineWidth',lwtu);
            plot(mu(ind(Iunstab)),-real(r1(ind(Iunstab))),'--','Color',color2,'LineWidth',lwtu);
        end
    end
    ind = find(abs(imag(r2)) < 1e-10);
    if ~isempty(ind)
        test = stability_test(sqrt(mu(ind)),real(r2(ind)),mu(ind),ep);
        Iunstab = find(test == 1);
        Istab = find(test == -1);
        if ~isempty(Istab)
            plot(mu(ind(Istab)),real(r2(ind(Istab))),'Color',color2,'LineWidth',lwt);
            plot(mu(ind(Istab)),-real(r2(ind(Istab))),'Color',color1,'LineWidth',lwt);
        end
        if ~isempty(Iunstab)
            plot(mu(ind(Iunstab)),real(r2(ind(Iunstab))),'--','Color',color1,'LineWidth',lwtu);
            plot(mu(ind(Iunstab)),-real(r2(ind(Iunstab))),'--','Color',color2,'LineWidth',lwtu);
        end
    end
    ind = find(abs(imag(r3)) < 1e-10);
    if ~isempty(ind)
        test = stability_test(sqrt(mu(ind)),real(r3(ind)),mu(ind),ep);
        Iunstab = find(test == 1);
        Istab = find(test == -1);
        if ~isempty(Istab)
            plot(mu(ind(Istab)),real(r3(ind(Istab))),'Color',color1,'LineWidth',lwt);
            plot(mu(ind(Istab)),-real(r3(ind(Istab))),'Color',color2,'LineWidth',lwt);
        end
        if ~isempty(Iunstab)
            plot(mu(ind(Iunstab)),real(r3(ind(Iunstab))),'--','Color',color1,'LineWidth',lwtu);
            plot(mu(ind(Iunstab)),-real(r3(ind(Iunstab))),'--','Color',color2,'LineWidth',lwtu);
        end
    end
    
    
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Equilibrium y','Fontsize',fsz)
end
%% Plot basins
if flag == 2
    col_attr = [0,0.5,0];
    mark_attr = '.'; sz_attr = 50;
    col_sad = [0.7,0,0];
    mark_sad = 'x'; sz_sad = 20;
    col_sour = [0.4,0,0.4];
    mark_sour = '*';sz_sour = 20;
    lws = 3;
    Tmax = 5e5;
    %
    ep = 0.1;
    roots = find_cubic_roots_mu(ep,lambda);
    fprintf("ep = %d\nlambda = %d\nroots = [%d,%d,%d]\n",ep,lambda,roots(1),roots(2),roots(3));
    %
    mu = [1e-4,2e-4];
    Nmu = length(mu);
    roots = find_cubic_roots(mu,ep,lambda);
    for jmu = 1 : Nmu
        r = roots(jmu,:);
        fprintf("mu = %d\n",mu(jmu));
        Ireal = find(abs(imag(r)) < 1e-10);
        Nreal = length(Ireal);
        r = real(r(Ireal));
        stab_test = sign(mu(jmu)+ep-3*r.^2);
        Istab = find(stab_test < 0);
        Iunstab = find(stab_test > 0);
        fig = figure;
        hold on
        grid on
        if ~isempty(Istab)
            plot(sqrt(mu(jmu))*ones(1,length(Istab)),r(Istab),...
                'Marker',mark_attr,'color',col_attr,'MarkerSize',sz_attr,...
                'Displayname','Attractors');
            plot(-sqrt(mu(jmu))*ones(1,length(Istab)),-r(Istab),...
                'Marker',mark_attr,'color',col_attr,'MarkerSize',sz_attr,...
                'Displayname','Attractors');
        end
        if ~isempty(Iunstab)
            plot(sqrt(mu(jmu))*ones(1,length(Iunstab)),r(Iunstab),...
                'Marker',mark_sad,'color',col_sad,'MarkerSize',sz_sad,...
                'LineWidth', lws,'Displayname','Saddles');
            plot(-sqrt(mu(jmu))*ones(1,length(Iunstab)),-r(Iunstab),...
                'Marker',mark_sad,'color',col_sad,'MarkerSize',sz_sad,...
                'LineWidth', lws,'Displayname','Saddles');
        end
        plot(0,0,'Marker',mark_sour,'color',col_sour,'MarkerSize',sz_sour,...
                'LineWidth', lws,'Displayname','Sources');
        ysad = sqrt(mu(jmu) + ep);
        plot([0,0],[ysad,-ysad],'Linestyle','none','Marker',mark_sad, 'LineWidth', lws,...
            'color',col_sad,'MarkerSize',sz_sad,...                
            'Displayname','Saddles');
        % saddles [0,0],[ysad,-ysad]
        % stable direction [0,1]
        % unstable direction [2\mu + 3\epsilon,\lambda]


        for j = 1 : 2
            sys1 = @(mu,ep,lam,y)[y(1)*mu - y(1).^3;y(2)*(mu + ep) - y(2).^3 - lam*y(1)];
            options = odeset('AbsTol',1e-10,'RelTol',1e-10,'events',@events);
            % backward in time along the stable direction of the saddles            
            fun = @(t,y)-sys1(mu(jmu),ep,lambda,y);
            v = [0;1];
            y0 = [0;sign(-3+2*j)*ysad]+1e-4*v;
            [~,Y] = ode45(fun,[0,Tmax],y0,options);
            plot(Y(:,1),Y(:,2),'Color','b','Linewidth',lws);
            y0 = [0;sign(-3+2*j)*ysad]-1e-4*v;
            [~,Y] = ode45(fun,[0,Tmax],y0,options);
            plot(Y(:,1),Y(:,2),'Color','b','Linewidth',lws);
            % forward in time along the unstable direction of the saddles            
            fun = @(t,y)sys1(mu(jmu),ep,lambda,y);
            v = [2*mu(jmu)+ep;lambda];
            y0 = [0;sign(-3+2*j)*ysad]+1e-4*v;
            [~,Y] = ode45(fun,[0,Tmax],y0,options);
            plot(Y(:,1),Y(:,2),'Color','r','Linewidth',lws);
            y0 = [0;sign(-3+2*j)*ysad]-1e-4*v;
            [~,Y] = ode45(fun,[0,Tmax],y0,options);
            plot(Y(:,1),Y(:,2),'Color','r','Linewidth',lws);
        end
        if ~isempty(Iunstab) % there is an unstable root
            for j = length(Iunstab)
                ysad = r(Iunstab);
                % integrate backward in time along the stable direction
                v_stable = [3*mu(jmu)+ep-3*ysad^2;lambda];
                fun = @(t,y)-sys1(mu(jmu),ep,lambda,y);
                for sig = 1:2
                    y0 = [sqrt(mu(jmu));ysad]+(-3+2*sig)*1e-4*v_stable;
                    [~,Y] = ode45(fun,[0,Tmax],y0,options);
                    % plot(Y(:,1),Y(:,2),'Color','b','Linewidth',lws);
                    if sig == 1
                        vpos1 = Y;
                    else
                        vpos2 = Y;
                    end
                    y0 = [-sqrt(mu(jmu));-ysad]+(-3+2*sig)*1e-4*v_stable;
                    [~,Y] = ode45(fun,[0,Tmax],y0,options);
                    % plot(Y(:,1),Y(:,2),'Color','b','Linewidth',lws);
                    if sig == 1
                        vneg1 = Y;
                    else
                        vneg2 = Y;
                    end
                end
                vrts = process_verts([vpos1(end:-1:1,:);vpos2]);
                verts_pos_major = [vrts;[1,1];[1,-1];[0,-1];[0,0]];
                verts_pos_minor = [vrts;[0,1];[0,0]];
                patch(verts_pos_minor(:,1),verts_pos_minor(:,2),'y');
                alpha(0.2);
                patch(-verts_pos_minor(:,1),-verts_pos_minor(:,2),[0,0.5,0]);
                alpha(0.2);
                patch(verts_pos_major(:,1),verts_pos_major(:,2),col(1,:));
                alpha(0.2);
                patch(-verts_pos_major(:,1),-verts_pos_major(:,2),col(2,:));
                alpha(0.2);
                plot(vpos1(:,1),vpos1(:,2),'Color','b','Linewidth',lws);
                plot(vpos2(:,1),vpos2(:,2),'Color','b','Linewidth',lws);
                plot(vneg1(:,1),vneg1(:,2),'Color','b','Linewidth',lws);
                plot(vneg2(:,1),vneg2(:,2),'Color','b','Linewidth',lws);
                % forward in time along the unstable direction of the saddles            
                fun = @(t,y)sys1(mu(jmu),ep,lambda,y);
                for sig = 1:2
                    v = [0;1]; % unstable direction
                    y0 = [sqrt(mu(jmu));ysad]+(-3+2*sig)*1e-4*v;
                    [~,Y] = ode45(fun,[0,Tmax],y0,options);
                    plot(Y(:,1),Y(:,2),'Color','r','Linewidth',lws);
                    y0 = [-sqrt(mu(jmu));-ysad]+(-3+2*sig)*1e-4*v;
                    [~,Y] = ode45(fun,[0,Tmax],y0,options);
                    plot(Y(:,1),Y(:,2),'Color','r','Linewidth',lws);
                end
            end
        else
            verts_pos_major = [[1,1];[1,-1];[0,-1];[0,1]];
            verts_neg_major = [[-1,-1];[-1,1];[0,1];[0,-1]];
            patch(verts_pos_major(:,1),verts_pos_major(:,2),col(1,:));
            alpha(0.2);
            patch(verts_neg_major(:,1),verts_neg_major(:,2),col(2,:));
            alpha(0.2);
        end

        if ~isempty(Istab)
            plot(sqrt(mu(jmu))*ones(1,length(Istab)),r(Istab),...
                'Marker',mark_attr,'color',col_attr,'MarkerSize',sz_attr,...
                'Displayname','Attractors');
            plot(-sqrt(mu(jmu))*ones(1,length(Istab)),-r(Istab),...
                'Marker',mark_attr,'color',col_attr,'MarkerSize',sz_attr,...
                'Displayname','Attractors');
        end
        if ~isempty(Iunstab)
            plot(sqrt(mu(jmu))*ones(1,length(Iunstab)),r(Iunstab),...
                'Marker',mark_sad,'color',col_sad,'MarkerSize',sz_sad,...
                'LineWidth', lws,'Displayname','Saddles');
            plot(-sqrt(mu(jmu))*ones(1,length(Iunstab)),-r(Iunstab),...
                'Marker',mark_sad,'color',col_sad,'MarkerSize',sz_sad,...
                'LineWidth', lws,'Displayname','Saddles');
        end
        plot(0,0,'Marker',mark_sour,'color',col_sour,'MarkerSize',sz_sour,...
                'LineWidth', lws,'Displayname','Sources');
        ysad = sqrt(mu(jmu) + ep);
        plot([0,0],[ysad,-ysad],'Linestyle','none','Marker',mark_sad, 'LineWidth', lws,...
            'color',col_sad,'MarkerSize',sz_sad,...                
            'Displayname','Saddles');
        xlabel('x','FontSize',fsz);
        ylabel('y','Fontsize',fsz);  
        set(gca,'Fontsize',fsz);
        axis([-0.015,0.015,-1,1])
    end

end

%% Plot possible jump in y as a result of mu jumping from 0 to a positive value

if flag == 3
    lambda = 1;
    Nmu = 100;
    aux = linspace(log(1e-8),0,Nmu);
    muvals = exp(aux);
    Ne = 10;
    epvals = linspace(-0.3,0.6,Ne);
    epvals = [-1e-1,-1e-2,0,1e-2,1e-1];
    Ne = length(epvals);
    yjump1 = zeros(Ne,Nmu);
    yjump2 = zeros(Ne,Nmu);
    yy = zeros(Ne,Nmu);
    xx = zeros(Ne,Nmu);
    for je = 1 : Ne
        ep = epvals(je);
        if ep <= 0
            y0 = 0;
        else
            y0 = sqrt(ep);
        end
        for jmu = 1 : Nmu
            mu = muvals(jmu);
            roots = find_cubic_roots(mu,ep,lambda);
            Ireal = abs(imag(roots))<1e-10;
            r = real(roots(Ireal));
            Istab = find(mu+ep-3*r.^2 < 0);
            rstab = r(Istab);
            if length(Istab) == 2
                rpos =rstab(rstab >= 0);
                rneg = rstab(rstab < 0);
                yjump1(je,jmu) = abs(rpos-y0);
                yjump2(je,jmu) = abs(-rneg-y0);
            else % only one sink
                yjump1(je,jmu) = abs(rstab-y0);
                yjump2(je,jmu) = abs(rstab+y0);
            end
        end
    end

    fig = figure;
    hold on
    lwt = 3;
    g = linspace(0,1,Ne)';
    z = zeros(size(g));
    col = [g,z,g(Ne:-1:1)];
    col = 0.8*parula(Ne);
    col = col(end:-1:1,:);
    col(3,:) = [1,0,0];
    
    for je = 1 : Ne
        roots = find_cubic_roots_mu(epvals(je),lambda);
        ind = find(roots > 0);
        mu1 = min(roots(ind));
        mu2 = max(roots(ind));
        fprintf("epsilon = %d, mu1 = %d, mu2 = %d\n",epvals(je),mu1,mu2);
        if mu1 > 1e-10
            plot([mu1,mu1],[min(min(yjump1)),max(max(yjump2))],'color',col(je,:),'Linewidth',3,'Linestyle',':');
        end
        sname = sprintf("%.0e",epvals(je));
        sname = strcat("\epsilon = ",sname);
        h(2*je-1) = plot(muvals,yjump1(je,:),'LineWidth',lwt,'color',col(je,:),'LineStyle','-','DisplayName',sname);
        h(je*2) = plot(muvals,yjump2(je,:),'LineWidth',lwt,'color',col(je,:),'LineStyle','-','DisplayName',sname);
    end
    set(gca,'Xscale','log','Yscale','log','Fontsize',fsz);
    % set(gca,'Fontsize',fsz);
    xlabel('\mu','FontSize',28)
    ylabel('Absolute jump in y','Fontsize',fsz)
    legend(h(1:2:end))
    axis([1e-8,1,min(min(yjump1)),max(max(yjump2))])
    grid on
        
end

%%
%% plot time required for y to reach a new equilibrium as mu changes from zero to a positive value

if flag == 4
    tic
    lambda = 1;
    sys1 = @(mu,ep,lam,y)[y(1)*mu - y(1).^3;y(2)*(mu + ep) - y(2).^3 - lam*y(1)];
    Nmu = 100;
    x0 = 1e-3;
    Tmax = 1e7;
    aux = linspace(log(1e-5),0,Nmu);
    muvals = exp(aux);
    Ne = 10;
    epvals = linspace(-0.3,0.6,Ne);
    epvals = [0,1e-1];
    Ne = length(epvals);
    T1 = zeros(Ne,Nmu);
    T2 = zeros(Ne,Nmu);
    for je = 1 : Ne
        ep = epvals(je);
        if ep <= 0
            y0 = 0;
        else
            y0 = sqrt(ep);
        end
        for jmu = 1 : Nmu
            mu = muvals(jmu);
            fun = @(t,y)sys1(muvals(jmu),epvals(je),lambda,y);
            roots = find_cubic_roots(mu,ep,lambda);
            Ireal = abs(imag(roots))<1e-10;
            r = real(roots(Ireal));
            Istab = mu+ep-3*r.^2 < 0;
            rstab = r(Istab);
            % if length(Istab) == 2
                % ystar = rstab(rstab >= 0);
                flag5events = @(t,y)Flag5Events(rstab,y);
                options = odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',flag5events);
                xy0 = [1e-3;y0];
                [T,~] = ode45(fun,[0,Tmax],xy0,options);
                T1(je,jmu) = T(end);
                %
                % ystar = rstab(rstab < 0);
                % flag5events = @(t,y)Flag5Events(ystar,y);
                % options = odeset('AbsTol',1e-9,'RelTol',1e-9,'Events',flag5events);
                % xy0 = [-1e-4;y0];
                % [T,~] = ode45(fun,[0,Tmax],xy0,options);
                % T2(je,jmu) = T(end);
            % else % only one sink
            %     ystar = rstab(rstab >= 0);
            %     flag5events = @(t,y)Flag5Events(rstab,y);
            %     options = odeset('AbsTol',1e-9,'RelTol',1e-9,'Events',flag5events);
            %     xy0 = [1e-4;-y0];
            %     [T,~] = ode45(fun,[0,Tmax],xy0,options);
            %     T1(je,jmu) = T(end);
            %     %
            %     ystar = rstab(rstab < 0);
            %     flag5events = @(t,y)Flag5Events(ystar,y);
            %     options = odeset('AbsTol',1e-9,'RelTol',1e-9,'Events',flag5events);
            %     xy0 = [1e-4;y0];
            %     [T,~] = ode45(fun,[0,Tmax],xy0,options);
            %     T2(je,jmu) = T(end);
            % end
        end
    end
    fig = figure;
    hold on
    lwt = 3;
    g = linspace(0,1,Ne)';
    z = zeros(size(g));
    col = [g,z,g(Ne:-1:1)];
    col = 0.8*parula(Ne);
    col = col(end:-1:1,:);
    col(1,:) = [1,0,0];
    col(2,:) = [0,0,0.5];
    for je = 1 : Ne
        roots = find_cubic_roots_mu(ep,lambda);
        ind = find(roots > 0);
        mu1 = min(roots(ind));
        mu2 = max(roots(ind));
        if mu1 > 1e-5
            plot([mu1,mu1],[1e-5,2e6],'color',col(je,:),'Linewidth',3,'Linestyle',':');
        end
        sname = sprintf("%.0e",epvals(je));
        sname = strcat("\epsilon = ",sname);
        h(je) = plot(muvals,T1(je,:),'LineWidth',lwt,'color',col(je,:),'LineStyle','-','DisplayName',sname);
        % h(je*2) = plot(muvals,T2(je,:),'LineWidth',lwt,'color',col(je,:),'LineStyle','--','DisplayName',sname);
    end
    set(gca,'Xscale','log','Yscale','log','Fontsize',fsz);
    % set(gca,'Fontsize',fsz);
    xlabel('\mu','FontSize',28)
    ylabel('Switching time','Fontsize',fsz)
    legend(h)
    axis([1e-5,1,10,1e6]);
    grid on
        
    toc
    xdata = log(muvals);
    ydata = log(T1(2,:));
    p = polyfit(xdata,ydata,1)

end
%% flag = 6. Plot bifurcation diagram at epsilon = 0.1
if flag == 6
    %% plot bifurcation diagram at 0 < epsilon < lambda
    color1 = [0,0,1]; %[0, 0.4470, 0.7410]; % Light blue 
    color2 = [1,0,0]; %[0.8500, 0.3250, 0.0980]; % Orange 
    color3 = [0.9290, 0.6940, 0.1250]; % Light yellow: 
    color4 = [0.4940, 0.1840, 0.5560]; % Purple: 
    color5 = [0.4660, 0.6740, 0.1880]; % Green: 
    color6 = [0.3010, 0.7450, 0.9330]; % Cyan: 
    color7 = [0.6350, 0.0780, 0.1840]; % Reddish-brown:

    ep = 0.1;
    roots = find_cubic_roots_mu(ep,lambda);
    ind = find(roots > 0);
    mu1 = min(roots(ind));
    mu2 = max(roots(ind));

    mu_min = -1e-4;
    mu_max = 2.5e-4;
    Nmu = 1000;
    aa = 1e-10;
    lwt = 3;
    lwtu = 1;

    fig_y1 = figure;
    hold on
    grid on
    % plot([mu_min,-ep],[0,0],'Color','k','LineWidth',lwt);
    plot([-ep,mu_max],[0,0],'--','Color','k','LineWidth',lwtu);
    mu = linspace(-ep,0,Nmu);
    plot(mu,sqrt(mu + ep),'Color',color4,'LineWidth',lwt);
    plot(mu,-sqrt(mu + ep),'Color',color3,'LineWidth',lwt);
    mu = linspace(0,mu_max,Nmu);
    plot(mu,sqrt(mu + ep),'--','Color',color4,'LineWidth',lwtu);
    plot(mu,-sqrt(mu + ep),'--','Color',color3,'LineWidth',lwtu);
    
    mu = linspace(0,mu1,Nmu);
    plot([mu1,mu1],[-2,2],'--','LineWidth',lwtu,'Color',[0.5,0,0.5]);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    mu = linspace(mu1+aa,mu2-aa,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    mu = linspace(mu2+aa,mu_max,Nmu);
    roots = find_cubic_roots(mu,ep,lambda);
    r1 = roots(:,1)';
    r3 = roots(:,2)'; % the root that always exists
    r2 = roots(:,3)';
    ind = find(abs(imag(r1)) < 1e-10);
    plot(mu(ind),real(r1(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r1(ind)),'Color',color2,'LineWidth',lwt);
    ind = find(abs(imag(r2)) < 1e-10);
    plot(mu(ind),real(r2(ind)),'--','Color',color1,'LineWidth',lwtu);
    plot(mu(ind),-real(r2(ind)),'--','Color',color2,'LineWidth',lwtu);
    ind = find(abs(imag(r3)) < 1e-10);
    plot(mu(ind),real(r3(ind)),'Color',color1,'LineWidth',lwt);
    plot(mu(ind),-real(r3(ind)),'Color',color2,'LineWidth',lwt);
    
    set(gca,'Fontsize',fsz)
    xlabel('\mu','FontSize',28)
    ylabel('Equilibrium y','Fontsize',fsz)
    


end

end % function system1
%%
function roots = find_cubic_roots(mu,ep,lambda)
roots = zeros(length(mu),3);
mu = mu(:);
a = mu + ep;
b = lambda*sqrt(mu);
aux = (1/3)*acos((-3*b./(2*a)).*sqrt(3./a));
roots(:,1) = 2*sqrt(a/3).*cos(aux);
roots(:,2) = 2*sqrt(a/3).*cos(aux + 2*pi/3);
roots(:,3) = 2*sqrt(a/3).*cos(aux + 4*pi/3);
end
%%
function roots = find_cubic_roots_mu(ep,lambda)
roots = zeros(length(ep),3);
ep = ep(:);
a = (1.5*sqrt(3)*lambda)^(2/3)*ones(size(ep));
b = ep;
aux = (1/3)*acos((-3*b./(2*a)).*sqrt(3./a));
roots(:,1) = 2*sqrt(a/3).*cos(aux);
roots(:,2) = 2*sqrt(a/3).*cos(aux + 2*pi/3);
roots(:,3) = 2*sqrt(a/3).*cos(aux + 4*pi/3);
roots = roots.^3;
end
%%
function test = stability_test(x,y,mu,ep)
test = max(sign(mu - 3*x.^2),sign(mu + ep -3*y.^2));
end
%%
function [value,isterminal,direction] = events(t,y)
        value = max(abs(y(1)) - 1,abs(y(2)) - 1);      % Detect when the height (y(1)) crosses zero
        isterminal = 1;    % Stop the integration at the event
        direction = 0;    % Detect only when height is decreasing (going from positive to negative)
end
%%
function vrts = process_verts(v)
N0 = length(v);
if N0 < 1000
    vrts = unique(v,'rows');
else
    v = unique(v,'rows');
    dv = v(2:end,:) - v(1:end-1,:);
    dd = [0;cumsum(sqrt(sum(dv.^2,2)))];
    g = linspace(0,dd(end),1000)';
    vrts = interp1(dd,v,g);
end
end
%%
function [value,isterminal,direction] = Flag5Events(ystar,y)
        value = min(abs(y(2)-ystar))-1e-3;      
        isterminal = 1;    % Stop the integration at the event
        direction = 0;    % Detect only when height is decreasing (going from positive to negative)
end
