function system4_amplitude_comparison()
lam = 1;
gam = 1;

fsz = 20;
Nmu = 1000;
muvals = linspace(1e-10,0.2,Nmu);

func = @(x,mu,ep,sig,gam) (mu+ep - x)^2*x + (sig - gam*x)^2*x - lam^2*mu;
detJ = @(x,mu,ep,sig,gam) (mu+ep - x).*(mu+ep - 3*x) + (sig - gam*x).*(sig - 3*gam*x);
trJ = @(x,mu,ep,sig,gam) 2*(mu+ep - 2*x);
%% sigma = 0
lwt = 3;

sig = 0;
figure(1); clf; hold on; grid
epvals = [-0.2,0,0.2];
Ne = length(epvals);
col = 0.7*winter(Ne);
h(1) = plot(muvals,sqrt(muvals),'color','k','Linewidth',lwt,'DisplayName','|z_1|');
for j = 1 : Ne
    ep = epvals(j);
    x0 = ep;
    func = @(x,mu,ep,sig,gam) (mu+ep - x)^2*x + (sig - gam*x)^2*x - lam^2*mu;
    y = zeros(size(muvals));
    for k = 1 : Nmu
        mu = muvals(k);
        f = @(x)func(x,mu,ep,sig,gam);
        y(k) = fzero(f,x0);
        x0 = y(k);
    end
    lbl = strcat("|z_2|, eps = ",sprintf("%.1f",ep));
    if ep == 0
        col(j,:) = [1,0,0];
    end
    ind = find(detJ(y,muvals,ep,sig,gam) > 0 & trJ(y,muvals,ep,sig,gam) < 0);
    h(j+1) = plot(muvals(ind),sqrt(y(ind)),'Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
    ind = find(detJ(y,muvals,ep,sig,gam) <= 0 | trJ(y,muvals,ep,sig,gam) > 0);
    plot(muvals(ind),sqrt(y(ind)),'--','Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
end
% legend(h(1:Ne+1),'NumColumns', 2);
set(gca,'Fontsize',fsz);
xlabel("\mu",'FontSize',fsz);
ylabel("Amplitude","fontsize",fsz);
axis([0,max(muvals),0,1])
%%
figure(2); clf; hold on; grid
h(1) = plot(muvals,sqrt(muvals),'color','k','Linewidth',lwt,'DisplayName','|z_1|');

sig = 0.2;
for j = 1 : Ne
    ep = epvals(j);
    x0 = ep;
    func = @(x,mu,ep,sig,gam) (mu+ep - x)^2*x + (sig - gam*x)^2*x - lam^2*mu;
    y = zeros(size(muvals));
    for k = 1 : Nmu
        mu = muvals(k);
        f = @(x)func(x,mu,ep,sig,gam);
        y(k) = fzero(f,x0);
        x0 = y(k);
    end
    lbl = strcat("|z_2|, eps = ",sprintf("%.1f",ep));
    if ep == 0
        col(j,:) = [1,0,0];
    end
    ind = find(detJ(y,muvals,ep,sig,gam) > 0 & trJ(y,muvals,ep,sig,gam) < 0);
    plot(muvals(ind),sqrt(y(ind)),'Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
    ind = find(detJ(y,muvals,ep,sig,gam) <= 0 | trJ(y,muvals,ep,sig,gam) > 0);
    plot(muvals(ind),sqrt(y(ind)),'--','Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
end
% legend('NumColumns', 2);
set(gca,'Fontsize',fsz);
xlabel("\mu",'FontSize',fsz);
ylabel("Amplitude","fontsize",fsz);
axis([0,max(muvals),0,1])

%%
figure(3); clf; hold on; grid
h(1) = plot(muvals,sqrt(muvals),'color','k','Linewidth',lwt,'DisplayName','|z_1|');


sig = 0.5;
for j = 1 : Ne
    ep = epvals(j);
    x0 = ep;
    func = @(x,mu,ep,sig,gam) (mu+ep - x)^2*x + (sig - gam*x)^2*x - lam^2*mu;
    y = zeros(size(muvals));
    for k = 1 : Nmu
        mu = muvals(k);
        f = @(x)func(x,mu,ep,sig,gam);
        y(k) = fzero(f,x0);
        x0 = y(k);
    end
    lbl = strcat("|z_2|, eps = ",sprintf("%.1f",ep));
    if ep == 0
        col(j,:) = [1,0,0];
    end
    ind = find(detJ(y,muvals,ep,sig,gam) > 0 & trJ(y,muvals,ep,sig,gam) < 0);
    plot(muvals(ind),sqrt(y(ind)),'Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
    ind = find(detJ(y,muvals,ep,sig,gam) <= 0 | trJ(y,muvals,ep,sig,gam) > 0);
    plot(muvals(ind),sqrt(y(ind)),'--','Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
end



% legend('NumColumns', 2);
set(gca,'Fontsize',fsz);
xlabel("\mu",'FontSize',fsz);
ylabel("Amplitude","fontsize",fsz);
axis([0,max(muvals),0,1])


% %% sigma = 0.2
% sig = -0.2;
% %figure(1); clf; hold on; grid
% epvals = [-0.2,-0.1,0,0.1,0.2];
% Ne = length(epvals);
% col = 0.7*winter(Ne);
% %plot(muvals,sqrt(muvals),'color','k','Linewidth',lwt,'DisplayName','|z_1|');
% for j = 1 : Ne
%     ep = epvals(j);
%     x0 = ep;
%     func = @(x,mu,ep,sig,gam) (mu+ep - x)^2*x + (sig - gam*x)^2*x - lam^2*mu;
%     y = zeros(size(muvals));
%     for k = 1 : Nmu
%         mu = muvals(k);
%         f = @(x)func(x,mu,ep,sig,gam);
%         y(k) = fzero(f,x0);
%         x0 = y(k);
%     end
%     lbl = strcat("|z_2|, eps = ",sprintf("%.1f",ep));
%     if ep == 0
%         col(j,:) = [1,0,0];
%     end
%     p(j+Ne) = plot(muvals,sqrt(y),'--','Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
% end
% legend(h(1:Ne+1),'NumColumns', 2);
% set(gca,'Fontsize',fsz);
% xlabel("\mu",'FontSize',fsz);
% ylabel("Amplitude","fontsize",fsz);
% %% sigma = -0.5
% sig = -0.5;
% %figure(1); clf; hold on; grid
% epvals = [-0.2,-0.1,0,0.1,0.2];
% Ne = length(epvals);
% col = 0.7*winter(Ne);
% %plot(muvals,sqrt(muvals),'color','k','Linewidth',lwt,'DisplayName','|z_1|');
% for j = 1 : Ne
%     ep = epvals(j);
%     x0 = ep;
%     func = @(x,mu,ep,sig,gam) (mu+ep - x)^2*x + (sig - gam*x)^2*x - lam^2*mu;
%     y = zeros(size(muvals));
%     for k = 1 : Nmu
%         mu = muvals(k);
%         f = @(x)func(x,mu,ep,sig,gam);
%         y(k) = fzero(f,x0);
%         x0 = y(k);
%     end
%     lbl = strcat("|z_2|, eps = ",sprintf("%.1f",ep));
%     if ep == 0
%         col(j,:) = [1,0,0];
%     end
%     p(j+Ne) = plot(muvals,sqrt(y),':','Linewidth',lwt,'color',col(j,:),'DisplayName',lbl);
% end
% legend(h(1:Ne+1),'NumColumns', 2);
% set(gca,'Fontsize',fsz);
% xlabel("\mu",'FontSize',fsz);
% ylabel("Amplitude","fontsize",fsz);
% 
end



