function System2_invariant_tori()
%close all

% if flag == 0: plot graphs
% if flag == 1: oscillator 2 in the co-rotating frame
% if flag == 2: compute Lyapunov exponents

flag = 0;
% Make sample plots
lwt = 3;
fsz = 20;
col(1,:) = [0, 0.4470, 0.7410]; % Light blue 
col(2,:) = [0.8500, 0.3250, 0.0980]; % Orange 

Tplot = 50;
Tskip = 1000;
omega = 1;
func = @(sigma,mu,y)[(mu + 1i*omega)*y(1) - y(1)*abs(y(1)).^2;...
    (mu + 1i*(omega+sigma))*y(2) - y(2)*abs(y(2)).^2 - y(1)];
options = odeset("AbsTol",1e-9,"RelTol",1e-9);
options1 = odeset("AbsTol",1e-9,"RelTol",1e-9,"events",@events);
par = [%1.2,2;
    % -1.2,2;
    % 1.25,2.75;
    1.4, 0.5;
    % -1.4,1.75;
    %1.4,1.75
    ];
npar = size(par,1);

%% plot graphs
if flag == 0

for j = 1:npar
    odefun = @(t,y)func(par(j,1),par(j,2),y);
    [~,Y1] = ode45(odefun,[0,Tskip],[1,0],options);
    y0 = Y1(end,:);
    [T2,Y2] = ode45(odefun,[Tskip,Tskip+Tplot],y0,options1);

    %
    fig = figure;
    hold on;
    h(1) = plot(T2,abs(Y2(:,1)),"LineWidth",3,"DisplayName","abs(z_1)","color",col(1,:));
    h(2) = plot(T2,abs(Y2(:,2)),"LineWidth",3,"DisplayName","abs(z_2)","color",col(2,:));
    h(3) = plot(T2,-abs(Y2(:,1)),"LineWidth",3,"DisplayName","abs(z_1)","color",col(1,:));
    h(4) = plot(T2,-abs(Y2(:,2)),"LineWidth",3,"DisplayName","abs(z_2)","color",col(2,:));
    % set(gca,'Fontsize',fsz);
    % xlabel("Time, t","FontSize",fsz);
    % ylabel("abs(z_1), abs(z_2)","FontSize",fsz);
    % tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    % title(tit,"FontSize",fsz);
    % legend()
    % figname = sprintf("sigma%.2fmu%.2f.eps",par(j,1),par(j,2)); 
    % 
    % fig = figure;
    % hold on;
    h(5) = plot(T2,real(Y2(:,1)),"LineWidth",1,"DisplayName","Re(z_1)","color",col(1,:));
    h(6) = plot(T2,real(Y2(:,2)),"LineWidth",1,"DisplayName","Re(z_2)","color",col(2,:));
    set(gca,'Fontsize',fsz);
    xlabel("Time, t","FontSize",fsz);
    ylabel("z_1, z_2","FontSize",fsz);
    tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    title(tit,"FontSize",fsz);
    legend(h([1,2,5,6]))
    figname = sprintf("Torus_sigma%.2fmu%.2f.eps",par(j,1),par(j,2)); 
    % saveas(fig,figname,"epsc");
    %
    [T2,Y2,te,ye,ie] = ode45(odefun,[Tskip,10*Tskip],y0,options1);
    fig = figure;
    plot(real(ye(:,2)),imag(ye(:,2)),'.')
    set(gca,'Fontsize',fsz);
    xlabel("Re(z_2)","FontSize",fsz);
    ylabel("Im(z_2)","FontSize",fsz);
    tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    title(tit,"FontSize",fsz);
    daspect([1,1,1])
    %
    fig = figure;
    % if j > 8
        [T2,Y2] = ode45(odefun,[0,Tskip],y0,options);
    % end
    plot3(real(Y2(:,2)),imag(Y2(:,2)),real(Y2(:,1)),"LineWidth",0.25);
    set(gca,'Fontsize',fsz);
    xlabel("Re(z_2)","FontSize",fsz);
    ylabel("Im(z_2)","FontSize",fsz);
    zlabel("Re(z_1)","FontSize",fsz);
    title(tit,"FontSize",fsz);
    figname = sprintf("sigma%.2fmu%.2f_phase.eps",par(j,1),par(j,2));
    daspect([1,1,1])
    grid
    % saveas(fig,figname,"epsc");

end
end
%% co-rotating frame
if flag == 1

func1 = @(sigma,mu,y)(mu + 1i*sigma)*y - y*abs(y)^2 - sqrt(mu);

for j = 1:npar
    odefun = @(t,y)func1(par(j,1),par(j,2),y);
    [~,Y1] = ode45(odefun,[0,Tskip],1+0i,options);
    y0 = Y1(end);
    [T2,Y2,te2,ye2,ie2] = ode45(odefun,[Tskip,Tskip+Tplot],y0,options1);
    %
    fig = figure;
    hold on;
    h1 = plot(T2,abs(Y2),"LineWidth",3,"DisplayName","abs(u)","color",col(2,:));
    h2 = plot(T2,-abs(Y2),"LineWidth",3,"DisplayName","abs(u)","color",col(2,:));
    % set(gca,'Fontsize',fsz);
    % xlabel("Time, t","FontSize",fsz);
    % ylabel("abs(z_1), abs(z_2)","FontSize",fsz);
    % tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    % title(tit,"FontSize",fsz);
    % legend()
    % figname = sprintf("sigma%.2fmu%.2f.eps",par(j,1),par(j,2)); 
    % 
    % fig = figure;
    % hold on;
    h3 = plot(T2,real(Y2),"LineWidth",1,"DisplayName","Re(u)","color",col(2,:));
    set(gca,'Fontsize',fsz);
    xlabel("Time, t","FontSize",fsz);
    ylabel("u","FontSize",fsz);
    tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    title(tit,"FontSize",fsz);
    legend([h1,h3])
    figname = sprintf("Torus_sigma%.2fmu%.2f.eps",par(j,1),par(j,2)); 
    % saveas(fig,figname,"epsc");

    fig = figure;
    hold on;
    h1 = plot(T2,angle(Y2),"LineWidth",3,"DisplayName","angle(u)","color",col(2,:));
    set(gca,'Fontsize',fsz);
    xlabel("Time, t","FontSize",fsz);
    ylabel("angle(u)","FontSize",fsz);
    tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    title(tit,"FontSize",fsz);
    %
    [T2,Y2,te,ye,ie] = ode45(odefun,[Tskip,2*Tskip],y0,options1);
    fig = figure;
    plot(real(ye),imag(ye),".","DisplayName","u","color",col(2,:))
    set(gca,'Fontsize',fsz);
    xlabel("Re(u)","FontSize",fsz);
    ylabel("Im(u)","FontSize",fsz);
    tit = strcat("\sigma = ",sprintf("%.2f",par(j,1)),", \mu = ",sprintf("%.2f",par(j,2)));
    title(tit,"FontSize",fsz);
    daspect([1,1,1])
    %
    fig = figure;
    % if j > 8
        [T2,Y2] = ode45(odefun,[Tskip,2*Tskip],y0,options);
    % end
    plot(real(Y2),imag(Y2),"LineWidth",3,"DisplayName","u","color",col(2,:));
    set(gca,'Fontsize',fsz);
    xlabel("Re(u)","FontSize",fsz);
    ylabel("Im(u)","FontSize",fsz);
    title(tit,"FontSize",fsz);
    figname = sprintf("sigma%.2fmu%.2f_phase.eps",par(j,1),par(j,2));
    daspect([1,1,1])
    % saveas(fig,figname,"epsc");

end
end
%% compute Lyapunov's exponents
if flag == 2
for j = 1:npar
    odefun = @(t,y)func(par(j,1),par(j,2),y);
    [T1,Y1] = ode45(odefun,[0,Tskip],[1,0],options);

    % fig = figure;
    % hold on;
    % plot(T1,real(Y1(:,1)),"LineWidth",3,"DisplayName","x1");
    % plot(T1,imag(Y1(:,1)),"LineWidth",3,"DisplayName","y1");
    % plot(T1,real(Y1(:,2)),"LineWidth",3,"DisplayName","x2");
    % plot(T1,imag(Y1(:,2)),"LineWidth",3,"DisplayName","y2");

    y0 = Y1(end,:);
    odefun_extended = @(t,y)SL_extended(y,par(j,1),par(j,2));
% % Initial conditions for the perturbation matrix (identity matrix)
    P0 = eye(4);
% 
% % Combine initial conditions
    X0 = [real(Y1(end,1));
        imag(Y1(end,1));
        real(Y1(end,2));
        imag(Y1(end,2));
    P0(:)];
% 
% % Time span for integration
    [T, X] = ode45(odefun_extended, [Tskip,Tskip + 5*Tplot], X0);
% 
% % Initialize Lyapunov exponents
    LEs = zeros(1, 4);
% 
% % Orthogonalization interval (adjust as needed)
    ortho_interval = 0.1; % Time interval for re-orthonormalization
% 
% % Loop for computing Lyapunov exponents
    num_steps = length(T);
    current_P = reshape(X(1, 5:end), 4,4);


    for i = 2:num_steps
%     % Get the current perturbation matrix
        next_P = reshape(X(i, 5:end), 4,4);
% 
%     % If enough time has passed for re-orthonormalization
        if T(i) - T(i-1) >= ortho_interval || i == num_steps
        % Perform QR decomposition
            [Q, R] = qr(next_P);

        % Update Lyapunov exponents
            LEs = LEs + log(abs(diag(R)))';
        % Reset perturbation matrix for next interval
            current_P = Q;
            X(i, 5:end) = current_P(:); % Update the state vector with the orthonormalized P
    end
end
% 
% % Normalize Lyapunov exponents by total time
    LEs = LEs / (T(end) - T(1));
% 
    disp('Lyapunov Exponents:');
    disp(LEs);
% 


end
end
end
%%
function [position,isterminal,direction] = events(t,y)
  position = mod(t,2*pi)-pi; % The value that we want to be zero
  isterminal = 0;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

%%
function dX = SL_extended(X,sig,mu)
om1 = 1;
om2 = om1 + sig;
lam = 1;
    x1 = X(1);
    y1 = X(2);
    x2 = X(3);
    y2 = X(4);


    % Variational matrix (Jacobian of the system)
    J = [mu - 3*x1^2 - y1^2, -om1-2*y1*x1, 0,0;
         om1 - 2*y1*x1, mu-3*y1^2 - x1^2, 0, 0;
         -lam, 0, mu - 3*x1^2 - y1^2, -om2-2*y1*x2;
         0, -lam, om2-2*x2*y2, mu-3*y2^2 - x2^2
         ];

    % Extract the perturbation matrix from the extended state vector
    % The extended state vector X will be [x; y; z; P11; P12; P13; P21; ... P33]
    % where P is the 3x3 perturbation matrix
    P = reshape(X(5:end), 4, 4);

    % System equations
    a1 = x1^2 + y1^2;
    a2 = x2^2 + y2^2;
    dx1 = mu*x1 - om1*y1 - a1*x1;
    dy1 = mu*y1 + om1*x1 - a1*y1;
    dx2 = mu*x2 - om2*y2 - a2*x2 - lam*x1;
    dy2 = mu*y2 + om2*x2 - a2*y2 - lam*y1;

    % Variational equation: dP/dt = J * P
    dP = J * P;

    % Combine into the output vector
    dX = [dx1; dy1; dx2; dy2; dP(:)];
end

% % Main script to compute Lyapunov exponents
% 
% % Initial conditions for the system
% x0 = [0.1; 0.1; 0.1];
% 
% % Initial conditions for the perturbation matrix (identity matrix)
% P0 = eye(3);
% 
% % Combine initial conditions
% X0 = [x0; P0(:)];
% 
% % Time span for integration
% tspan = [0 100];
% 
% % Solve the extended system of ODEs
% [t, X] = ode45(@lorenz_extended, tspan, X0);
% 
% % Initialize Lyapunov exponents
% LEs = zeros(1, 3);
% 
% % Orthogonalization interval (adjust as needed)
% ortho_interval = 0.1; % Time interval for re-orthonormalization
% 
% % Loop for computing Lyapunov exponents
% num_steps = length(t);
% current_P = reshape(X(1, 4:12), 3, 3);
% for i = 2:num_steps
%     % Get the current perturbation matrix
%     next_P = reshape(X(i, 4:12), 3, 3);
% 
%     % If enough time has passed for re-orthonormalization
%     if t(i) - t(i-1) >= ortho_interval || i == num_steps
%         % Perform QR decomposition
%         [Q, R] = qr(next_P);
% 
%         % Update Lyapunov exponents
%         LEs = LEs + log(abs(diag(R)))';
% 
%         % Reset perturbation matrix for next interval
%         current_P = Q;
%         X(i, 4:12) = current_P(:); % Update the state vector with the orthonormalized P
%     end
% end
% 
% % Normalize Lyapunov exponents by total time
% LEs = LEs / (t(end) - t(1));
% 
% disp('Lyapunov Exponents:');
% disp(LEs);
% 

