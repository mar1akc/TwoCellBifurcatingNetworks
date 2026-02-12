function system2_HBpts()
% The reduced equation is
% \dot{v}_R & = \tilde{\mu}(1 - |v|^2)v_R - \tilde{\sigma} v_I - 1,\\
% \dot{v}_I & = \tilde{\mu}(1 - |v|^2)v_I + \tilde{\sigma} v_R.  
% Plot the level sets of |v|^2
% Let a: = 1 - |v|^2
% Then v_R = mu*a/(mu^2*a^2 + s^2), v_I = - s/(mu^2*a^2 + s^2).
% The equation for a is 1 - a = 1/(mu^2*a^2 + s^2), or 
% (mu^2*a^2 + s^2)(1 - a) = 1
% Hence, mu and s satisfying for a given |v|^2 lie on the ellipse
% mu^2*(1 - |v|^2)^2 + s^2 = 1/|v|^2
% Stability conditions: 
% det J = mu^2*(1-|v|^2)*(1 - 3*|v|^2) + s^2 > 0
% tr J = mu*(1 - 2*|v|^2) < 0
% Hence, |v|^2 > 1/2 for stability
% |v|^2 takes all values from 1/2 to +\infty
%
% Stability boundary. Trace. |v|^2 = 0.5.
% Then the ellipse is 0.25*mu^2 + s^2 = 2.
%
% Stability boundary. Determinant. 
% s^2 > mu^2*(3|v|^2 -1)*(1-|v|^2)
%%
fsz = 20;
flag = 0;
% flag = 0 ==> plot phase diagram
% flag = 1 ==> plot bifuration diagrams
s = 0.75;
%================================================

detJ = @(mu,s,v2) mu.^2*(1-v2)*(1 - 3*v2) + s.^2;
trJ = @(mu,s,v2) mu*(1 - 2*v2);

if flag == 0
    Nt = 2000;
    lwt_bdry = 3;
    lwt1 = 1.5;
    lwt2 = 0.5;
    fsz = 20;
    xmin = -1.5;
    xmax = 1.5;
    ymin = 0;
    ymax = 5;
    
    figure(2); clf; hold on
    grid
       
    
        % trace boundary
    v2 = 0.5;
    ang = get_angle(v2);
    t = linspace(0,pi,Nt);
    [r1,r2] = get_radius(v2);
    x = r1*cos(t);
    y = r2*sin(t);
    axy = angle(x + 1i*y);%atan(r1/r2*tan(angle(x + 1i*y))); %
    ind = find(axy < ang);
    plot(r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt1, 'color','r');
    plot(-r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt1, 'color','r');
    ind = find(axy > ang & axy < pi - ang);
    plot(r1*cos(t(ind)),r2*sin(t(ind)),'--','LineWidth', lwt2, 'color','r');
    
    % det boundary
    v2vals = 1-1./linspace(2,100,2000);
    Nev = length(v2vals);
    bdry = zeros(Nev,2);
    for j = 1 : Nev
        v2 = v2vals(j);
        [r1,r2] = get_radius(v2);
        t = linspace(0,pi,Nt);
        x = r1*cos(t);
        y = r2*sin(t);
        axy = atan(r1/r2*tan(angle(x + 1i*y)));%angle(x + 1i*y);
        ang = get_angle(v2);
        tstar = interp1(axy,t,ang);
        bdry(j,:) = [r1*cos(tstar),r2*sin(tstar)];
    end
    
    % plot(bdry(:,1),bdry(:,2),'LineWidth',lwt_bdry,'color',[0.5,0,0]);
    % plot(-bdry(:,1),bdry(:,2),'LineWidth',lwt_bdry,'color',[0.5,0,0]);

    % plot(bdry(:,1),bdry(:,2),'-.','LineWidth',lwt1,'color',[0.5,0,0]);
    % plot(-bdry(:,1),bdry(:,2),'-.','LineWidth',lwt1,'color',[0.5,0,0]);

    % three roots boundary
    v2 = linspace(1/3,0.999,100);
    s = sqrt((3*v2 - 1)./(2*v2.^2));
    m = sqrt(0.5./(v2.^2.*(1-v2)));
    plot(s,m,'-.','LineWidth',lwt1,'color',[0,0.5,0]);
    plot(-s,m,'-.','LineWidth',lwt1,'color',[0,0.5,0]);



    % sync boundary
    v2star = 0.75; %(13+sqrt(41))/32;
    v2vals = 1-1./linspace(1/(1-v2star),100,2000);
    Nev = length(v2vals);
    bdry = zeros(Nev,2);
    for j = 1 : Nev
        v2 = v2vals(j);
        [r1,r2] = get_radius(v2);
        t = linspace(0,pi,Nt);
        x = r1*cos(t);
        y = r2*sin(t);
        axy = angle(x + 1i*y);
        ang = get_angle(v2);
        tstar = interp1(axy,t,ang);
        bdry(j,:) = [r1*cos(tstar),r2*sin(tstar)];
    end
    
    plot(bdry(:,1),bdry(:,2),'LineWidth',lwt_bdry,'color','m');
    plot(-bdry(:,1),bdry(:,2),'LineWidth',lwt_bdry,'color','m');

    v2 = 0.5;
    ang = get_angle(v2star);
    t = linspace(0,pi,Nt);
    [r1,r2] = get_radius(v2);
    x = r1*cos(t);
    y = r2*sin(t);
    axy = angle(x + 1i*y);
    ind = find(axy < ang);
    plot(r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt_bdry, 'color','r');
    plot(-r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt_bdry, 'color','r');
    % ind = find(axy > ang & axy < pi - ang);
    % plot(r1*cos(t(ind)),r2*sin(t(ind)),'--','LineWidth', lwt2, 'color','r');

    
    set(gca,'Fontsize',fsz);
    % xlabel(s,'Fontsize',fsz);
    % ylabel(\tilde{\mu},'Fontsize',fsz);
    % plot([0.99,1.08,1.08,0.99,0.99],[1.8,1.8,2.05,2.05,1.8],'color','k')
%% Bifurcation point
gam = 0;
m = 1.5*sqrt(3);
s = gam*m/3;
plot(s,m,'.','MarkerSize',50,'color','k');
plot(s,m,'.','MarkerSize',30,'color','c');

%% Hysteresis points
a = 1/sqrt(3);
m = 1.5*sqrt(1.5)/(sqrt(1+gam^2))*[(1+gam*a)^1.5,(1-gam*a)^1.5];
s = 1.5*sqrt(1.5)/(sqrt(1+gam^2))*[(gam-a)*sqrt(1+a*gam),(gam+a)*sqrt(1-a*gam)];
plot(s,m,'.','MarkerSize',50,'color','k');
plot(s,m,'.','MarkerSize',30,'color','y');



    axis([xmin,xmax,ymin,ymax])
    % axis([0.99,1.08,1.8,2.05])
end
%% flag == 1
if flag == 1
    mumin = 0;
    mumax = 5;
    vmax = min(3,1/s^2);

    figure(1); clf; hold on
    mu = @(x) (1-s^2*x)./(x.*(1-x).^2);
    v2 = linspace(1.01,vmax,100);
    m = sqrt(mu(v2));
    plot(m,v2,'LineWidth',3,'Color','r');
    v2 = linspace(0.001,0.999,100);
    m = sqrt(mu(v2));

    ind = find(m.^2.*((v2-1).*(3*v2 - 1))+ s^2 <= 0 | v2 <= 0.5);
    if ~isempty(ind)
        plot(m(ind),v2(ind),'--','LineWidth',3,'Color','r');
    end
    ind = find(m.^2.*((v2-1).*(3*v2 - 1))+ s^2 > 0 & v2 > 0.5);
    plot(m(ind),v2(ind),'LineWidth',3,'Color','r');
    axis([mumin,mumax,0,3])
    plot([mumin,mumax],[1,1],'--','color','k')
    set(gca,'Fontsize',fsz);
end
end
%%
function ang = get_angle(v2)
ang = atan(1./sqrt((3*v2 - 1).*(1-v2)));
end
%%
function [r1,r2] = get_radius(v2)
    r2 = 1./(sqrt(v2).*abs((1 - v2)));
    r1 = 1./sqrt(v2);
end

%%
%%
function [r1,r2] = get_radius_minus(v2)
    r2 = 1./(sqrt(v2).*abs((1 + v2)));
    r1 = 1./sqrt(v2);
end