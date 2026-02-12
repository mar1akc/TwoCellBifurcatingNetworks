function system3_reduced()
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
    
    figure(1); clf; hold on
    grid
       
    % plot solutions with |v|^2 < 1
    v2vals = 0.5 : 0.08 : 0.98; % used to be 0.08
    Nv = length(v2vals);
    col = 0.7*autumn(Nv);
    for j = 1 : Nv
        v2 = v2vals(j);
        [r1,r2] = get_radius(v2);
        ang = get_angle(v2);
        t = linspace(0,pi/2,Nt);
        [r1,r2] = get_radius(v2);
        x = r1*cos(t);
        y = r2*sin(t);
        dJ = detJ(y,x,v2);
        tJ = trJ(y,x,v2);
        ind = find(dJ > 0 & tJ < 0);
        plot(r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt1, 'color',col(j,:));
        plot(-r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt1, 'color',col(j,:));
        ind = find((dJ > 0 & tJ > 0) | dJ < 0);
        plot(r1*cos(t(ind)),r2*sin(t(ind)),'--','LineWidth', lwt2, 'color',col(j,:));
        if j <=3 
             t1 = t(round(j*Nt/25));
        else
             t1 = t(round((9-j)*Nt/20));
        end
        text(r1*cos(t1),abs(r2)*sin(t1),sprintf("%.3f",v2),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    end
    for j = 1 : Nv
        v2 = v2vals(j);
        [r1,r2] = get_radius(v2);
        ang = get_angle(v2);
        t = linspace(pi/2,pi,Nt);
        [r1,r2] = get_radius(v2);
        x = r1*cos(t);
        y = r2*sin(t);
        dJ = detJ(y,x,v2);
        tJ = trJ(y,x,v2);
        ind = find(dJ > 0 & tJ < 0);
        plot(r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt1, 'color',col(j,:));
        plot(-r1*cos(t(ind)),r2*sin(t(ind)),'LineWidth', lwt1, 'color',col(j,:));
        ind = find((dJ > 0 & tJ >= 0) | dJ <= 0);
        plot(r1*cos(t(ind)),r2*sin(t(ind)),'--','LineWidth', lwt2, 'color',col(j,:));
        % if j <=3 
        %      t1 = t(round(j*Nt/25));
        % else
        %      t1 = t(round((9-j)*Nt/20));
        % end
        % text(r1*cos(t1),abs(r2)*sin(t1),sprintf("%.3f",v2),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    end

        % plot solutions with |v|^2 < 1
    v2vals = 0.05 : 0.1 : 0.45;
    Nv = length(v2vals);
    col = 0.7*gray(Nv);
    for j = 1 : Nv
        v2 = v2vals(j);
        [r1,r2] = get_radius(v2);
        t = linspace(0,pi,Nt);
        plot(r1*cos(t),r2*sin(t),'--','LineWidth', lwt2, 'color',col(j,:));
        t1 = t(round((j+8.5)*Nt/20));
        text(r1*cos(t1),abs(r2)*sin(t1),sprintf("%.3f",v2),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    end

    
    % plot |v|^2 = 1
    plot([-1,-1],[0,ymax],'LineWidth', lwt1, 'color','k');
    plot([1,1],[0,ymax],'LineWidth', lwt1, 'color','k');
    
    % plot solutions with |v|^2 > 1
    v2vals = exp(0.01 : 0.02 : 0.5).^4;
    Nv = length(v2vals);
    col = winter(Nv);
    for j = 1 : Nv
        v2 = v2vals(j);
        [r1,r2] = get_radius(v2);
        t = linspace(0,pi,Nt);
        plot(r1*cos(t),abs(r2)*sin(t),'LineWidth', lwt1, 'color',col(j,:));
        if v2 < 2.7
            t1 = t(round(Nt/2));
            text(r1*cos(t1),abs(r2)*sin(t1),sprintf("%.3f",v2),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        end
    end
    
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
    % ind = find(axy > ang & axy < pi - ang);
    % plot(r1*cos(t(ind)),r2*sin(t(ind)),'--','LineWidth', lwt2, 'color','r');
    
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
    
    plot(bdry(:,1),bdry(:,2),'LineWidth',lwt_bdry,'color','r');
    plot(-bdry(:,1),bdry(:,2),'LineWidth',lwt_bdry,'color','r');

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

    % axis([xmin,xmax,ymin,ymax])

    % plot solutions with |v|^2 < 1
    v2vals = 0.1 : 0.1 : 3; % used to be 0.08
    Nv = length(v2vals);
    col = 0.7*autumn(Nv);
    for j = 1 : Nv
        v2 = v2vals(j);
        ang = get_angle(v2);
        t = linspace(0,pi,Nt);
        [r1,r2] = get_radius_minus(v2);
        plot(r1*cos(t),-r2*sin(t),'LineWidth', lwt1, 'color',col(j,:));
        if v2 < 0.7
            t1 = t(round(Nt/2));
            text(r1*cos(t1),-abs(r2)*sin(t1),sprintf("%.3f",v2),...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        end
    end
    plot([xmin,xmax],[0,0],'Linewidth',lwt_bdry,'Color','r')



    axis([xmin,xmax,-3,ymax])
    % axis([0.99,1.08,1.8,2.05])
%%
% 
figure(2); clf; hold on; grid on;

    v2vals = 0.1 : 0.1 : 3; % used to be 0.08
    Nv = length(v2vals);
    col = 0.7*autumn(Nv);
    for j = 1 : Nv
        v2 = v2vals(j);
        ang = get_angle(v2);
        t = linspace(0,pi,Nt);
        [r1,r2] = get_radius_zero(v2);
        plot(r1*cos(t),r2*sin(t),'LineWidth', lwt1, 'color',col(j,:));
        if v2 < 1.7
            t1 = t(round(Nt/2));
            text(r1*cos(t1),abs(r2)*sin(t1),sprintf("%.3f",v2),...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        end
    end
    set(gca,'Fontsize',fsz);
    axis([xmin,xmax,ymin,ymax])

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
function [r1,r2] = get_radius_minus(v2)
    r2 = 1./(sqrt(v2).*abs((1 + v2)));
    r1 = 1./sqrt(v2);
end
%%
function [r1,r2] = get_radius_zero(v2)
    r2 = 1./sqrt(v2.^3);
    r1 = 1./sqrt(v2);
end