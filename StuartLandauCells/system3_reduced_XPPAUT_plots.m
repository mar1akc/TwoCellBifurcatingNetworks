function system3_reduced_XPPAUT_plots()
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

    %% patches

    % pink
    v2 = linspace(0.5,0.75,100);
    s = sqrt((3*v2 - 1)./(2*v2.^2));
    m = sqrt(0.5./(v2.^2.*(1-v2)));

    mvals = linspace(2,4*sqrt(2)/3,100);
    svals = sqrt(2-0.25*mvals.^2);

    verts_pink = [[s,svals]',[m,mvals]'];
    Nverts_pink = size(verts_pink,1);

    patch('Vertices', verts_pink, 'Faces', (1:Nverts_pink), ...
      'FaceColor', [0.5,0,0], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');

    verts_pink = [-[s,svals]',[m,mvals]'];
    Nverts_pink = size(verts_pink,1);

    patch('Vertices', verts_pink, 'Faces', (1:Nverts_pink), ...
      'FaceColor', [0.5,0,0], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');


    % yellow
    mvals1 = linspace(0,4*sqrt(2)/3,100);
    svals1 = real(sqrt(2-0.25*mvals1.^2));

    v2 = linspace(1/3,0.75,1000);
    s = sqrt((3*v2 - 1)./(2*v2.^2));
    m = sqrt(0.5./(v2.^2.*(1-v2)));

    verts_yellow = [[svals1,s(end:-1:1),0]',[mvals1,m(end:-1:1),0]'];
    Nverts_yellow = size(verts_yellow,1);
    patch('Vertices', verts_yellow, 'Faces', (1:Nverts_yellow), ...
      'FaceColor', [1,252/255,121/255], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');
    verts_yellow = [-[svals1,s(end:-1:1),0]',[mvals1,m(end:-1:1),0]'];
    patch('Vertices', verts_yellow, 'Faces', (1:Nverts_yellow), ...
      'FaceColor', [1,252/255,121/255], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');

    % green
    v2 = linspace(1/3,0.5,100);
    s = sqrt((3*v2 - 1)./(2*v2.^2));
    m = sqrt(0.5./(v2.^2.*(1-v2)));
    
    v2 = linspace(0.75,0.999,100);
    s2 = sqrt((3*v2 - 1)./(2*v2.^2));
    m2 = sqrt(0.5./(v2.^2.*(1-v2)));

    verts_green = [[s,svals,s2,0]',[m,mvals,m2,ymax]'];
    Nverts_green = size(verts_green,1);
    patch('Vertices', verts_green, 'Faces', (1:Nverts_green), ...
      'FaceColor', [0,143/255,0], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');
    verts_green = [-[s,svals,s2,0]',[m,mvals,m2,ymax]'];
    Nverts_green = size(verts_green,1);
    patch('Vertices', verts_green, 'Faces', (1:Nverts_green), ...
      'FaceColor', [0,143/255,0], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');

    % blue
    verts_blue = [[svals1,s2,xmax,xmax]',[mvals1,m2,ymax,ymin]'];
    Nverts_blue = size(verts_blue,1);
    patch('Vertices', verts_blue, 'Faces', (1:Nverts_blue), ...
      'FaceColor', [118/255,214/255,1], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');
    verts_blue = [-[svals1,s2,xmax,xmax]',[mvals1,m2,ymax,ymin]'];
    Nverts_blue = size(verts_blue,1);
    patch('Vertices', verts_blue, 'Faces', (1:Nverts_blue), ...
      'FaceColor', [118/255,214/255,1], ...
      'FaceAlpha', 0.2, ...
      'EdgeColor', 'none');
    

%%


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

    axis([0.8,xmax,ymin,ymax])

    %%
    ep = 0.2;
    s = 0.98;
    lam = 1;
    Nmu = 1000;
    muvals = linspace(1e-3,5,Nmu);
    mu = (((muvals+ep).^1.5)./sqrt(muvals))/lam;
    sigma = (s/lam)*sqrt((muvals+ep)./muvals);
    plot(sigma,mu,'Linewidth',lwt_bdry,'color','k');

    rectangle('Position',[1.01,2.5,0.02,0.6],'Edgecolor','k','LineWidth',0.5);
    rectangle('Position',[0.99,1.8,0.09,0.25],'Edgecolor','k','LineWidth',0.5);
    


    % axis([1.01,1.03,2.5,3.1])
    % axis([0.99,1.08,1.8,2.05])
    % axis([0.95,1.5,ymin,ymax])

end

%% Compute the values of \mu and \tilde{mu} at singular points
s = 0.98;
e = 0.2;
l = 1;
f = @(x)(x+e)^3/(8*l^2*x) + s^2*(x+e)/(2*l^2*x) - 1;
x0 = 1.4;
mu = fzero(f,x0);
fprintf("mu = %d, mu_tilde = %d\n",mu,(mu+e)^1.5/(l*sqrt(mu)));
x0 = 0.5;
mu = fzero(f,x0);
fprintf("mu = %d, mu_tilde = %d\n",mu,(mu+e)^1.5/(l*sqrt(mu)));
%%
f = @(x) fun(x,s,e,l);
x0 = 7/12;
x = fzero(f,x0);
mu = s/sqrt((1-x)*(3*x-1)) - e;
fprintf("mu = %d, mu_tilde = %d\n",mu,(mu+e)^1.5/(l*sqrt(mu)));

x0 = 5/6;
x = fzero(f,x0);
mu = s/sqrt((1-x)*(3*x-1)) - e;
fprintf("mu = %d, mu_tilde = %d\n",mu,(mu+e)^1.5/(l*sqrt(mu)));

x0 = 0.95;
x = fzero(f,x0);
mu = s/sqrt((1-x)*(3*x-1)) - e;
fprintf("mu = %d, mu_tilde = %d\n",mu,(mu+e)^1.5/(l*sqrt(mu)));

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

%%
%%
function y = fun(x,s,e,l)
mu = s/sqrt((1-x)*(3*x-1)) - e;
mue = mu+e;
y = mue^3*x*(1-x)^2 + s^2*mue*x - l^2*mu;
end
