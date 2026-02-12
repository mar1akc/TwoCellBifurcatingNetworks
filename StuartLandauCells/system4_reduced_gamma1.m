function system4_reduced_gamma1()
gam = 1;
Nt = 400;
t = linspace(0,2*pi,2000);

detJfun = @(v2,s,m,g) m.^2*(1-v2).*(1-3*v2) + (s-m.*g.*v2).*(s-3*m.*g.*v2);
trJfun = @(v2,m) 2*m.*(1-2*v2);

% % check that the det and tr are correct
% for k = 1 : 10
%     vr = rand;
%     vi = rand;
%     m = rand;
%     s = rand;
%     g = gam;
%     v2rand = vr^2 + vi^2;
%     J = Jac(vr,vi,s,m,g);
%     fprintf("det J = %d, det Jfun = %d\n", det(J),detJfun(v2rand,s,m,g) );
%     fprintf("tr J = %d, tr Jfun = %d\n", trace(J),trJfun(v2rand,m) );
% end
% 
%%
figure(1); clf; hold on;
smin = -1.5; 
smax = 5;
mmin = 0;
mmax = 5;
lwt1 = 1.5;
lwt2 = 0.5;
lwt_bdry = 3;

v2 = (0.05:0.05:0.45);

Nv2 = length(v2);
col = 0.7*gray(Nv2);

%% small v2
% [ss,mm] = meshgrid(linspace(smin,smax,50),linspace(mmin,mmax,50));
% elfun = @(s,m,v)m.^2*(1-v).^2*v + (s-m*gam*v).^2*v;
for j = 1 : Nv2
    xy = par_ellipse(v2(j),gam,t);
    plot(xy(1,:),xy(2,:),'LineWidth',lwt2,'color',col(j,:),'Linestyle','--');
    % contour(ss,mm,elfun(ss,mm,v2(j)),[1,1],'LineWidth',lwt2,'color','r','Linestyle','--');
    i = 800 - 20*j;
    text(xy(1,i),xy(2,i),sprintf("%.3f",v2(j)),'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

end

%% v2 > 0.5

v2 = 0.58 : 0.08 : 0.98;

Nv2 = length(v2);
col = 0.7*autumn(Nv2);

for j = 1 : Nv2
    xy = par_ellipse(v2(j),gam,t);
    detJ = detJfun(v2(j),xy(1,:),xy(2,:),gam);
    ind = find(detJ > 0 & afun(xy(1,:),xy(2,:),gam,v2(j)) > 0 & xy(2,:) > 0);
    plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',col(j,:));
    ind = find(detJ > 0 & afun(xy(1,:),xy(2,:),gam,v2(j)) < 0 & xy(2,:) > 0);
    plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',col(j,:));
    ind = find(detJ < 0 & xy(2,:) > 0);
    plot(xy(1,ind),xy(2,ind),'LineWidth',lwt2,'color',col(j,:),'Linestyle','--');
    % plot(xy(1,:),xy(2,:),'LineWidth',lwt2,'color',col(j,:),'Linestyle','--');
    % contour(ss,mm,elfun(ss,mm,v2(j)),[1,1],'LineWidth',lwt2,'color','r','Linestyle','--');
    i = 530 + 10*round((j+1));
    if j == Nv2
        i = 520;
    end
    text(xy(1,i),xy(2,i),sprintf("%.3f",v2(j)),'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

end

%% v2 > 1
v2 = exp(0.01 : 0.02 : 0.5).^4;
Nv2 = length(v2);
col = winter(Nv2);
% [ss,mm] = meshgrid(linspace(smin,smax,50),linspace(mmin,mmax,50));
% elfun = @(s,m,v)m.^2*(1-v).^2*v + (s-m*gam*v).^2*v;

for j = 1 : Nv2
    xy = par_ellipse(v2(j),gam,t);
    detJ = detJfun(v2(j),xy(1,:),xy(2,:),gam);
    ind = find(detJ > 0 & afun(xy(1,:),xy(2,:),gam,v2(j)) > 0 & xy(2,:) > 0);
    plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',col(j,:));
    ind = find(detJ > 0 & afun(xy(1,:),xy(2,:),gam,v2(j)) < 0 & xy(2,:) > 0);
    plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',col(j,:));
    ind = find(detJ < 0 & xy(2,:) > 0);
    plot(xy(1,ind),xy(2,ind),'LineWidth',lwt2,'color',col(j,:),'Linestyle','--');
    % contour(ss,mm,elfun(ss,mm,v2(j)),[1,1],'LineWidth',lwt2,'color','r','Linestyle','--');
    if v2(j) < 2.0
        i = 540+20*j;
        text(xy(1,i),xy(2,i),sprintf("%.3f",v2(j)),'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    end
    % contour(ss,mm,elfun(ss,mm,v2(j)),[1,1],'LineWidth',lwt2,'color','r','Linestyle','--');
end

%% lines

plot([1+mmin*gam,1+mmax*gam],[mmin,mmax],'Linewidth',lwt1,'color','k');
plot([-1+mmin*gam,-1+mmax*gam],[mmin,mmax],'Linewidth',lwt1,'color','k');

%% det J = 0
v2 = linspace(1/3,0.99,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
itermax = 100;
tol = 1e-9;
m0 = 3*sqrt(3)*0.5;
a0 = 0;
s(1) = gam*m0*v2(1);
m(1) = m0;
w = [m0,a0]';
for j = 1:Nv2
    rJ = @(w)Res_Jac(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    s(j) = w(2) + gam*w(1)*v2(j);
    m(j) = w(1);
    fprintf("j = %d, s = %d, m = %d\n",j,s(j),m(j));
end
plot(s,m,'LineWidth',lwt_bdry,'color','m');

m0 = 3*sqrt(3)/(2*sqrt(2));
a0 = 2*gam*m0/3;
w = [m0,a0]';
for j = 1:Nv2
    rJ = @(w)Res_Jac(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    s(j) = w(2) + gam*w(1)*v2(j);
    m(j) = w(1);
    fprintf("j = %d, s = %d, m = %d\n",j,s(j),m(j));
end
plot(s,m,'LineWidth',lwt_bdry,'color','m');

v2 = linspace(1.01,2,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
for j = 1:Nv2
    rJ = @(w)Res_Jac(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    s(j) = w(2) + gam*w(1)*v2(j);
    m(j) = w(1);
    fprintf("j = %d, s = %d, m = %d\n",j,s(j),m(j));
end
plot(s,m,'LineWidth',lwt_bdry,'color','m');

v2 = linspace(1.01,2,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
m0 = 3.73;
s0 = 4.93;
w = [m0,s0]';

for j = 1:Nv2
    rJ = @(w)Res_Jac(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    s(j) = w(2) + gam*w(1)*v2(j);
    m(j) = w(1);
    fprintf("j = %d, s = %d, m = %d\n",j,s(j),m(j));
end
plot(s,m,'LineWidth',lwt_bdry,'color','m');

%%
%% v2 = 0.5
xy = par_ellipse(0.5,gam,t);
detJ = detJfun(0.5,xy(1,:),xy(2,:),gam);
ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color','r');
ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color','r');
ind = find(detJ < 0 & xy(2,:) > 0);
plot(xy(1,ind),xy(2,ind),'LineWidth',lwt2,'color','r','Linestyle','--');

%% three roots boundary
v2 = linspace(1/3-4.05e-2,1/3,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
m0 = 2.41458;
s0 = 1.41434;
w = [m0,s0]';
itermax = 100;
for j = 1:Nv2
    rJ = @(w)Res_Jac_3sol(w,gam,v2(j));
    [w,iter] = LevenbergMarquardt_iter(rJ,w,itermax,tol);
    if iter < itermax && abs(detJfun(v2(j),s(j),m(j),gam)) < 1e-9
        s(j) = w(2);
        m(j) = w(1);
    else
        s(j) = 0;
        m(j) = 0;
    end
    fprintf("j = %d, v2 = %d, s = %d, m = %d, det = %d\n",...
        j,v2(j),s(j),m(j),detJfun(v2(j),s(j),m(j),gam));
end
ind = find(m > 0);
plot(s(ind),m(ind),'LineWidth',lwt_bdry,'color','m');
plot(s(ind),m(ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');

sstar = s(1);
mstar = m(1);


v2 = linspace(1/3,0.99,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
itermax = 100;
for j = 1:Nv2
    rJ = @(w)Res_Jac_3sol(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    if iter < itermax && abs(detJfun(v2(j),s(j),m(j),gam)) < 1e-9
        s(j) = w(2);
        m(j) = w(1);
    else
        s(j) = 0;
        m(j) = 0;
    end
    fprintf("j = %d, v2 = %d, s = %d, m = %d, det = %d\n",...
        j,v2(j),s(j),m(j),detJfun(v2(j),s(j),m(j),gam));
end
ind = find(m > 0);
plot(s(ind),m(ind),'LineWidth',lwt_bdry,'color','m');
plot(s(ind),m(ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');
%% three roots boundary
%v2start = (2/3)*(1+gam^2)/(1-gam/3);
% 
v2 = linspace(1/3-4.e-2,1/3,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
m0 = 2.392;
s0 = 1.514;
w = [m0,s0]';
itermax = 100;
for j = 1:Nv2
    rJ = @(w)Res_Jac_3sol(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    s(j) = w(2);
    m(j) = w(1);
    if iter < itermax && abs(detJfun(v2(j),s(j),m(j),gam)) < 1e-9
        s(j) = w(2);
        m(j) = w(1);
    else
        s(j) = 0;
        m(j) = 0;
    end
    fprintf("j = %d, v2 = %d, s = %d, m = %d, det = %d\n",...
        j,v2(j),s(j),m(j),detJfun(v2(j),s(j),m(j),gam));
end
m = [mstar;m];
s = [sstar;s];
ind = find(m > 0);
plot(s(ind),m(ind),'LineWidth',lwt_bdry,'color','m');
plot(s(ind),m(ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');


%% three roots boundary
v2 = linspace(1/3,1,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
itermax = 100;
for j = 1:Nv2
    rJ = @(w)Res_Jac_3sol(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    if iter < itermax && abs(detJfun(v2(j),s(j),m(j),gam)) < 1e-9
        s(j) = w(2);
        m(j) = w(1);
    else
        s(j) = 0;
        m(j) = 0;
    end
    fprintf("j = %d, v2 = %d, s = %d, m = %d, det = %d\n",...
        j,v2(j),s(j),m(j),detJfun(v2(j),s(j),m(j),gam));
end
ind = find(m > 0);
plot(s(ind),m(ind),'LineWidth',lwt_bdry,'color','m');
plot(s(ind),m(ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');

%% three roots boundary
v2 = linspace(1.1,2,100);
Nv2 = length(v2);
s = zeros(Nv2,1);
m = zeros(Nv2,1);
m0 = 3.73;
s0 = 4.93;
w = [m0,s0]';
itermax = 100;
for j = 1:Nv2
    rJ = @(w)Res_Jac_3sol(w,gam,v2(j));
    w = LevenbergMarquardt(rJ,w,itermax,tol);
    if iter < itermax && abs(detJfun(v2(j),s(j),m(j),gam)) < 1e-9
        s(j) = w(2);
        m(j) = w(1);
    else
        s(j) = 0;
        m(j) = 0;
    end
    fprintf("j = %d, v2 = %d, s = %d, m = %d, det = %d\n",...
        j,v2(j),s(j),m(j),detJfun(v2(j),s(j),m(j),gam));
end
ind = find(m > 0);
plot(s(ind),m(ind),'LineWidth',lwt_bdry,'color','m');
plot(s(ind),m(ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');


fsz = 20;
set(gca,'Fontsize',fsz);
grid;
axis([smin,smax,mmin,mmax])
end
%%
function J = Jac(vr,vi,s,m,g)
J = zeros(2);
v2 = vr.^2 + vi.^2;
J(1,1) = m.*(1-v2) - 2*m.*(vr.^2 - g*vr.*vi);
J(1,2) = -(s-m.*g.*v2) - 2*m.*(vr.*vi - g*vi.^2);
J(2,1) = (s-m.*g.*v2) - 2*m.*(vr.*vi + g*vr.^2);
J(2,2) = m.*(1-v2) - 2*m.*(vi.^2 + g*vr.*vi);
end

%%
function M = Mtx(v2,g)
M = zeros(2);
M(1,1) = v2;
M(1,2) = -g*v2.^2;
M(2,1) = M(1,2);
M(2,2) = (1-v2).^2*v2 + g*v2.^3;
end
%%
function xy = par_ellipse(v2,g,t)
    M = Mtx(v2,g);
    [V,E] = eig(M);
    xy = V*[(1/sqrt(E(1,1)))*cos(t);(1/sqrt(E(2,2)))*sin(t)];
end
%%
function f = afun(s,m,g,v2)
f = s - g*m*v2;
end
%%
function [r,J] = Res_Jac(ma,g,v2)
m = ma(1);
a = ma(2);
r = zeros(2,1);
J = zeros(2);
r(1) = m^2*(1-v2)^2 + a^2 - 1/v2;
r(2) = m^2*(1-v2)*(1-3*v2) + a^2 - 2*a*g*m*v2;
J(1,1) = 2*m*(1-v2)^2;
J(1,2) = 2*a;
J(2,1) = 2*m*(1-v2)*(1-3*v2) - 2*a*g*v2;
J(2,2) = 2*a - 2*m*g*v2;
end

%%
%%
function [r,J] = Res_Jac_3sol(ma,g,v2)
m = ma(1);
s = ma(2);
a = s - g*m*v2;
r = zeros(2,1);
J = zeros(2);
r(1) = m^2*(1-v2)^2 + a^2 - 1/v2;
r(2) = 2*m^2*(1-v2) + 2*a*g*m - 1/(v2^2);
J(1,1) = 2*m*(1-v2)^2 - 2*a*g*v2;
J(1,2) = 2*a;
J(2,1) = 4*m*(1-v2) + 2*a*g - 2*g^2*m*v2;
J(2,2) = 2*m*g;
end



%%
function w = LevenbergMarquardt(Res_and_Jac,w,iter_max,tol)
fsz = 16; % fontsize
%%
Rmax = 1e-0;
Rmin = 1e-14;
rho_good = 0.75;
rho_bad = 0.25;
eta = 0.01;
% iter_max = 10000;
% tol = 5e-3;
%% setup training mesh
[r,J] = Res_and_Jac(w);
f = F(r);
g = J'*r;
nor = norm(g);
R = Rmax/5; % initial trust region radius

% fprintf('Initially: f = %d, nor(g) = %d\n',f,nor); 
%%
tic

iter = 0;
flag = 1;
I = eye(length(w));
% quadratic model: m(p) = (1/2)||r||^2 + p'*J'*r + (1/2)*p'*J'*J*p;
norg = zeros(iter_max+1,0);
fall = zeros(iter_max+1,0);
norg(1) = nor;
fall(1) = f;
rho = 1;
while nor > tol && iter < iter_max
    % solve the constrained minimization problem using dogleg strategy
    B = J'*J + (1e-6)*I;
    pstar = -B\g;
    % fprintf('iter %d: ',iter);
    if norm(pstar) <= R
        p = pstar;
        % fprintf('Global min of quad model\n');
    else % solve constrained minimization problem
        lam = 1;
        isub = 0;
        while 1
            B1 = B + lam*I;
            C = chol(B1);
            p = -C\(C'\g);
            np = norm(p);
            dd = abs(np - R);
            if dd < 1e-6
                break
            end
            q = C'\p;
            nq = norm(q);
            lamnew = lam +(np/nq)^2*(np - R)/R;
            if lamnew < 0
                lam = 0.5*lam;
            else
                lam = lamnew;
            end
            isub = isub + 1;
        end
        % fprintf('Contraint minimization: %d substeps\n',isub);
    end
    iter = iter + 1;  
    if flag == 0
        break;
    end
    % assess the progress
    wnew = w + p;
    [rnew, Jnew] = Res_and_Jac(wnew);
    mnew = 0.5*r'*r + g'*p + 0.5*p'*B*p;
    fnew = F(rnew);
    rho = (f - fnew + 1e-14)/(f - mnew + 1e-14);
    

    % adjust the trust region
    if rho < rho_bad
        R = max([0.25*R,Rmin]);
    else
        if rho > rho_good && abs(norm(p) - R) < tol
            R = min([Rmax,2*R]);
        end
    end
    % accept or reject step
    if rho > eta            
        w = wnew;
        r = rnew;
        J = Jnew;
        f = fnew;
        g = J'*r;
        nor = norm(g);        
         % fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
    end
    norg(iter+1) = nor;
    fall(iter+1) = f;
end
% fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
cputime = toc;
% fprintf('CPUtime = %d, iter = %d\n',cputime,iter);
end

%%
function p = cauchy_point(B,g,R)
    ng = norm(g);
    ps = -g*R/ng;
    aux = g'*B*g;
    if aux <= 0
        p = ps;
    else
        p = min(ng^3/(R*aux),1);
    end
end
%%
function f = F(r)
    f = 0.5*r'*r;
end

%%
%%
function [w,iter] = LevenbergMarquardt_iter(Res_and_Jac,w,iter_max,tol)
fsz = 16; % fontsize
%%
Rmax = 1e-0;
Rmin = 1e-14;
rho_good = 0.75;
rho_bad = 0.25;
eta = 0.01;
% iter_max = 10000;
% tol = 5e-3;
%% setup training mesh
[r,J] = Res_and_Jac(w);
f = F(r);
g = J'*r;
nor = norm(g);
R = Rmax/5; % initial trust region radius

% fprintf('Initially: f = %d, nor(g) = %d\n',f,nor); 
%%
tic

iter = 0;
flag = 1;
I = eye(length(w));
% quadratic model: m(p) = (1/2)||r||^2 + p'*J'*r + (1/2)*p'*J'*J*p;
norg = zeros(iter_max+1,0);
fall = zeros(iter_max+1,0);
norg(1) = nor;
fall(1) = f;
rho = 1;
while nor > tol && iter < iter_max
    % solve the constrained minimization problem using dogleg strategy
    B = J'*J + (1e-6)*I;
    pstar = -B\g;
    % fprintf('iter %d: ',iter);
    if norm(pstar) <= R
        p = pstar;
        % fprintf('Global min of quad model\n');
    else % solve constrained minimization problem
        lam = 1;
        isub = 0;
        while 1
            B1 = B + lam*I;
            C = chol(B1);
            p = -C\(C'\g);
            np = norm(p);
            dd = abs(np - R);
            if dd < 1e-6
                break
            end
            q = C'\p;
            nq = norm(q);
            lamnew = lam +(np/nq)^2*(np - R)/R;
            if lamnew < 0
                lam = 0.5*lam;
            else
                lam = lamnew;
            end
            isub = isub + 1;
        end
        % fprintf('Contraint minimization: %d substeps\n',isub);
    end
    iter = iter + 1;  
    if flag == 0
        break;
    end
    % assess the progress
    wnew = w + p;
    [rnew, Jnew] = Res_and_Jac(wnew);
    mnew = 0.5*r'*r + g'*p + 0.5*p'*B*p;
    fnew = F(rnew);
    rho = (f - fnew + 1e-14)/(f - mnew + 1e-14);
    

    % adjust the trust region
    if rho < rho_bad
        R = max([0.25*R,Rmin]);
    else
        if rho > rho_good && abs(norm(p) - R) < tol
            R = min([Rmax,2*R]);
        end
    end
    % accept or reject step
    if rho > eta            
        w = wnew;
        r = rnew;
        J = Jnew;
        f = fnew;
        g = J'*r;
        nor = norm(g);        
         % fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
    end
    norg(iter+1) = nor;
    fall(iter+1) = f;
end
% fprintf('iter # %d: f = %.14f, |df| = %.4e, rho = %.4e, R = %.4e\n',iter,f,nor,rho,R);
cputime = toc;
% fprintf('CPUtime = %d, iter = %d\n',cputime,iter);
end


% %% v2 = 1/3
% xy = par_ellipse(1/3,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');
% 
% %% v2 = 1/3
% v2 = 1/3-1e-2;
% xy = par_ellipse(v2,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0,0.5,0],'Linestyle','-.');
% 
% v2 = 1/3-2e-2;
% xy = par_ellipse(v2,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0,0.5,0.5],'Linestyle','-.');
% 
% v2 = 1/3-3e-2;
% xy = par_ellipse(v2,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0.5,0.5,0],'Linestyle','-.');
% 
% v2 = 1/3-4e-2;
% xy = par_ellipse(v2,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0.5,0.1,0],'Linestyle','-.');
% 
% v2 = 1/3-5e-2;
% xy = par_ellipse(v2,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0.5,0.1,0.3],'Linestyle','-.');
% 
% v2 = 1/3-6e-2;
% xy = par_ellipse(v2,gam,t);
% % detJ = detJfun(1/3,xy(1,:),xy(2,:),gam);
% % ind = find(detJ > 0 & xy(1,:) > 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% % ind = find(detJ > 0 & xy(1,:) < 1 & xy(2,:) > 0);
% % plot(xy(1,ind),xy(2,ind),'LineWidth',lwt_bdry,'color',[0,0.5,0],'Linestyle','-.');
% ind = find(detJ < 0 & xy(2,:) > 0);
% plot(xy(1,ind),xy(2,ind),'LineWidth',lwt1,'color',[0,0.1,0.7],'Linestyle','-.');
% 
