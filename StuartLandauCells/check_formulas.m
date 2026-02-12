function check_formulas()
G = @(x,m,e,s,g,l)x^3*(1+g^2) - 2*(m+e+s*g)*x^2 +((m+e)^2 + s^2)*x - l^2*m;
Gx = @(x,m,e,s,g,l)3*x^2*(1+g^2) - 4*(m+e+s*g)*x + (m+e)^2 + s^2;
Gxx = @(x,m,e,s,g,l)6*x*(1+g^2) - 4*(m+e+s*g);

r3 = 1/sqrt(3);
ep1 = @(m,g,l)(1.5*(1 + r3*g)*l^(2/3)*m^(1/3))/((1+g^2)^(1/3)) - m;
ep2 = @(m,g,l)(1.5*(1 - r3*g)*l^(2/3)*m^(1/3))/((1+g^2)^(1/3)) - m;
sig1 = @(m,g,l)(1.5*(g - r3)*l^(2/3)*m^(1/3))/((1+g^2)^(1/3));
sig2 = @(m,g,l)(1.5*(g + r3)*l^(2/3)*m^(1/3))/((1+g^2)^(1/3));
u2 = @(m,g,l)l^(2/3)*m^(1/3)*(1+g^2)^(-1/3);

for j = 1 : 10
    if j == 1
        g = 0;
    else
        if j == 2
            g = 1;
        else
            if j == 3
                g = -1;
            else
                g = sqrt(3)*(rand);
            end
        end
    end
    l = rand;
    m = 3*rand;

    e1 = ep1(m,g,l);
    e2 = ep2(m,g,l);
    s1 = sig1(m,g,l);
    s2 = sig2(m,g,l);
    x = u2(m,g,l);
    f1 = G(x,m,e1,s1,g,l);
    fx1 = Gx(x,m,e1,s1,g,l);
    fxx1 = Gxx(x,m,e1,s1,g,l);
    f2 = G(x,m,e2,s2,g,l);
    fx2 = Gx(x,m,e2,s2,g,l);
    fxx2 = Gxx(x,m,e2,s2,g,l);
    fprintf("j = %d\n",j);
    fprintf("m = %d, g = %d, l = %d\n",m,g,l);
    fprintf("e1 = %d, s1 = %d, x = %d\n",e1,s1,x);
    fprintf("G = %d, Gx = %d, Gxx = %d\n",f1,fx1,fxx1);
    fprintf("e2 = %d, s2 = %d, x = %d\n",e2,s2,x);
    fprintf("G = %d, Gx = %d, Gxx = %d\n",f2,fx2,fxx2);
    [x,s,A] = Hfun(m,g,l);
    fprintf("x = %d, s1 = %d, s2 = %d, e1 = %d, e2 = %d\n",x,s(1),s(2),A(1)-m,A(2)-m);

    fprintf("\n");
end

%% test the reduced system formulas
mt1 = @(g)1.5*sqrt(1.5)*(1+g*r3)^(1.5)/(1+g^2)^(0.5);
mt2 = @(g)1.5*sqrt(1.5)*(1-g*r3)^(1.5)/(1+g^2)^(0.5);
st1 = @(g)1.5*sqrt(1.5)*(g-r3)*(1+g*r3)^(0.5)/(1+g^2)^(0.5);
st2 = @(g)1.5*sqrt(1.5)*(g+r3)*(1-g*r3)^(0.5)/(1+g^2)^(0.5);
v1 = @(g)(2/3)/(1+g*r3);
v2 = @(g)(2/3)/(1-g*r3);
Eq = @(m,s,g,x)m^2*(1-x)^2*x + (s-g*m*x)^2*x - 1;

for j = 1 : 10
    if j == 1
        g = 0;
    else
        if j == 2
            g = 1;
        else
            if j == 3
                g = -1;
            else
                g = sqrt(3)*(rand);
            end
        end
    end
    m1 = mt1(g);
    s1 = st1(g);
    x1 = v1(g);
    m2 = mt2(g);
    s2 = st2(g);
    x2 = v2(g);

    fprintf("j = %d, m = %d, s = %d, g = %d, x = %d, Eq = %d\n",j,m1,s1,g,x1,Eq(m1,s1,g,x1));
    fprintf("j = %d, m = %d, s = %d, g = %d, x = %d, Eq = %d\n",j,m2,s2,g,x2,Eq(m2,s2,g,x2));

    l = rand;
    m = 3*rand;

    e1 = ep1(m,g,l);
    e2 = ep2(m,g,l);
    s1 = sig1(m,g,l);
    s2 = sig2(m,g,l);
    x = u2(m,g,l);

    fprintf("mtilde1 = %d, mtilde2 = %d\n",(m + e1)^(3/2)/(l*sqrt(m)), (m + e2)^(3/2)/(l*sqrt(m)));
    fprintf("stilde1 = %d, stilde2 = %d\n",s1*(m + e1)^(1/2)/(l*sqrt(m)), s2*(m + e2)^(1/2)/(l*sqrt(m)));
    fprintf("v1 = %d, v2 = %d\n",x/(m+e1),x/(m+e2));
    fprintf("\n");
end
end
%%

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