clear all
clf

%Setup variables
dt = 60 * 60 * 24;
tmax = 365 * 60 * 60 * 24;
clockmax = ceil(tmax/dt);

G = 6.6738E-11;
M = 5.9736E24;
Ra = 384400;
Rb = 384400;
Pearth = 1.5e11;
DP = 1.2*max(Ra, Rb);

xsave = zeros(1, clockmax);
ysave = zeros(1, clockmax);
X = Ra + Pearth;
Y = 0;
U = 0;
V = sqrt(2*G*M*Rb/(Ra*(Ra + Rb)));

%Setup animation
set(gcf, 'double', 'on');
hsun = plot(Pearth,0,'r*');
hold on;

axis manual;
axis equal;

hearth = plot(X,Y,'bo');
htrail = plot(X,Y);

axis([Pearth - DP,Pearth + DP,0 - DP,0 + DP]);

%hold off;
for clock=1:clockmax
    
    R = sqrt((X - Pearth)^2 + Y^2);
    U = U-dt*G*M*(X - Pearth)/R^3;
    V = V-dt*G*M*Y/R^3;
    X = X+dt*U;
    Y = Y+dt*V;
    xsave(clock) = X;
    ysave(clock) = Y;
    
    %Y = Y + 1000;
    set(hearth,'xdata',X,'ydata',Y);
    set(htrail,'xdata',xsave(1:clock),'ydata',ysave(1:clock));
    disp(X - Pearth);
    
    drawnow;
end
%close all