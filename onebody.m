clear all
clf
G = 6.67384E-11;
M = 5.9736E24;
Ra = 384400;
Rb = 384400;
DP = 1.2*max(Ra,Rb);
tmax = 365.25*60*60*24;
dt = 1;
clockmax = ceil(tmax/dt);

X = Ra;
Y = 0;
U = 0;
V = sqrt(2*G*M*Rb/(Ra*(Ra+Rb)));

set(gcf, 'double', 'on');
hold on
hearth = plot(0,0, 'go');
hmoon = plot(X,Y, 'bo');
axis equal; %Scale between x and y directions
axis manual; %To force the axis to be the same 
axis([-DP,DP,-DP,DP]); % Set xlim and ylim

for clock = 1:clockmax
	R = sqrt(X^2 + Y^2) % Update the distance between earth and sun
	U = U - dt*G*M*X/R^3 % Update the x component of velocity
	V = V - dt*G*M*Y/R^3 % Update the y component of velocity
    X = X + dt*U % Update the x coordinate of Earth
    Y = Y + dt*V % Update the y coordinate of Earth
    set(hmoon,'xdata',X,'ydata',Y);
    drawnow;
end
%movie2avi(Mov,'GravityAnimation.avi','compression','None');
%mpgwrite(Mov,jet,'movie.mpg')