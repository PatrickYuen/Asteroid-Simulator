clear all
clf

%Setup variables
N = 4; %Set this to 3 and you get the normal orbits of a Sun Earth Moon system
dt = 30*60*24;
tmax = 365 * 60 * 60 * 24;
clockmax = ceil(tmax/dt);

G = 6.6738E-11;
M = [1.989E30, 5.9736E24, 7.349E22, 2.31E21]; %2.31E21
Ra = zeros(1,N);
Rb = zeros(1,N);
Ra(1) = 0;
Rb(1) = 0;
Ra(2) = 1.5e11;
Rb(2) = 1.5e11;
Ra(3) = 384400000;
Rb(3) = 384400000;
DP = 2*max(Rb(2) + Rb(3), Ra(2) + Ra(3));
DAS = 2.2 * 149597870700;

X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);
U = zeros(1,N);
V = zeros(1,N);
W = zeros(1,N);

%positions: Sun, Moon, Earth
X(1) = 0;
Y(1) = 0;
Z(1) = 0;

X(2) = Ra(2);
Y(2) = 0;
Z(2) = 0;

X(3) = Ra(2) + Ra(3);
Y(3) = 0;
Z(3) = 0;

%speeds: Sun Moon Earth
U(1) = 0;
U(2) = 0;
U(3) = 0;

V(1) = 0;
V(2) = sqrt(2*G*M(1)*Rb(2)/(Ra(2)*(Ra(2) + Rb(2))));
DMS = Rb(3) + Rb(2);
V(3) =sqrt(2*G*M(2)*Rb(3)/(Ra(3)*(Ra(3) + Rb(3)))) + V(2);

W(1) = 0;
W(2) = 0;
W(3) = 0;

%asteroid
X(4) = DAS;
Y(4) = -643118413961.729;
Z(4) = 0;
U(4) = 17361.9304307217;
V(4) = 11270.2402846395;
W(4) = 0;
%{
X(4) = -418961267386.32;
Y(4) = -643118413961.729;
Z(4) = 0;
U(4) = 17361.9304307217;
V(4) = 11270.2402846395;
W(4) = 0;
%}

xsave = zeros(1, clockmax);
ysave = zeros(1, clockmax);

%Setup animation
%Subplot 1
subplot(2,2,1);
set(gcf, 'double', 'on');

hold on;

axis manual;
axis equal;

hsun = plot(X(1),X(1),'r*');
hearth = plot(X(2),Y(2),'bo');
hmoon = plot(X(3),Y(3),'go');
hast = plot(X(4),Y(4),'mo');
htrail = plot(X(4),Y(4));

axis([-DP,DP,-DP,DP]);

%Subplot2
subplot(2,2,2);
hold on;

axis manual;
axis equal;

hearth2 = plot(Y(2),Z(2),'bo');
hmoon2 = plot(Y(3),Z(3),'go');
hast2 = plot(Y(4),Z(4),'mo');

axis([-DP,DP,-DP,DP]);

%Subplot3
subplot(2,2,3);

hold on;

axis manual;
axis equal;

hsun3 = plot3(X(1),Y(1),Z(1),'r*');
hearth3 = plot3(X(2),Y(2),Z(2),'bo');
hmoon3 = plot3(X(3),Y(3),Z(3),'go');
hast3 = plot3(X(4),Y(4),Z(4),'mo');

%axis([-DP,DP,-DP,DP,-DP,DP]);
AX = 1.2 * DAS;
axis([-AX,AX,-AX,AX,-AX,AX]);

%Subplot4
s4 = subplot(2,2,4);

hold on;

axis manual;
axis equal;

hearthlocal = plot(X(2),Y(2),'bo');
hmoonlocal = plot(X(3),Y(3),'go');
hastlocal = plot(X(4),Y(4),'mo');
axis([X(2) - 1.2*Ra(3),X(2) + 1.2*Ra(3),Y(2)-1.2*Ra(3),Y(2) + 1.2*Ra(3)]);

for clock=1:clockmax   
    if(sqrt((X(2)-X(4))^2+(Y(2)-Y(4))^2+(Z(2)-Z(4))^2) < 2 * Ra(3))
        dt = 2*60*24;
    end
    for i = 1:N
        for j = 1:N
            if(i ~= j)
                DX = X(i) - X(j);
                DY = Y(i) - Y(j);
                DZ = Z(i) - Z(j);
                R = sqrt(DX^2 + DY^2 + DZ^2);
                U(i) = U(i)-dt*G*M(j)*DX/R^3;
                V(i) = V(i)-dt*G*M(j)*DY/R^3;
                W(i) = W(i)-dt*G*M(j)*DZ/R^3;
            end
        end
    end
    
    for i = 1:N
        X(i) = X(i)+dt*U(i);
        Y(i) = Y(i)+dt*V(i);
        Z(i) = Z(i)+dt*W(i);
    end
    
    xsave(clock) = X(4);
    ysave(clock) = Y(4);

    set(hsun,'xdata',X(1),'ydata',Y(1));
    set(hsun3,'xdata',X(1),'ydata',Y(1),'zdata',Z(1));
    
    set(hearth,'xdata',X(2),'ydata',Y(2));
    set(hearth2,'xdata',Y(2),'ydata',Z(2));
    set(hearth3,'xdata',X(2),'ydata',Y(2),'zdata',Z(2));
    set(hearthlocal,'xdata',X(2),'ydata',Y(2));
        
    set(hmoon,'xdata',X(3),'ydata',Y(3));
    set(hmoon2,'xdata',Y(3),'ydata',Z(3));
    set(hmoon3,'xdata',X(3),'ydata',Y(3),'zdata',Z(3));
    set(hmoonlocal,'xdata',X(3),'ydata',Y(3));
    
    set(hast,'xdata',X(4),'ydata',Y(4));
    set(hast2,'xdata',Y(4),'ydata',Z(4));
    set(hast3,'xdata',X(4),'ydata',Y(4),'zdata',Z(4));
    set(hastlocal,'xdata',X(4),'ydata',Y(4));
    set(htrail,'xdata',xsave(1:clock),'ydata',ysave(1:clock));
    
    axis(s4,[X(2) - 2*Ra(3),X(2) + 2*Ra(3),Y(2)-2*Ra(3),Y(2) + 2*Ra(3)]);
    drawnow;
end
%close all