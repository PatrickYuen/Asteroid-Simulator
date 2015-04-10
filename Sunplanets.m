clear all
clf

%Setup variables
N = 7; %Set this to 6 and you get the normal orbits of a Sun Earth Moon system
dt = 30*60*24;
tmax = 365 * 60 * 60 * 24;
clockmax = ceil(tmax/dt);

G = 6.6738E-11;
M = [1.989E30, 3.285E23, 4.867E24, 5.9736E24, 7.349E22, 6.39E23, 2.31E21]; %2.31E21
Ra = [0, 5.791E10, 1.082E11, 1.5e11, 384400000, 2.279E11];
Rb = [0, 5.791E10, 1.082E11, 1.5e11, 384400000, 2.279E11];
DAS = 2.2 * 149597870700;
DP = 1.2*DAS;

%positions: Sun, Mercury, Venus, Earth, Moon, Mars

%angles = [0,259.32,63.71,166.59,180,29.09];
angles = [0,90,180,0,0,180];
%angles = [0,0,0,0,0,0];

X = [0, Ra(2) * cosd(angles(2)), Ra(3) * cosd(angles(3)), Ra(4) * cosd(angles(4)), Ra(4) * cosd(angles(4)) + Ra(5), Ra(6) * cosd(angles(6))];
Y = [0, Ra(2) * sind(angles(2)), Ra(3) * sind(angles(3)), Ra(4) * sind(angles(4)), Ra(4) * sind(angles(4)), Ra(6) * sind(angles(6))];
Z = [0, 0, 0, 0, 0, 0];

%speeds
V = [0, 
    -1 * cosd(angles(2)) * sqrt(G*M(1)/Ra(2)),
    cosd(angles(3)) * sqrt(G*M(1)/Ra(3)),
    cosd(angles(4)) * sqrt(G*M(1)/Ra(4)),
    cosd(angles(5)) * sqrt(G*M(4)/Ra(5)),
    cosd(angles(6)) * sqrt(G*M(1)/Ra(6))];

U = [0, 
    -1 * sind(angles(2)) * sqrt(G*M(1)/Ra(2)),
    sind(angles(3)) * sqrt(G*M(1)/Ra(3)),
    sind(angles(4)) * sqrt(G*M(1)/Ra(4)),
    sind(angles(5)) * sqrt(G*M(4)/Ra(5)),
    sind(angles(6)) * sqrt(G*M(1)/Ra(6))];

V(5) = V(4) + V(5);
U(5) = U(4) + U(5);

W = [0, 0, 0, 0, 0, 0];

distsun = zeros(1, clockmax);
distearth = zeros(1, clockmax);
distes = zeros(1, clockmax);

%asteroid
 X(7) = -418851267386.32;
 Y(7) = -643118413961.729;
 Z(7) = 0;
 U(7) = 17361.9304307217;
 V(7) = 11270.2402846395;
 W(7) = 0;

xmercurysave = zeros(1, clockmax);
ymercurysave = zeros(1, clockmax);
xvenussave = zeros(1, clockmax);
yvenussave = zeros(1, clockmax);
xearthsave = zeros(1, clockmax);
yearthsave = zeros(1, clockmax);
xmoonsave = zeros(1, clockmax);
ymoonsave = zeros(1, clockmax);
xmarssave = zeros(1, clockmax);
ymarssave = zeros(1, clockmax);

xastsave = zeros(1, clockmax);
yastsave = zeros(1, clockmax);

%Setup animation
%Subplot 1
subplot(2,2,1);
set(gcf, 'double', 'on');

hold on;

axis manual;
axis equal;

hsun = plot(X(1),Y(1),'r*');
hmercury = plot(X(2),Y(2),'bo');
hvenus = plot(X(3),Y(3),'bo');
hearth = plot(X(4),Y(4),'go');
hmoon = plot(X(5),Y(5),'bo');
hmars = plot(X(6),Y(6),'bo');

hast = plot(X(7),Y(7),'mo');

hasttrail = plot(X(7),Y(7));
hmercurytrail = plot(X(2),Y(2));
hvenustrail = plot(X(3),Y(3));
hmarstrail = plot(X(6),Y(6));
hearthtrail = plot(X(4),Y(4));
hmoontrail = plot(X(5),Y(5));

axis([-DP,DP,-DP,DP]);

%Subplot2
subplot(2,2,2);
hold on;

axis manual;
axis equal;

hmercury2 = plot(Y(2),Z(2),'bo');
hvenus2 = plot(Y(3),Z(3),'bo');
hearth2 = plot(Y(4),Z(4),'go');
hmoon2 = plot(Y(5),Z(5),'bo');
hmars2 = plot(Y(6),Z(6),'bo');

hast2 = plot(Y(7),Z(7),'mo');

axis([-DP,DP,-DP,DP]);

%Subplot3
subplot(2,2,3);

hold on;

axis manual;
axis([0,clockmax,0,1E12]);

hdistsun= plot(0,0);
hdistearth= plot(0,0);
hdistes= plot(0,0);

%Subplot4
s4 = subplot(2,2,4);

hold on;

axis manual;
axis equal;

hearthlocal = plot(X(4),Y(4),'bo');
hmoonlocal = plot(X(5),Y(5),'go');
hastlocal = plot(X(7),Y(7),'mo');
hasttraillocal = plot(X(7),Y(7));
hearthtraillocal = plot(X(4),Y(4));
hmoontraillocal = plot(X(5),Y(5));

axis([X(4) - 1.2*Ra(5),X(4) + 1.2*Ra(5),Y(4)-1.2*Ra(5),Y(4) + 1.2*Ra(5)]);

for clock=1:clockmax   
    if(sqrt((X(4)-X(7))^2+(Y(4)-Y(7))^2+(Z(4)-Z(7))^2) < 2 * Ra(5))
        dt = 2 * 60*24;
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
    xmercurysave(clock) = X(2);
    ymercurysave(clock) = Y(2);
    xvenussave(clock) = X(3);
    yvenussave(clock) = Y(3);
    xmarssave(clock) = X(6);
    ymarssave(clock) = Y(6);
    xearthsave(clock) = X(4);
    yearthsave(clock) = Y(4);
    xmoonsave(clock) = X(5);
    ymoonsave(clock) = Y(5);
    xastsave(clock) = X(7);
    yastsave(clock) = Y(7);
    
    distsun(clock) = sqrt((X(1)-X(7))^2+(Y(1)-Y(7))^2+(Z(1)-Z(7))^2);
    distearth(clock) = sqrt((X(4)-X(7))^2+(Y(4)-Y(7))^2+(Z(4)-Z(7))^2);
    distes(clock) = sqrt((X(1)-X(4))^2+(Y(1)-Y(4))^2+(Z(1)-Z(4))^2);

    set(hsun,'xdata',X(1),'ydata',Y(1));
    
    set(hmercury,'xdata',X(2),'ydata',Y(2));
    set(hvenus,'xdata',X(3),'ydata',Y(3));
    set(hmars,'xdata',X(6),'ydata',Y(6));
    
    set(hearth,'xdata',X(4),'ydata',Y(4));
    set(hearth2,'xdata',Y(4),'ydata',Z(4));
    set(hearthlocal,'xdata',X(4),'ydata',Y(4));
        
    set(hmoon,'xdata',X(5),'ydata',Y(5));
    set(hmoon2,'xdata',Y(5),'ydata',Z(5));
    set(hmoonlocal,'xdata',X(5),'ydata',Y(5));
    
    set(hast,'xdata',X(7),'ydata',Y(7));
    set(hast2,'xdata',Y(7),'ydata',Z(7));
    set(hastlocal,'xdata',X(7),'ydata',Y(7));
    
    set(hdistsun,'xdata',1:clock,'ydata',distsun(1:clock));
    set(hdistearth,'xdata',1:clock,'ydata',distearth(1:clock));
    set(hdistes,'xdata',1:clock,'ydata',distes(1:clock));
    
    set(hmercurytrail,'xdata',xmercurysave(1:clock),'ydata',ymercurysave(1:clock));
    set(hvenustrail,'xdata',xvenussave(1:clock),'ydata',yvenussave(1:clock));
    set(hearthtrail,'xdata',xearthsave(1:clock),'ydata',yearthsave(1:clock));
    set(hmoontrail,'xdata',xmoonsave(1:clock),'ydata',ymoonsave(1:clock));
    set(hmarstrail,'xdata',xmarssave(1:clock),'ydata',ymarssave(1:clock));
    set(hasttrail,'xdata',xastsave(1:clock),'ydata',yastsave(1:clock));
    
    set(hearthtraillocal,'xdata',xearthsave(1:clock),'ydata',yearthsave(1:clock));
    set(hmoontraillocal,'xdata',xmoonsave(1:clock),'ydata',ymoonsave(1:clock));
    set(hasttraillocal,'xdata',xastsave(1:clock),'ydata',yastsave(1:clock));
    
    axis(s4,[X(4) - 2*Ra(5),X(4) + 2*Ra(5),Y(4)-2*Ra(5),Y(4) + 2*Ra(5)]);
    drawnow;
end
%close all