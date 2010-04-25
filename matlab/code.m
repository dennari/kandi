%% Numeerisia integroimismenetelmiä

f=@(x) exp(-x.^2);
a=0;
b=1;
sum_s = quad(f,a,b); % oikea arvo (Simpsonin säännöllä)
multiplier = 2; % used in symmetrical cases 



%% Ylä ja alasummat

n=20; % intervallien lukumäärä
x=linspace(a,b,300);
p=linspace(a,b,n+1);
w = p(2)-p(1);

sum_u = 0;
sum_l = 0;

hold on;

for i=1:n

    if(f(p(i)) < f(p(i+1)))
       l = f(p(i));
       u = f(p(i+1));
    else
       l = f(p(i+1));
       u = f(p(i));
    end
    
    h_u = fill([p(i) p(i) p(i+1) p(i+1)]',[0 u u 0]',[1 0 0],...
            'FaceAlpha',0.25,'EdgeColor',[0 0 1],'EdgeAlpha',0.0); % upper
    h_l = fill([p(i) p(i) p(i+1) p(i+1)]',[0 l l 0]',[0 0 1],...
            'FaceAlpha',0.3,'EdgeColor',[0 0 1],'EdgeAlpha',0.1); % lower
    
    sum_u = sum_u+w*f(p(i));
    sum_l = sum_l+w*f(p(i+1));
    %get(h_l);

  
end;
h = plot(x,f(x),'k','LineWidth',2);

hold off;
legend([h h_u h_l],'e^{-x^2}','yläsumma','alasumma');


%% print figure in pdf
figLength = 6;
xyRatio = 16/10;
figHeight = (1/xyRatio)*figLength;

set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]*1.7);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figLength figHeight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 -0.03 figLength figHeight]);
print(gcf, '-dpdf', '../fig_lusums.pdf');

%% Puolisuunnikassääntö

n=2; % intervallien lukumäärä
p=linspace(a,b,n+1);
w = p(2)-p(1);
sum_trapezoidal = 0;
figure;
hold on;

for i=1:n
    h_t = area([p(i) p(i+1)],[f(p(i)) f(p(i+1))]);
    sum_trapezoidal = 0.5*w*(f(p(i))+f(p(i+1)));  
    set(h_t,'FaceColor',[0.9 0.9 1]);
    set(h_t,'EdgeColor',[0 0 1]);
end

h = plot(x,f(x),'b','LineWidth',2);

hold off;
legend([h h_t],'e^{-x^2}','puolisuunnikas');

%% Simpsonin sääntö

n=1; % intervallien lukumäärä
p=linspace(a,b,2*n+1);
w = p(2)-p(1);
sum_simpson = 0;

for i=1:n
    sum_simpson = sum_simpson + (w/3)*(f(p(i))+4*f(p(i+1))+f(p(i+2)));  
end

disp(sprintf('simpson: %.7f',multiplier*sum_simpson));
disp(sprintf('error: %.7f',abs(multiplier*(sum_simpson-sum_s))));

%% print figure in pdf
figLength = 6;
xyRatio = 16/10;
figHeight = (1/xyRatio)*figLength;
axis([-2 2 0 1.02]);
set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]*1.7);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figLength figHeight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0.01 figLength figHeight]);
print(gcf, '-dpdf', '../fig_trapezoid.pdf');

%% Kolmipisteinen Gaussin kvadratuuri välillä [a,b] ja painofunktiolla w=1

Q3=@(f) 0.5*(b-a)*((5/9)*f(0.5*(b-a)*-sqrt(0.6)+0.5*(b+a))...
        +(8/9)*f(0.5*(b+a))...
        +(5/9)*f(0.5*(b-a)*sqrt(0.6)+0.5*(b+a)));
disp(sprintf('quadrature: %.7f',multiplier*Q3(f)));
disp(sprintf('error: %.7f',abs(multiplier*(Q3(f)-sum_s))));


%% Koeasetelma

a = -1;
b = 1;
step = 0.04;
[X,Y] = meshgrid(a*2:step:b*2);
% integroimisalue
X1 = [-1 -1 1 -1; 1 -1 1 1; 1 -1 1 1; -1 -1 1 -1];
Y1 = [-1 -1 -1 1; -1 1 1 1; -1 1 1 1; -1 -1 -1 1];
Z1 =[0 0 0 0; 0 0 0 0; 1 1 1 1; 1 1 1 1];
results = [];
results2 = [];
%figure;


%% funktio 1 - gaussinen
subplot(1,2,1);

c = 0.7;
f2=@(x,y,c) exp(-1*(2*(1-c^2))^-1*(x.^2+y.^2-2*c*x.*y));

[cdf,cErr] = mvncdf([a a],[b b],[],[1 c; c 1],statset('TolFun',1e-14));
correct = cdf*2*pi*sqrt(1-c^2);
Z = f2(X,Y,c);
surf(X,Y,Z,'FaceAlpha',0.8,'EdgeAlpha',0.3);
%axis([-2 2 -2 2 -1 1]);
%alpha(0.7);
hold on;
fill3(X1,Y1,max(max(Z))*Z1,'r','FaceAlpha',0.4,'EdgeAlpha',0.1);
hold off;
axis square;
title('$$F_1=e^{-\frac{1}{1.02}\left(x^2+y^2-1.4xy\right)}$$','interpreter','latex','FontSize',14);
%t=text(0.5,0.5,'$$e^{-\frac{1}{1.02}\left(x^2+y^2-1.4xy\right)}$$','interpreter','latex');
%t=text(0.5,0.5,'e^{-1/1.02\left(x^2+y^2-1.4xy\right)}','interpreter','tex');
%get(t);
%% funktio 2 - itseisarvo

subplot(1,2,2);
f3=@(x,y) abs(x.*y).^0.5;%+abs(x)-0.5*y;
correct2 = dblquad(f3,a,b,a,b,1e-8);
Z = f3(X,Y);
h = surf(X,Y,Z,'FaceAlpha',0.8,'EdgeAlpha',0.3,'FaceColor','interp');
shading faceted;
%get(h);
hold on;
h = fill3(X1,Y1,max(max(Z))*Z1,'r','FaceAlpha',0.4,'EdgeAlpha',0.1);

% line = [a a:step:b b];
% const = ones(1,length(line));
% X1 = [line; const; -const; line];
% Y1 = [const; line; line; -const];
% fV = f2(X1,Y1,c);
% fV(:,1) = zeros(4,1);
% fV(:,end)=zeros(4,1);
% fill3(X1',Y1',1.4*fV','r','FaceAlpha',0.8,'EdgeAlpha',0.3);

hold off;
axis square;
title('$$F_2=\sqrt{\vert xy \vert}$$','interpreter','latex','FontSize',14);
%get(h);

%% print figure in pdf

%previewfig(gcf,'Color','rgb');
%exportfig(gcf,'../fig_3dfunc.png','Format','png','Color','rgb','Renderer','opengl');
print(gcf, '-djpeg100','-r350','../fig_3dfunc.jpg');
%%
figLength = 6;
xyRatio = 16/10;
figHeight = (1/xyRatio)*figLength;
%axis([-2 2 0 1.02]);
%set(gca, 'Position', get(gca, 'OuterPosition') - ...
%    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]*1.7);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figLength figHeight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0.01 figLength figHeight]);
print(gcf, '-dpdf', '../fig_3dfunc.pdf');

%% 9-pisteinen tulosääntö

points = 3;
prodRule2 = gauss2D(@(x,y)f2(x,y,c),points-1);
err = abs(correct-prodRule2);
relErr = (err/correct)*100;
fprintf('gauss2D-%d: %.14f, err: %.14f\n',points^2,prodRule2,err);
results(1,1:2)=[prodRule2 relErr];


prodRule2 = gauss2D(f3,points-1);
err = abs(correct2-prodRule2);
relErr = (err/correct)*100;
fprintf('gauss2D-%d: %.14f, err: %.14f\n',points^2,prodRule2,err);
results(1,3:4)=[prodRule2 relErr];


%% radon 7-pisteinen

radRule = radonSquare(@(x,y)f2(x,y,c));
err = abs(correct-radRule);
relErr = (err/correct)*100;
fprintf('radon: %.5f, err: %.5f\n',radRule,abs(correct-radRule));
results(2,1:2)=[radRule relErr];

radRule = radonSquare(f3);
err = abs(correct2-radRule);
relErr = (err/correct)*100;
fprintf('radon: %.5f, err: %.5f\n',radRule,abs(correct2-radRule));
results(2,3:4)=[radRule relErr];

%% minimaalinen d=19,N=68,S_2
minim = minimal19(@(x,y)f2(x,y,c));
err = abs(correct-minim);
relErr = (err/correct)*100;
fprintf('minimal19: %.12f, err: %.12f\n',minim,abs(correct-minim));
results(3,1:2)=[minim relErr];

minim = minimal19(f3);
err = abs(correct2-minim);
relErr = (err/correct)*100;
fprintf('minimal19: %.12f, err: %.12f\n',minim,abs(correct2-minim));
results(3,3:4)=[minim relErr];

%% exporttaa tulokset

matrix2latex(results,'~/Documents/School/Kandi/matlab/results1.tex','rowLabels',{'$Q_{3^2}$' '$Q_{7}$' '$Q_{68}$'},...
    'columnLabels',{'$Q[F_1]$' '$E[F_1]$' '$Q[F_2]$' '$E[F_2]$'},'alignment','c','format',{'$%.5f$' '$%.1f\\%%$' '$%.5f$' '$%.1f\\%%$'});







