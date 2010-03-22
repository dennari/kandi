%% Numeerisia integroimismenetelmiä

f=@(x) exp(-x.^2);
a=0;
b=2;
sum_s = quad(f,a,b); % oikea arvo (Simpsonin säännöllä)




%% Ylä ja alasummat

n=20; % intervallien lukumäärä
x=linspace(a,b,300);
p=linspace(a,b,n+1);
w = p(2)-p(1);

sum_u = 0;
sum_l = 0;

hold on;

for i=1:n
    if(p(i+1) <= 0)
        h_r = area([p(i) p(i+1)],[f(p(i+1)) f(p(i+1))]); % right
        h_l = area([p(i) p(i+1)],[f(p(i)) f(p(i))]); % left        
        sum_u = sum_u+w*f(p(i+1));
        sum_l = sum_l+w*f(p(i));
        set(h_l,'FaceColor',[0.9 0.9 1]);
        set(h_l,'EdgeColor',[0 0 1]);
        set(h_r,'FaceColor',[0.8 0.8 1]);
        set(h_r,'EdgeColor',[0 0 1]);
    else
        h_l = area([p(i) p(i+1)],[f(p(i)) f(p(i))]); % left
        h_r = area([p(i) p(i+1)],[f(p(i+1)) f(p(i+1))]); % right                
        sum_u = sum_u+w*f(p(i));
        sum_l = sum_l+w*f(p(i+1));
        set(h_l,'FaceColor',[0.8 0.8 1]);
        set(h_l,'EdgeColor',[0 0 1]);
        set(h_r,'FaceColor',[0.9 0.9 1]);
        set(h_r,'EdgeColor',[0 0 1]);
    end    
end;
h = plot(x,f(x),'b','LineWidth',2);

hold off;
legend([h h_l h_r],'e^{-x^2}','yläsumma','alasumma');


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
disp(sprintf('quadrature: %.7f',Q3(f)));
disp(sprintf('correct: %.7f',sqrt(pi)*erf(b)));
disp(sprintf('error: %.7f',abs(sqrt(pi)*erf(b)-Q3(f))));

