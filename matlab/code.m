%% Ylä ja alasummat

%f=@(x) x.^3-x.^2-x+4;
f=@(x) exp(-x.^2);
a=-2;
b=2;
n=21;
x=linspace(a,b,300);

p=linspace(a,b,n);


hold on;

for i=1:n-1    
    if(p(i+1) <= 0)
        h_r = area([p(i) p(i+1)],[f(p(i+1)) f(p(i+1))]); % right
        h_l = area([p(i) p(i+1)],[f(p(i)) f(p(i))]); % left        
        set(h_l,'FaceColor',[0.9 0.9 1]);
        set(h_l,'EdgeColor',[0 0 1]);
        set(h_r,'FaceColor',[0.8 0.8 1]);
        set(h_r,'EdgeColor',[0 0 1]);
    else
        h_l = area([p(i) p(i+1)],[f(p(i)) f(p(i))]); % left
        h_r = area([p(i) p(i+1)],[f(p(i+1)) f(p(i+1))]); % right                
        set(h_l,'FaceColor',[0.8 0.8 1]);
        set(h_l,'EdgeColor',[0 0 1]);
        set(h_r,'FaceColor',[0.9 0.9 1]);
        set(h_r,'EdgeColor',[0 0 1]);
    end    
end;
h = plot(x,f(x),'b','LineWidth',2);

hold off;
legend([h h_l h_r],'e^{-x^2}','upper sum','lower sum');

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