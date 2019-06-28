phi1=310;
ind=phi1; %index of phi coordinate
r=20;

figure(32153)

cla
X=squeeze(rmesh(ind,:,:));
Y=squeeze(zmesh(ind,:,:));
I=squeeze(I_cyl(ind,:,:));

x0_1 = sin(chi_rad_1m1m1(phi1))*qmax_Fit*100;
y0_1 = cos(chi_rad_1m1m1(phi1))*qmax_Fit*100;
x0_2 = sin(chi_rad_m11m1(phi1))*qmax_Fit*100;
y0_2 = cos(chi_rad_m11m1(phi1))*qmax_Fit*100;

pcolor(X, Y, I); shading('interp')
hold on
circle(x0_1, y0_1, r);
circle(x0_2, y0_2, r);

axis equal

mask1 = (X-x0_1).^2+(Y-y0_1).^2 <r^2;
mask2 = (X-x0_2).^2+(Y-y0_2).^2 <r^2;
mask12 = logical(mask1 .* mask2);

Int1 = sum(I(mask))+sum(I(mask2))-sum(I(mask12))


Int = zeros(1,360);
omega = 0:1:360;
plot_fns = 0;

for phi1=1:361
    ind=phi1; %index of phi coordinate
    X=squeeze(rmesh(ind,:,:));
    Y=squeeze(zmesh(ind,:,:));
    I=squeeze(I_cyl(ind,:,:));

    x0_1 = sin(chi_rad_1m1m1(phi1))*qmax_Fit*100;
    y0_1 = cos(chi_rad_1m1m1(phi1))*qmax_Fit*100;
    x0_2 = sin(chi_rad_m11m1(phi1))*qmax_Fit*100;
    y0_2 = cos(chi_rad_m11m1(phi1))*qmax_Fit*100;
  
    mask1 = (X-x0_1).^2+(Y-y0_1).^2 <r^2;
    mask2 = (X-x0_2).^2+(Y-y0_2).^2 <r^2;
    mask12 = (mask1 & mask2);

    Int(phi1) = sum(I(mask))+sum(I(mask2))-sum(I(mask12));
end

figure(312)
plot(omega, Int)

x = [0, 35.2, 90, 90+35.2, 180-35.2, 180, 180+35.2, 270, 270+35.2, 360-19.5, 360];
for i=1:length(x)
    line([x(i),x(i)],[0.9*min(Int),1.1*max(Int)],'LineStyle', '--', 'color', [0.6,0.6,0.6] , 'LineWidth', 0.5)
end

xlim([0, 360])
ylim([0.9*min(Int),1.1*max(Int)])
xlabel('\omega (deg)')
ylabel('Scattering intensity a.u.')


function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
    hold off
end