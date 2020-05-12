function plot_circle(h,x,y,r)
th = 0:pi/50:2*pi;
for i = 0.001:0.01:r;
    xunit = i * cos(th) + x;
    yunit = i * sin(th) + y;
    plot(xunit, yunit,'k');
end
end