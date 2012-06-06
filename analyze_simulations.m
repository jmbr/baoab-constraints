clear all
close all

common;                                 % Load global variables.

for d = 1:num_dts
    dt = dts(d);
    nsteps = ceil(total_time / dt);
    outfile = sprintf('result-%g-%g.dat', temperature, dt);
    value = load(outfile);
    results(d) = value(end);
end

x = log10(dts);
y = log10(results);

%n = num_dts-1;
n = 6;
xx = x(end-n:end);
yy = y(end-n:end);

p = polyfit(xx, yy, 1);
slope = p(1);

plot(x, y, 'ko', ...
     xx, polyval(p, xx), 'k--');
legend('Observations', 'Least-squares line', 'Location', 'SouthEast');
xlabel('Time step length (Log. base 10)');
ylabel('Error (Log. base 10)');
title(['Slope: ', num2str(slope)]);
%axis([-1.05 -1 -4.8 -4.4]);
