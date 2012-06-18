function f = analyze_simulations(K, first, last)

config;                                 % Load configuration file.

dts

for d = 1:num_dts
    dt = dts(d);
    nsteps = ceil(total_time / dt);
    outfile = sprintf('result-%02.04g-%02.04g-%02.04g.dat', K, temperature, dt);
    value = load(outfile);
    results(d) = value(end);
end

x = log10(dts);
y = log10(results);

xx = x(first:last);
yy = y(first:last);

p = polyfit(xx, yy, 1);
slope = p(1)

f = figure;
plot(x, y, 'ko', ...
     xx, polyval(p, xx), 'k--');
legend('Observations', 'Least-squares line', 'Location', 'SouthEastOutside');
xlabel('Time step length (Log. base 10)');
ylabel('Error (Log. base 10)');
title(['Slope: ', num2str(slope)]);
%axis tight;
%axis([-1.05 -1 -4.8 -4.4]);
