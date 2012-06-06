clear all
close all

temperature = 1;
beta = 1/temperature;

a = 5;

D = @(theta) sqrt(a) * cos(theta).^2 + 1 / sqrt(a) * sin(theta).^2;
distribution = @(theta) sqrt(D(theta)) .* exp(-beta * sin(theta));
partition_function = quadgk(distribution, -pi, pi);

for nbins = 10:20
    intervals = linspace(-pi, pi, nbins+1);
    for n = 2:nbins+1
        histogram(n-1) = quadgk(distribution, intervals(n-1), intervals(n)) / partition_function;
    end

    x = intervals(2:end);
    plot(x, histogram, 'bx-');
    xlabel('Angle');
    ylabel('Frequency');
    grid on;
    drawnow;
end
