function hst = histogram(temperature, nbins)

format long

if nargin < 3
    nbins = 500;
end

d = 1;
a = 2;
r0 = sqrt(2.0);
r = @(theta) sqrt((cos(theta) - 1).^2 + sin(theta).^2);
U = @(theta) d * (1 - exp(-a * (r(theta) - r0))).^2;
f = @(theta) exp(-U(r(theta)) ./ temperature);

partition_function = quadgk(f, -pi, pi);

hst = zeros(1, nbins);
intervals = linspace(-pi, pi, nbins + 1);

for r = 1:nbins
    a = intervals(r);
    b = intervals(r+1);
    hst(r) = quadgk(f, a, b) / partition_function;
end

figure;
plot(intervals(2:end), hst);
% data = [intervals(1:end-1)' hst'];
% save('histogram.dat', '-ascii', 'data')
