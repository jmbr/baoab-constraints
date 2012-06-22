function data = histogram(K, temperature, nbins)

format long

if nargin < 3
    nbins = 500;
end

D = @(theta) sqrt(K) .* cos(theta).^2 + sin(theta).^2 ./ sqrt(K);
d = 1.0;
a = 2.0;
r0 = 1.0;
r = @(theta) sqrt(cos(theta).^2 / K + sin(theta).^2);
U = @(theta) d * (1 - exp(-a * (r(theta) - r0))).^2;
% U = @(theta) cos(theta).^2 .* sin(theta).^2 / K;
%U = @(theta) sin(theta);
f = @(theta) sqrt(D(theta)) .* exp(-U(theta) ./ temperature);

N = 1e7;

z1 = quadgk(f, -pi, pi);

hst1 = zeros(1, nbins);
intervals = linspace(-pi, pi, nbins + 1);

for r = 1:nbins
    a = intervals(r);
    b = intervals(r+1);
    hst1(r) = quadgk(f, a, b) / z1;
end

figure;
plot(intervals(2:end), hst1);
data = [intervals(1:end-1)' hst1'];
save('histogram.dat', '-ascii', 'data')
