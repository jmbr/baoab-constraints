clear all
close all
format long

K = 0.5
temperature = 0.1
nbins = 500

D = @(theta) sqrt(K) .* cos(theta).^2 + sin(theta).^2 ./ sqrt(K);
U = @(theta) cos(theta).^2 .* sin(theta).^2 / K;
%U = @(theta) sin(theta);
f = @(theta) sqrt(D(theta)) .* exp(-U(theta) ./ temperature);

N = 1e7;

z1 = quadgk(f, -pi, pi)

% ts = linspace(-pi, pi, N*10);
% h = ts(2) - ts(1);
% z2 = h*trapz(f(ts))

% zerr = log10(abs(z1-z2))

hst1 = zeros(1, nbins);
% hst2 = zeros(1, nbins);
intervals = linspace(-pi, pi, nbins + 1);

for r = 1:nbins
    % disp(['Step ' num2str(r)]);
    
    a = intervals(r);
    b = intervals(r+1);
    hst1(r) = quadgk(f, a, b) / z1;

    % ts = linspace(a, b, N);
    % dt = ts(2) - ts(1);
    % hst2(r) = dt * trapz(f(ts)) / z2;
end

% err = log10(norm(hst1 - hst2, 'inf'))

% subplot(1, 2, 1);
plot(intervals(2:end), hst1);
data = [intervals(1:end-1)' hst1'];
save('histogram.dat', '-ascii', 'data')
% subplot(1, 2, 2);
% hist(hst2);

% hst3 = [0.06335038037960476 0.07846018979248426 0.08547751219967743 0.0803223327528606 0.06946234526594158 0.06946234526594158 0.08032233275286062 0.08547751219967741 0.07846018979248426 0.06335038037960475 0.04648530341829146 0.03194415120522183 0.02106348636216899 0.01370652625218747 0.009727772371561617 0.009727772371561617 0.01370652625218747 0.02106348636216898 0.03194415120522182 0.04648530341829145];
% %hst3 = [0.0676574 0.0650374 0.0551709 0.0418376 0.031793 0.031793 0.0418376 0.0551709 0.0650374 0.0676574 0.0637844 0.0565438 0.0479972 0.0387983 0.0313801 0.0313801 0.0387983 0.0479972 0.0565438 0.0637844];
% err13 = log10(norm(hst1 - hst3, 'inf'))
% err23 = log10(norm(hst2 - hst3, 'inf'))
