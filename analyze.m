clear all
close all

% function slope = analyze(directory)

% if nargin == 0
directory = pwd;
% end

files = dir(fullfile(directory, 'result-*.dat'));
num_dts = length(files);

dts = zeros(num_dts, 1);
results = {};

total_time = +Inf;

for k = 1:num_dts
    name = files(k).name;
    out = textscan(name, 'result-0.5-0.1-%f.dat');
    % out = textscan(name, 'result-dt-%f.dat');
    dts(k) = out{1};

    results{k} = load(name);

    r = results{k};
    total_time = min(total_time, r(end, 1));
end

N = 1000;
errors = {};
ts = linspace(0, total_time, N);
for k = 1:num_dts
    r = results{k};
    errors{k} = interp1(r(:, 1), r(:, 2), ts);
end

slopes = zeros(N, 1);
residuals = zeros(N, 1);
for n = 1:N
    t = ts(n);

    x = zeros(1, num_dts);
    y = zeros(1, num_dts);
    for k = 1:num_dts
        x(k) = log10(dts(k));
        err = errors{k};
        y(k) = log10(err(n));
    end

    p = polyfit(x, y, 1);
    slopes(n) = p(1);
    residuals(n) = mean(abs(err(n) - p));

    if mod(n, 500) == 0
        disp(['Progress: ' num2str(round(t/total_time*100))]);
    end
end

subplot(3, 1, 1);
plot(ts(1:n), slopes(1:n), '-');
%  axis([0 total_time]);
%set(gca, 'YTick', -6:1:6);
grid on;
xlabel('Time');
ylabel('Slope of least-squares line');

subplot(3, 1, 2);
plot(x, y, 'bo', ...
     x, polyval(p, x), 'k--');
% axis([log10(dts(1)) log10(dts(end)) -5 -2]);
legend('Observations', ...
       sprintf('Least-squares line\n(Slope %f)', slopes(n)), ...
       'Location', 'NorthWest');
xlabel('Time step length (Log. base 10)');
ylabel('Error (Log. base 10)');
drawnow;
