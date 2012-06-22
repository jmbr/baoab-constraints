function slope = analyze(directory)

if nargin == 0
    directory = pwd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read input files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = dir(fullfile(directory, 'result-dt-*.dat'));
num_dts = length(files);

dts = zeros(num_dts, 1);
results = {};

total_time = +Inf;

for k = 1:num_dts
    name = files(k).name;
    out = textscan(name, 'result-dt-%f.dat');
    dts(k) = out{1};

    results{k} = load(name);

    r = results{k};
    total_time = min(total_time, r(end, 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain slopes of linear regressions over time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000;
ts = linspace(0, total_time, N);

errors = {};
for k = 1:num_dts
    r = results{k};
    errors{k} = interp1(r(:, 1), r(:, 2), ts);
end

slopes = zeros(N, 1);
% residuals = zeros(N, 1);
h = waitbar(0, 'Starting computation...');
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
    % residuals(n) = max(abs(err(n) - p));

    if mod(n, 100) == 0
        val = t/total_time;
        waitbar(val, h, sprintf('%u%% completed', round(val*100)));
    end
end
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts_slopes = [transpose(ts) slopes];
save('time-vs-slopes.dat', '-ascii', 'ts_slopes');

last_regression = transpose([x; y; polyval(p, x)]);
save('last-regression-slope.dat', '-ascii', 'last_regression');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

subplot(1, 2, 1);
plot(ts, slopes, '-');
grid on;
xlabel('Time');
ylabel('Slope of least-squares line');

subplot(1, 2, 2);
plot(x, y, 'bo', x, polyval(p, x), 'k--');
% axis([log10(dts(1)) log10(dts(end)) -5 -2]);
legend('Observations', ...
       sprintf('Least-squares line\n(Slope %f)', slopes(n)), ...
       'Location', 'NorthWest');
xlabel('Time step length (Log. base 10)');
ylabel('Error (Log. base 10)');

% subplot(3, 1, 2);
% plot(ts, residuals, '-');
% grid on;
% xlabel('Time');
% ylabel('Residuals of least-square regression.');

slope = slopes(end);
