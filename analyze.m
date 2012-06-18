function slope = analyze(directory)

if nargin == 0
    directory = pwd;
end

files = dir(fullfile(directory, 'result-*.dat'));
num_files = length(files);

time_step_sizes = zeros(num_files, 1);
results = zeros(num_files, 2);

for k = 1:num_files
    name = files(k).name;
    out = textscan(name, 'result-dt-%f.dat');
    time_step_sizes(k) = out{1};

    data = load(name);

    results(k, :) = data(end, :);
end

[time_step_sizes indices] = sort(time_step_sizes);
results = results(indices, :);

x = log10(time_step_sizes);
y = log10(results(:, 2));

p = polyfit(x, y, 1);
slope = p(1);

plot(x, y, 'bx-', ...
     x, polyval(p, x), 'k--');
legend('Observations', ...
       sprintf('Least-squares line (slope %f)', slope), ...
       'Location', 'SouthEastOutside');
xlabel('Time step length (Log. base 10)');
ylabel('Error (Log. base 10)');
axis tight;
