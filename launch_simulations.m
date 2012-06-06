clear all
close all

common;                                 % Load global variables.

% Make sure we're using the latest revision of the code.
status = system('make');
if status ~= 0
    error('Unable to compile source code');
end

%rng(now, 'twister');

script = 'simulation-script.sh';
f = fopen(script, 'w');

for d = 1:num_dts
    dt = dts(d);
    nsteps = ceil(total_time / dt);

    seed = floor(rand * 1e7);
    outfile = sprintf('result-%g-%g.dat', temperature, dt);
    cmd = sprintf('./baoab %g %g %g %g %u > %s', ...
                  temperature, friction, dt, nsteps, seed, outfile);

    fprintf(f, '%s\n', cmd);
end

fclose(f);

tic; system(['cat ' script '  | parallel']); toc;

analyze_simulations;

system('xmessage -center "Simulation finished."');
