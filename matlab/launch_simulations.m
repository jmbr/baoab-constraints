clear all
close all

config;                                 % Load configuration file.

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

    outfile = sprintf('result-%02.04g-%02.04g-%02.04g.dat', K, temperature, dt);
    logfile = sprintf('log-%02.04g-%02.04g-%02.04g.dat', K, temperature, dt);
    cmd = sprintf(['./simul --K %02.04g --temperature %02.04g --friction %02.04g' ...
                   ' --dt %02.04g --steps %02.04g --seed %u --bins %u > %s 2> %s'], ...
                  K, temperature, friction, dt, nsteps, ...
                  seed, bins, outfile, logfile);

    fprintf(f, '%s\n', cmd);
end

fclose(f);

tic; system(['cat ' script '  | parallel']); toc;

system('xmessage -center "Simulation finished."');
