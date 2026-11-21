% Writes the expected outputs of the MATLAB implementation of TriGer for the
% fixture inputs in tests/testthat/fixtures/, which make_fixtures.R creates.
% The output files and their format are those described in make_fixtures.R.
%
% Run from the R package root, with the MATLAB functions of TriGer
% (https://github.com/jy2941/TriGer, folder matlab/methods) on the path:
%     addpath('<path to TriGer>/matlab/methods')
%     run(fullfile('inst', 'scripts', 'make_fixtures.m'))
%
% In MATLAB, max_round = 5 returns at most 4 blocks, which corresponds to
% max_blocks = 4 in the R and Python packages. The block densities of
% greedy_peeling_XY_all, which the MATLAB function does not return, are
% computed from the definition sum(XY(x, y)) / (|x| * |y|)^(lambda / 2).

out = fullfile('tests', 'testthat', 'fixtures');

XY = readmatrix(fullfile(out, 'xy.csv'));
blocks = greedy_peeling_XY_all(XY, 1.5, 5);
blocks = blocks(~cellfun(@isempty, blocks(:, 1)), :);
dens = zeros(size(blocks, 1), 1);
for k = 1:size(blocks, 1)
    x = blocks{k, 1};
    y = blocks{k, 2};
    dens(k) = sum(sum(XY(x, y))) / (numel(x) * numel(y))^(1.5 / 2);
end
write_blocks(blocks, {'X', 'Y'}, fullfile(out, 'xy_blocks.csv'));
write_density(dens, fullfile(out, 'xy_density.csv'));

S = readmatrix(fullfile(out, 'x.csv'));
comms = greedy_peeling_X_all(S, 1.5);
comms = comms(~cellfun(@isempty, comms(:, 1)), :);
write_blocks(comms(:, 1), {'X'}, fullfile(out, 'x_blocks.csv'));
write_density(double(cell2mat(comms(:, 2))), fullfile(out, 'x_density.csv'));

P = readmatrix(fullfile(out, 'three_step_P.csv'));
R = readmatrix(fullfile(out, 'three_step_R.csv'));
modules = greedy_3step(P, R, 30, 40, 1.1, 1.5, 1.5, [5 10], [5 15]);
write_blocks(modules, {'X', 'Y'}, fullfile(out, 'three_step_blocks.csv'));

function write_blocks(blocks, sides, file)
    fid = fopen(file, 'w');
    fprintf(fid, 'block,side,index\n');
    for k = 1:size(blocks, 1)
        for s = 1:numel(sides)
            idx = blocks{k, s};
            for i = 1:numel(idx)
                fprintf(fid, '%d,%s,%d\n', k, sides{s}, idx(i));
            end
        end
    end
    fclose(fid);
end

function write_density(dens, file)
    fid = fopen(file, 'w');
    fprintf(fid, 'block,density\n');
    for k = 1:numel(dens)
        fprintf(fid, '%d,%.17g\n', k, dens(k));
    end
    fclose(fid);
end
