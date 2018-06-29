import MatrixGenerator.*;
% hardcoded directory with naive and recommended implementation
algorithms_dir = fullfile(fileparts(mfilename('fullpath')), 'reference');
addpath(algorithms_dir);
matrices = operand_generator();
naive_ = @() naive(matrices{{:}});
recommended_ = @() recommended(matrices{{:}});
naive_mat = naive_();
recomd_mat = recommended_();
%fprintf('Norm(naive - recomd) %f\\n', norm(naive_mat - recomd_mat));
%fprintf('Naive(0, 0): %f\\n', naive_mat(1, 1));
%fprintf('Naive(0, 0): %f\\n', recomd_mat(1, 1));

benchmarker = Benchmarker();
benchmarker.benchmark('naive_matlab', 10, naive_);
benchmarker.benchmark('recommended_matlab', 10, recommended_);
benchmarker.save('matlab_results_{0}.txt');