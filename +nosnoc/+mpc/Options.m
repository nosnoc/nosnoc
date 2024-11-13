classdef Options < handle
    properties
        fullmpcc_fast_sigma_0(1,1) double {mustBeNonnegative} = 1e-9;
        fullmpcc_fast_N_homotopy(1,1) double {mustBeNonnegative, mustBeInteger} = 2;
        fullmpcc_do_shift_initialization(1,1) logical = true;
        fullmpcc_progressive_relaxation(1,1) double = 0;
    end
end
