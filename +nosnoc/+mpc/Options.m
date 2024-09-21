classdef Options < handle
    properties
        fullmpcc_fast_sigma_0(1,1) double {mustBeNonnegative} = 1e-9;
        fullmpcc_do_shift_initialization(1,1) logical = true;
    end
end
