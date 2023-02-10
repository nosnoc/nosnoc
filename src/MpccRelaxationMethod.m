classdef MpccRelaxationMethod
    enumeration
        equality   % [0,0]
        inequality % [-inf, 0]
        two_sided  % Only for elastic modes
    end
end
