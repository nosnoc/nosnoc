classdef Base < matlab.mixin.Scalar & handle & matlab.mixin.CustomDisplay
% Base class for Dynamic Complementarity Systems of the (most generic) form.
%
% .. math::
%     :nowrap:
%
%     \begin{align*}
%        \dot{x}&= f(x,z) \\
%        0 &\le h(x,z) \perp z \ge 0
%     \end{align*}
    properties
        model

        f_x_fun % casadi.Function: Right hand side of the DCS.
        f_q_fun % casadi.Function: Lagrange cost term function.
        g_z_fun % casadi.Function: User algebraic constraints function.
        g_alg_fun % casadi.Function: Constraints for algorithmically defined algebraics. 
        g_path_fun % casadi.Function: Function for non-box path constraints. Only defined in an OCP.
        G_path_fun % casadi.Function: Function for one half of path complementarity. Only defined in an OCP.
        H_path_fun % casadi.Function: Function for one half of path complementarity. Only defined in an OCP.
        g_terminal_fun % casadi.Function: Function for terminal constraints. Only defined in an OCP.
        f_q_T_fun % casadi.Function: Function for Mayer cost term. Only defined in an OCP.
        f_lsq_x_fun % casadi.Function: Function for least squares cost on the differential state. Only defined in an OCP.
        f_lsq_u_fun % casadi.Function: Function for least squares cost on the controls. Only defined in an OCP.
        f_lsq_T_fun % casadi.Function: Function for least squares cost on the terminal state. Only defined in an OCP.
    end

    methods(Abstract)
        generate_variables(obj, opts)
        generate_equations(obj, opts)
    end
    % Print dcs output
    methods(Access=protected)
        function propgrp = getPropertyGroups(obj)
            gTitle1 = 'Populated Functions';
            propList1 = struct;
            names = properties(obj);
            for ii=1:length(names)
                name = names{ii};
                if ~endsWith(name, '_fun')
                    continue
                end
                if any(obj.(name).size_out(0) == 0)
                    continue
                end
                % some custom handling for objective functions:
                if strcmp(name, 'f_q_fun')
                    if obj.model.f_q == 0
                        continue
                    end
                end
                if strcmp(name, 'f_q_T_fun')
                    if obj.model.f_q_T == 0
                        continue
                    end
                end
                propList1.(names{ii}) = obj.(name).repr;% TODO(@anton) better custom display here
            end
            propgrp(1) = matlab.mixin.util.PropertyGroup(propList1,gTitle1);
        end

        function displayScalarObject(obj)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            scalarHeader = [className ' DCS'];
            header = sprintf('%s\n',scalarHeader);
            disp(header)
            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
        end
    end
end
