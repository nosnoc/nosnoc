function [varargout] = unfold_struct(input,mode)
% mode = 'caller'; % copy fields locally in the function where unfold_struct was called
% mode = 'base'; % copy fields  into the main workspace

%% 
import casadi.*
names = fieldnames(input);
for iii=1:length(names)
    eval([names{iii} '=input.' names{iii} ';']);
%     if eval(['~exist( '''  names{i} ''', ''var'')'])
        %         Check if a variable exists in the workspace, within a function
        assignin(mode, names{iii}, eval(names{iii}));
%     end
end
end




% for i=1:length(names)
%     eval([names{i} '=input.' names{i} ]);
%     eval(['varargout{' num2str(i) '} =' names{i}])
% end
