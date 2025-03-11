% Load the matrices from files
actual = casadi.Sparsity.from_file('debug_fatrop_actual.mtx');
A = casadi.Sparsity.from_file('debug_fatrop_A.mtx');
B = casadi.Sparsity.from_file('debug_fatrop_B.mtx');
C = casadi.Sparsity.from_file('debug_fatrop_C.mtx');
D = casadi.Sparsity.from_file('debug_fatrop_D.mtx');
I = casadi.Sparsity.from_file('debug_fatrop_I.mtx');
errors = row(casadi.Sparsity.from_file('debug_fatrop_errors.mtx'));

% Create a new figure
figure;

% Plot the sparsity patterns
[iA,jA] = A.get_triplet;iA=iA+1;jA=jA+1; % index-one correction
[iB,jB] = B.get_triplet;iB=iB+1;jB=jB+1;
[iC,jC] = C.get_triplet;iC=iC+1;jC=jC+1;
[iD,jD] = D.get_triplet;iD=iD+1;jD=jD+1;
[iI,jI] = I.get_triplet;iI=iI+1;jI=jI+1;

[iActual,jActual] = actual.get_triplet;iActual=iActual+1;jActual=jActual+1;

% Plotting with custom markers
plot(jA, iA, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
hold on;
plot(jB, iB, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
plot(jC, iC, 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
plot(jD, iD, 'co', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
plot(jI, iI, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
spy(sparse(casadi.DM(actual,1)),'ko',2)

% Highlight offending rows
for idx = 1:length(errors)% index-one correction
    yline(errors(idx)+1, 'Color', '#999', 'LineStyle', '-');
end

% Add title and legend
title('Debug view of fatrop interface structure detection');
legend('expected A','expected B','expected C','expected D','expected I','actual','offending rows');

% Display the plot
hold off;
