

clear; clc; close all;

% --- Fixed parameters ---%
q_r = 1;
q_w = 100;
a_r_values = logspace(0,5,5000);
a_w = 1;
b_w = 1;
n = 5;
E_position = 3;

p = 5;   % change range as needed
MFPT_ratio = zeros(size(a_r_values));

s = s_matrix(n, E_position);

initial_state = find(s(:, E_position-1) == 1 & ...
                    (s(:, E_position) == 1 | s(:, E_position) == 2) & ...
                     all(s(:, E_position+1:end) == 0, 2));

final_state = find(s(:, E_position-1) == 1 & ...
                   all(s(:, E_position:end) == 0, 2));

for k = 1:length(a_r_values)
    a_r = a_r_values(k);

    % Transition rate matrix (non-absorbing)
    Q = non_absorbing_t_matrix(s, p, q_r, q_w, a_r, 1/a_r, a_w, b_w);

    % Remove absorbing (final) state rows/columns
    Q(final_state, :) = [];
    Q(:, final_state) = [];

    % Fundamental matrix
    F = -inv(Q);

    % Mean first-passage times
    T = sum(F, 2);
    MFPT_ratio(k) = T(12) / T(20);
end

% --- Plot  ---%
figure;
loglog(a_r_values, MFPT_ratio, 'b-','LineWidth', 3);
xlabel('\alpha\rightarrow','FontSize',18);
ylabel('MFPT_{11100\rightarrow 11000}/MFPT_{11200\rightarrow 11000}','FontSize',12);
hold on;
target_a_r = 5e3;
[~, idx] = min(abs(a_r_values - target_a_r));  
plot(a_r_values(idx), MFPT_ratio(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');
ylim([10.283  1500]);
