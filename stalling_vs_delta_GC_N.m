clear; clc; close all;

q_r = 1;
E_position = 3;
ther_dis_fac = 100;
a_r = 5000; 
b_r = 1/5000;
p = linspace(1,400,2*10^3);
w = 1e6;

N_values = [6, 7, 8];     % <<<<<< MULTIPLE n VALUES
colors = {[0, 0, 1],[0, 0.75, 0.75],[0.8500, 0.3250, 0.0980]};            
styles = {'-','-','-'};
legend_names = {};

figure; hold on;

for idx = 1:length(N_values)

    n = N_values(idx);

    ratio = zeros(1,length(p));

    for i = 1:length(p)
        [T_C, T_W] = mfpt_ratio(n, E_position, ther_dis_fac, a_r, b_r, p(i), q_r, w);
        ratio(i) = T_W / T_C;
    end

    plot(log(p/q_r), ratio, ...
         'Color', colors{idx}, ...
         'LineStyle', styles{idx}, ...
         'LineWidth', 3);

    legend_names{idx} = ['N = ' num2str(n)];
end

%% ===== Formatting ===== %%
xlabel('\DeltaG_{C} (in units of k_{B}T) \rightarrow');
ylabel('\tau_{K+T} \rightarrow');
axis tight;
set(gca, 'YLim', [35 200]);
ax = gca;
ax.YAxis.Exponent = 2;
legend(legend_names, 'Location', 'best', 'FontSize', 12);

set(gca, 'Box', 'on', 'LineWidth', 1, 'FontSize', 14, ...
         'FontName', 'Times New Roman', 'TickDir', 'out', ...
         'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on');

hold off;

function [T_C,T_W]=mfpt_ratio(n,E_position,ther_dis_fac,a_r,b_r,p,q_r,w)
s=covalent_s_matrix(n,E_position);
q_w = ther_dis_fac * q_r;  
a_w = 1; b_w= 1;

initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));
final_state = find(all(s~=0,2) & all(s~=1,2) & all(s~=2,2));

Q = covalent_t_matrix(s, p, q_r, q_w, a_r, b_r, a_w, b_w, w);
Q_temp = Q;
Q_temp(final_state, :) = [];
Q_temp(:, final_state) = [];

transient_states = setdiff(1:size(Q,1), final_state);
R = Q(transient_states,final_state); 
F = -inv(Q_temp);

B = F * R; 

numerator_C = sum(F(initial_state,:) .* B(:,1)');
T_C = numerator_C / (B(initial_state,1) + eps);

numerator_W = sum(F(initial_state,:) .* B(:,2)');
T_W = numerator_W / (B(initial_state,2) + eps);
end
