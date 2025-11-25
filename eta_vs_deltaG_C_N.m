clear; clc; close all;

E_position = 3;                   % Error position
q_r = 1;                                % dissociation rate of correct pair
ther_dis_fac = 100;             % thermodynamic discrimination factor
a_r = 5000; b_r = 1/5000;  
q_w = ther_dis_fac * q_r;
a_w = 1;  b_w = 1;

p = linspace(q_r, 200, 1e3);    % p-range
w = 1e6;                                    % thermodynamic discrimination parameter

N_values = [6, 7, 8];                % <<<<<< Multiple n values
colors = {[0, 0, 1],[0, 0.75, 0.75],[0.8500, 0.3250, 0.0980]};            
styles = {'-','-','-'};

figure; hold on;
legend_names = {};

for idx = 1:length(N_values)

    n = N_values(idx);

    %% ===== STATE MATRIX ===== %%
    s = covalent_s_matrix(n, E_position);

    %% ==== INITIAL / FINAL STATE ===== %%
    initial_state = find(s(:,1)==1 & all(s(:,2:end)==0,2));
    final_state   = find(all(s~=0 & s~=1 & s~=2, 2));

    %% ===== PROBABILITY RATIO ===== %%
    prob_ratio = zeros(1, length(p));

    for i = 1:length(p)
        Q = covalent_t_matrix(s, p(i), q_r, q_w, a_r, b_r, a_w, b_w, w);

        Q_temp = Q;
        Q_temp(final_state,:) = [];
        Q_temp(:,final_state) = [];

        F = -inv(Q_temp);
        t_all = sum(F, 2);
        t_ab = t_all(initial_state);

        [~, ~, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);

        prob_ratio(i) = steady_state(end-1) / steady_state(end);
    end

    %% ===== PLOT FOR THIS VALUE OF n ===== %%
    plot(log(p/q_r), prob_ratio, ...
        'Color', colors{idx}, ...
        'LineStyle', styles{idx}, ...
        'LineWidth', 3);

    legend_names{idx} = ['N = ' num2str(n)];
end

%% ===== FINAL FIGURE FORMATTING ===== %%
xlabel('\DeltaG_{C} (in units of k_{B}T)\rightarrow');
axis tight;
ax = gca;
ax.YAxis.Exponent = 3;
ylabel('\eta_{K+T}\rightarrow');
ylim([350 3700]);
legend(legend_names, 'Location','best', 'FontSize', 12);

set(gca, 'Box', 'on', 'LineWidth',1, 'FontName','Times New Roman', 'FontSize', 14);
set(gca, 'TickDir','out', 'TickLength',[0.02 0.02], 'XMinorTick','on', 'YMinorTick','on');

hold off;
