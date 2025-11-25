clear; clc; close all;

E_position = 3;                             % Error position
q_r = 1;                                         % dissociation rate of right base pair
ther_dis_fac = 100;                      % thermodynamic discrimination factor
a_r = 5000; b_r = 1/5000;          % kinetic factors (right)
q_w = ther_dis_fac * q_r;  
a_w = 1; b_w = 1;
p = 2;                            
v = logspace(-20, 15, 500);        

N_values = [6, 7, 8]; 
colors = {[0, 0, 1],[0, 0.75, 0.75],[0.8500, 0.3250, 0.0980]};            
styles = {'-','-','-','-'};


legend_names = {};

figure; hold on;

for idx = 1:length(N_values)

    n = N_values(idx);

    %% ====== COVALENT STATE MATRIX ======%%
    s = covalent_s_matrix(n, E_position);
    initial_state = find(s(:,1)==1 & all(s(:,2:end)==0,2));
    final_state = find(all(s~=0 & s~=1 & s~=2,2));

    %% ====== MFPT ======%%
    mfpt_c_r = zeros(1, length(v));

    for i = 1:length(v)
        Q = covalent_t_matrix(s, p, q_r, q_w, a_r, b_r, a_w, b_w, v(i));
        Q_temp = Q;
        Q_temp(final_state,:) = [];
        Q_temp(:,final_state) = [];

        F = -inv(Q_temp);
        transient_states = setdiff(1:size(Q,1), final_state);
        R = Q(transient_states, final_state);

        B = F * R;
        numerator_C = sum(F(initial_state,:) .* B(:,1)');
        T_C = numerator_C / B(initial_state,1);
        mfpt_c_r(i) = T_C;
    end

    w = v .* mfpt_c_r;

    %% ====== PROBABILITY RATIO =====%%
    prob_ratio = zeros(1,length(w));

    for i = 1:length(w)
        Q = covalent_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w, w(i));
        Q_temp = Q;
        Q_temp(final_state,:) = [];
        Q_temp(:,final_state) = [];
        F = -inv(Q_temp);
        t_all = sum(F,2);

        t_ab = t_all(initial_state);
        [~, ~, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);

        prob_ratio(i) = steady_state(end-1) / steady_state(end);
    end

    %% ====== PLOT ======%%
    loglog(w, prob_ratio, 'Color', colors{idx}, ...
           'LineStyle', styles{idx}, 'LineWidth', 3);

    legend_names{idx} = ['N = ' num2str(n)];
end

%% Formatting
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XLim', [1e-6 1e9]);
ylim([90 2000]);
ylabel('\eta_{K+T}\rightarrow');

yticks([1e2 2e2 5e2 1e3 2e3 5e3 1e4]);
yticklabels({'10^{2}','2×10^{2}','5×10^{2}', '10^{3}','2×10^{3}','5×10^{3}','10^{4}'});

plot([1,1], ylim, ':k', 'LineWidth', 2);

xticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);
xticklabels({'10^{-6}','10^{-3}','10^{0}','10^{3}','10^{6}','10^{9}'});

xlabel('w \rightarrow');

legend(legend_names, 'Location','best', 'FontSize',12);

set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
hold off;
