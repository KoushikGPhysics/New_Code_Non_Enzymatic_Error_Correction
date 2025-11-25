
n = 5;                                        % Number of base pairs
E_position = 3;                        % Error position
q_r = 1;                                    % Dissociation rate of the right base pair
ther_dis_fac = 100;                 % Thermodynamic discrimination factor
q_w = ther_dis_fac * q_r;
p = 2;

a_w = 1; 
b_w = 1;

s = covalent_s_matrix(n, E_position);

initial_state = find(s(:,1) == 1 & all(s(:,2:end) == 0, 2));
final_state = find(all(s ~= 0 & s ~= 1 & s ~= 2, 2));


v = logspace(-15, 20, 100);                 
a_r_values = logspace(1, 4, 100);    

prob_ratio_mat = zeros(length(a_r_values), length(v));
w_mat = zeros(length(a_r_values), length(v));

for j = 1:length(a_r_values)
    a_r = a_r_values(j);
    mfpt_c_r = zeros(1, length(v));

    for i = 1:length(v)
        % --- Transition Matrix ---%
        Q = covalent_t_matrix(s, p, q_r, q_w, a_r, 1/a_r, a_w, b_w, v(i));
        Q_temp = Q;
        Q_temp(final_state,:) = [];
        Q_temp(:,final_state) = [];
        
        warning('off', 'MATLAB:nearlySingularMatrix');
        if rcond(Q_temp) < 1e-12
            F = -pinv(Q_temp);   
        else
            F = -inv(Q_temp);
        end

        % --- Absorption Probability ---%
        R = Q(setdiff(1:size(Q,1), final_state), final_state);
        B = F * R;

        % ----Conditional absorption time into correct state----%
        numerator_C = sum(F(initial_state, :) .* B(:, 1)');           
        T_C = numerator_C / (B(initial_state, 1) + eps); 
        mfpt_c_r(i)=T_C;
    end

    % Effective rate
    w = v .*mfpt_c_r;
    w_mat(j, :) = w;

    % ---- Probability Ratio Calculation ----%
    for i = 1:length(w)
        Q = covalent_t_matrix(s, p, q_r, q_w, a_r, 1/a_r, a_w, b_w, w(i));
        Q_temp = Q;
        Q_temp(final_state,:) = [];
        Q_temp(:,final_state) = [];
        F = -inv(Q_temp);
        t_all = sum(F, 2);
        t_ab = t_all(initial_state);
        [~, ~, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);
        prob_ratio_mat(j,i) = steady_state(end-1) ./ steady_state(end);
    end
end


% Interpolate each row of prob_ratio_mat to a common w grid (log–log)

w_common = logspace(-9, 9, 100);  

% Initialize matrix for interpolated values
prob_ratio_interp = nan(length(a_r_values), length(w_common));

for j = 1:length(a_r_values)
    w_row = w_mat(j, :);
    pr_row = prob_ratio_mat(j, :);

    % Remove invalid points
    valid = ~isnan(w_row) & ~isnan(pr_row) & (w_row > 0) & (pr_row > 0);
    if nnz(valid) < 3
        continue
    end

    % Sort and prepare for interpolation
    [w_sorted, idx] = sort(w_row(valid));
    pr_sorted = pr_row(valid);
    pr_sorted = pr_sorted(idx);

    % Interpolate in log–log space for numerical stability
    logw_src = log10(w_sorted);
    logpr_src = log10(pr_sorted);
    logw_tgt = log10(w_common);

    % Linear interpolation (with extrapolation)
    logpr_tgt = interp1(logw_src, logpr_src, logw_tgt, 'linear', 'extrap');

    % Back-transform to normal scale
    prob_ratio_interp(j, :) = 10.^logpr_tgt;
end


figure;
pcolor(log10(w_common), log10(a_r_values), log10(prob_ratio_interp));
shading interp;
set(gca, 'YDir', 'normal');
xlabel('log_{10}(w)\rightarrow');
ylabel('log_{10}(\alpha)\rightarrow');
cb = colorbar;
ylabel(cb, 'log_{10}(\eta_{K+T})\rightarrow');
set(gca, 'FontSize', 14, ...
         'LineWidth', 1, ...
         'TickDir', 'out', ...
         'Box', 'on', ...
         'Layer', 'top', ...
         'TickLabelInterpreter', 'latex');
colormap(jet);
xlim([-9 9]);
ylim([min(log10(a_r_values)) max(log10(a_r_values))]);

