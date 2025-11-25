
function[]=eta_as_function_of_w()

%%% Function to genarate the figure 7

Type='C';                                    % Discrimination type('K','T','C')
n=5;                                            % Number of base pairs
E_position=3;                            % Error position
q_r = 1;                                       % dissocaition rate of the right base pair
ther_dis_fac=100;                       % thermodynamic discrimiantion factor
a_r=5000; b_r=1/5000;  % Kinetic modulating factors fot the correct base pair


%% ==== THERMODYNAMIC DISCRIMINATION ===== %%%
% Dissociation rate of ther incorrect base pairs according to the
% 'ther_dis_fac' according to the Type
switch Type
    case 'K'  
        q_w = q_r;
    case 'T'  
        q_w = ther_dis_fac * q_r; 
    case 'C'  
        q_w = ther_dis_fac * q_r;  
    otherwise
        error('Unknown Type. Use ''K'', ''T'', or ''C''.');
end
% base pair formation rate 
p = 2;      

%% ==== KINETIC DISCRIMINATION=====%%%
% Setting the kinetic modulating factor for the incorrect base pair
% according to the Type
switch Type
    case 'K'  
        a_w = 1;
        b_w= 1;
    case 'T'  
        a_w = a_r;  
        b_w = b_r;
    case 'C'  
        a_w = 1;
        b_w= 1;
    otherwise
        error('Unknown Type. Use ''K'', ''T'', or ''C''.');
end

%% ====== COVALENT BOND STATE MATRIX ======%%%

s=covalent_s_matrix(n,E_position);


%% ==== INTIAL STATE / FINAL STATE INDEX =====%%
initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));

final_state=find(all(s~=0 & s~=1 & s~=2,2));

%% ------- convalent bond formation rate -----%%
v=logspace(-20,15,500);

%% ====== MFPT =====%%%
mfpt_c_r=zeros(1,length(v));

for i=1:length(v)
    Q = covalent_t_matrix(s, p, q_r, q_w, a_r, b_r, a_w, b_w,v(i));
    % % Q=absorbing_t_matrix(s, p, q_r, q_w, a_r, b_r, a_w, b_w);
    Q_temp = Q;
    Q_temp(final_state, :) = [];
    Q_temp(:, final_state) = [];


    %%===== Extract transient states =====%%
    transient_states = setdiff(1:size(Q,1), final_state);
    R = Q(transient_states,final_state); 
    %%%===== Fundamental matrix =====%%
    warning('off', 'MATLAB:nearlySingularMatrix');
    F = -inv(Q_temp);
    B = F * R; 


    numerator_C = sum(F(initial_state, :) .* B(:, 1)');           
    T_C = numerator_C / (B(initial_state, 1)); 
    mfpt_c_r(i)=T_C;
end


w=v.*mfpt_c_r;

% % % %%  ====== PROBABILITY RATIO ======%%%
prob_ratio=zeros(1,length(w));
for i=1:length(w)
    Q=covalent_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w, w(i));
    Q_temp=Q;
    Q_temp(final_state,:)=[ ];
    Q_temp(:,final_state)=[ ];
    F=-inv(Q_temp);
    t_all=sum(F,2);
    t_ab=t_all(initial_state);
    [~, ~, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);
    prob_ratio(i)=(steady_state(end-1) ./ steady_state( end));
end
%% ====== VISUALIZATION ======%%%% % % 
figure(1);
loglog(w, prob_ratio,'b-','LineWidth', 3);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XLim', [1e-6 1e9]);
ax=gca;
switch Type
    case 'K'
        set(gca, 'YLim', [0.6 1000]);
        ylabel('\eta_{K}\rightarrow');
    case 'T'
        set(gca, 'YLim', [40 110]);
        ylabel('\eta_{T}\rightarrow');
        ax.YAxis.Exponent = 2; 
    case 'C'
        set(gca, 'YLim', [70 4000]); 
        ylabel('\eta_{K+T}\rightarrow');
        yticks([1e2 2e2 5e2 1e3 2e3 5e3 1e4]);
       yticklabels({'10^{2}','2×10^{2}','5×10^{2}', ...
             '10^{3}','2×10^{3}','5×10^{3}','10^{4}'});
end
hold on;

plot([1, 1], ylim, ':', 'Color', 'k', 'LineWidth', 3);

xticks([1e-6 1e-3 1e0 1e3 1e6 1e9]);
xticklabels({'10^{-6}', '10^{-3}', '10^{0}', '10^{3}', '10^{6}','10^{9}'});
xlabel('w\rightarrow');
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');
hold off;

end