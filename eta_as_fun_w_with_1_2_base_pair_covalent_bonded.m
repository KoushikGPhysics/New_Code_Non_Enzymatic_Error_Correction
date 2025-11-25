
p=2; 
q_r=1; q_w=100;
a_r=5000; b_r=1/5000;
a_w=1; b_w=1;
w=logspace(-17,15,500); 


s=first_covalent_bonded_s_matrix();

initial_state=1;

final_state=[49,50];


mfpt_c_r=zeros(1,length(w));

for i=1:length(w)
    [Q]=first_covalent_bonded_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w,w(i));
    Q_temp=Q;
    Q_temp(final_state,:)=[];
    Q_temp(:,final_state)=[];
    % ===== Extract transient states =====%%
    transient_states = setdiff(1:size(Q,1), final_state);
    R = Q(transient_states,final_state); 
    % ===== Fundamental matrix =====%%
     warning('off', 'MATLAB:nearlySingularMatrix');
    F = -inv(Q_temp);
    B = F * R; 


    numerator_C = sum(F(initial_state, :) .* B(:, 1)');           
    T_C = numerator_C / (B(initial_state, 1) + eps); 
    mfpt_c_r(i)=T_C;
end

ratio=zeros(1,length(w));

w_w=w.*mfpt_c_r;

for i=1:length(w)
    [Q]=first_covalent_bonded_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w,w_w(i));
    Q_temp=Q;
    Q_temp(final_state,:)=[ ];
    Q_temp(:,final_state)=[ ];
    F=-inv(Q_temp);
    t_all=sum(F,2);
    t_ab=t_all(initial_state);
    [prob_all, prob_ter, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);
    ratio(i)=steady_state(end-1)/steady_state(end);
end

loglog(w_w, ratio,'b-','LineWidth', 3);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XLim', [1e-5 1e8]);
axis tight;
set(gca, 'YLim', [70 4000]); 
        ylabel('\eta_{K+T}\rightarrow');
        yticks([1e2 2e2 5e2 1e3 2e3 5e3 1e4]);
       yticklabels({'10^{2}','2×10^{2}','5×10^{2}', ...
             '10^{3}','2×10^{3}','5×10^{3}','10^{4}'});
hold on;

plot([1, 1], ylim, ':', 'Color', 'k', 'LineWidth', 3);

xlabel('w\rightarrow');
ax = gca;
ax.XAxis.FontSize = 11;   
ax.YAxis.FontSize = 11;   
hold off;
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');


