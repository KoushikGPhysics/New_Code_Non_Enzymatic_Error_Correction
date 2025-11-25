
p=linspace(1,100,500); 
q_r=1; q_w=100;
a_r=5000; b_r=1/5000;
a_w=1; b_w=1;
w=10^9; v=0;

s=first_covalent_bonded_s_matrix();

initial_state=1;

final_state=[49,50];

ratio=zeros(1,length(p));

for i=1:length(p)
    [Q]=first_covalent_bonded_t_matrix(s,p(i),q_r,q_w,a_r,b_r,a_w,b_w,w);
    Q_temp=Q;
    Q_temp(final_state,:)=[ ];
    Q_temp(:,final_state)=[ ];
    F=-inv(Q_temp);
    t_all=sum(F,2);
    t_ab=t_all(initial_state);
    [prob_all, prob_ter, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);
    ratio(i)=steady_state(end-1)/steady_state(end);
end
figure;
semilogy(log(p/q_r),ratio,'b-','lineWidth',3);
axis tight;
xlabel('\Delta G_{C} (in units of k_{B}T)\rightarrow');
ylabel('\eta_{K+T}\rightarrow')
ax=gca;
ax.YAxis.Exponent = 3;  
ax.XAxis.FontSize = 11;   
ax.YAxis.FontSize = 11;   
ylim([1800 4000]);
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');
