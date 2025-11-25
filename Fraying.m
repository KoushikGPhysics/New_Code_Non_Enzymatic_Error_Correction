function[]=Fraying()
%% input
p=2;
q_r=1;
q_w=1;
a_r=5000;
b_r=(5000)^(-1);
a_w=1;
b_w=1;
n=5;
E_position=3;
%% main body of the code
s=s_matrix(n,E_position);
initial_state=find(s(:, E_position-1) == 1 & (s(:, E_position) == 1|s(:, E_position) == 2) & all(s(:, E_position+1:end) == 0,2));
final_state=1;

[Q]=non_absorbing_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w);

t=linspace(0,0.0005,5000)./q_r;

P=[];
for i=1:length(initial_state)
    [~,prob_ter]=P_Evolution(t,initial_state(i),final_state,Q);
    P=[P;prob_ter];
end

%% visulaization 
figure;
color1 = [0.0, 0.45, 0.74];   
color2 = [0.8500, 0.3250, 0.0980];  
axis tight;
plot(t, P(1,:), 'b-','LineWidth', 3); 
hold on;
plot(t, P(2,:), 'LineWidth', 3, 'Color', color2);

legend('11100\rightarrow 10000','11200\rightarrow 10000','location','best');
xlabel('t ( in units of 1/q_r )\rightarrow');
ylabel('P(t)\rightarrow');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');     
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
ax=gca;
ax.YAxis.Exponent = -2;  
ylim([0  0.06])
end
