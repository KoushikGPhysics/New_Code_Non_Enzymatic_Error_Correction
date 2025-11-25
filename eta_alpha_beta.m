
function[]=eta_alpha_beta()
%% ===== input ========%%%
Type='K';                                  % Discrimination type('K','T','C')
n=5;                                           % Number of base pairs
E_position=3;                           % Error position
q_r = 1;                                      % dissocaition rate of the right base pair
ther_dis_fac=100;                     % thermodynamic discrimiantion factor
a_r=logspace(-2,3,50);
b_r=logspace(-2,3,50) ;         % Kinetic modulating factors fot the correct base pair

s=unidirectional_s_matrix(n,E_position);
switch Type
    case 'K'  
        q_w = q_r;
         a_w = 1;
        b_w= 1;
    case 'T'  
        q_w = ther_dis_fac * q_r;  
        a_w = a_r;  
        b_w = b_r;
    case 'C'  
        q_w = ther_dis_fac * q_r;
         a_w = 1;
        b_w= 1;
    otherwise
        error('Unknown Type. Use ''K'', ''T'', or ''C''.');
end
% base pair formation rate 
p = 2;      

initial_state = find(all(s(:, :) == 0, 2));
final_state=find(all(s~=0,2));

probability_ratio= zeros(length(b_r),length(a_r));

for i=1:length(b_r)
    for j=1:length(a_r)
        Q=absorbing_t_matrix(s,p,q_r,q_w,a_r(j),b_r(i),a_w,b_w);
        Q_temp=Q;
        Q_temp(final_state,:)=[ ];
        Q_temp(:,final_state)=[ ];
        F=-inv(Q_temp);
        t_all=sum(F,2);
        t_ab=t_all(initial_state);
        % to find the steady state we use the 'absorbing_time*10^8
        [~, prob_ter] = P_Evolution(t_ab*10^8, initial_state, final_state, Q);
        probability_ratio(i,j) = prob_ter(1, end) / prob_ter(2, end);
    end
end
%% ======= VISUALIZATION =========%%%
figure;
pcolor(log10(a_r), log10(b_r), (probability_ratio));
shading interp 
xlabel('\alpha\rightarrow');
ylabel('\beta\rightarrow')
x_axis = log10(a_r);
min_x_axis = floor(min(x_axis));
max_x_axis = ceil(max(x_axis));
xticks = min_x_axis:max_x_axis;
xtick_labels = arrayfun(@(x) sprintf('10^{%d}', x), xticks, 'UniformOutput', false);
set(gca, 'XTick', xticks, 'XTickLabel', xtick_labels);
y_axis = log10(b_r);
min_y_axis = floor(min(y_axis));
max_y_axis = ceil(max(y_axis));
yticks = min_y_axis:max_y_axis;
ytick_labels = arrayfun(@(x) sprintf('10^{%d}', x), yticks, 'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', ytick_labels);
axis xy;
c=colorbar;
switch Type
    case 'K'
         ylabel(c,'\eta_{K}\rightarrow');
    case 'T'
       ylabel(c,'\eta_{T}\rightarrow');
    otherwise
        ylabel(c,'\eta_{K+T}\rightarrow')
end
set(gca, 'Box', 'on', 'LineWidth',0.8, 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');  
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
end

