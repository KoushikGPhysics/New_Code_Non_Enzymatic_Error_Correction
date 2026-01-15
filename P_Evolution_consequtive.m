
Type='K';                                     % Discrimination type('K','T','C')
n=5;                                              % Number of base pairs
E_position=2;                              % Error position
q_r = 1;                                         % dissocaition rate of the right base pair
ther_dis_fac=100;                        % thermodynamic discrimiantion factor
a_r=5000; b_r=1/5000 ;             % Kinetic modulating factors fot the correct base pair

%% ===== STATE MATRIX ======%%%
s=s_matrix_full(n);

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
p =2;      

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

%% ==== INTIAL STATE / FINAL STATE INDEX =====%%

initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));
final_state=find(all(s~=0,2));

%% ====== PROBABILITY RATIO WITH RESPECT TO TIME =====%%%

Q=absorbing_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w);
t=linspace(0,3,5000);    

[~, prob_ter] = P_Evolution(t, initial_state, final_state, Q);

figure(1);

% -------- Main plot --------%
semilogy(t, prob_ter(5,:), '-', ...
         'Color', '#50F', ...
         'LineWidth', 2.2);
hold on;

semilogy(t, prob_ter(7,:), '-', ...
         'Color', '#F00', ...
         'LineWidth', 2.2);

xlabel('t (in units of 1/q_{r})\rightarrow');
ylabel('P(t)\rightarrow');
ylim([0 1.2*1e-3]);

legend('11211','11221','Location','best');

set(gca, 'Box', 'on', 'LineWidth', 0.8, ...
         'TickDir', 'out', 'TickLength', [0.02, 0.02], ...
         'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'FontSize', 12, 'FontName', 'Times New Roman');

     
   
% % % % % % ---- For Inset axes (Figure 2 inside Figure 1) ----%%%%%

% % % % % ax_inset = axes('Position',[0.55 0.55 0.35 0.35]);  
% % % % % % [x y width height] — normalized to figure
% % % % % 
% % % % % plot(t, prob_ter(5,:), '-', ...
% % % % %          'Color', '#50F', ...
% % % % %          'LineWidth', 2.2);
% % % % % hold on;
% % % % % 
% % % % % plot(t, prob_ter(7,:), '-', ...
% % % % %          'Color', '#F00', ...
% % % % %          'LineWidth', 2.2);
% % % % % 
% % % % % xlim([0 0.12])
% % % % % ylim([0 1e-3])
% % % % % 
% % % % % set(ax_inset, 'Box', 'on', ...
% % % % %                'TickDir', 'out', ...
% % % % %                'FontSize', 10, ...
% % % % %                'FontName', 'Times New Roman');
