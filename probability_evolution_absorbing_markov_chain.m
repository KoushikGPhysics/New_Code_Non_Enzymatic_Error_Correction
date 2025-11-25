function[]=probability_evolution_absorbing_markov_chain()
%%% Function to generate figure 5a/6a

Type='K';                                     % Discrimination type('K','T','C')
n=5;                                              % Number of base pairs
E_position=3;                              % Error position
q_r = 1;                                         % dissocaition rate of the right base pair
ther_dis_fac=100;                        % thermodynamic discrimiantion factor
a_r=5000; b_r=1/5000 ;   % Kinetic modulating factors fot the correct base pair

%% ===== STATE MATRIX ======%%%
s=s_matrix(n,E_position);
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
%%% initial state
initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));
%%% final state
final_state=find(all(s~=0,2));
%% ====== PROBABILITY RATIO WITH RESPECT TO TIME =====%%%

Q=absorbing_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w);
t=linspace(10^(-10)/q_r,3,100000);    %% time scale 

%%% Finding the probability
[~, prob_ter] = P_Evolution(t, initial_state, final_state, Q);
probability_ratio=prob_ter(1,:)./prob_ter(2,:);

etaK = probability_ratio(end);
exp_val = floor(log10(etaK));
mantissa = etaK / 10^exp_val;
switch Type
    case 'K'
         formatted_etaK = ['\eta_{K} = ', num2str(mantissa, '%.2f'), ' \times 10^{', num2str(exp_val), '}'];
    case 'T'
         formatted_etaK = ['\eta_{T} = ', num2str(mantissa, '%.2f'), ' \times 10^{', num2str(exp_val), '}'];
    otherwise
        formatted_etaK = ['\eta_{K+T} = ', num2str(mantissa, '%.2f'), ' \times 10^{', num2str(exp_val), '}'];
end

disp(formatted_etaK);

%% ======= Figure ======%%%

figure(1);
semilogy(t,prob_ter(1,:),'b-','lineWidth',3);
hold on;
semilogy(t,prob_ter(2,:),'lineWidth',3,'Color',[0.8500, 0.3250, 0.0980]	);
xlabel('t (in units of 1/q_{r} )\rightarrow');
ylabel('P(t)\rightarrow');
ylim([1e-8 1e1]); 
xlim([0 3])
yticks([1e-6 1e-4 1e-2 1e0 ]);
yticklabels({'10^{-6}','10^{-4}' '10^{-2}','10^{0}'});
legend('P_{R}(t)','P_{W}(t)','Location', 'Best');
set(gca, 'Box', 'on', 'LineWidth', 1, 'FontSize', 11, 'FontName', 'Times New Roman');

set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');
     

 end
