
function[]=stalling_as_fun_delta_G_C_with_1_2_base_cov_bonded()


q_r=1;                                                         % dissocation rate od incorrect base pair
n=5;                                                            % Number of base pairs
E_position=3;                                            % Error position                                                            
ther_dis_fac=100;                                      % thermodynamic discrimiantion factor
a_r=5000^1; b_r=5000^(-1) ;                 % Kinetic modulating factors fot the correct base pair
p=linspace(1,50,1000);
w=10^9;

ratio=zeros(1,length(p));

for i=1:length(p)
     [T_C,T_W]=mfpt_ratio(n,E_position,ther_dis_fac,a_r,b_r,p(i),q_r,w);
     ratio(i)=T_W/T_C;
end


figure(1);
plot(log(p/q_r),ratio,'b-','linewidth',3);
axis tight;
xlabel('\DeltaG_{C} (in units of k_{B}T)\rightarrow');
ylabel('\tau_{K+T}\rightarrow');
ax=gca;
ax.XAxis.FontSize = 11;   
ax.YAxis.FontSize = 11;
ylim([118 250]);
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');
end


function [T_C,T_W]=mfpt_ratio(~,~,ther_dis_fac,a_r,b_r,p,q_r,w)
s=first_covalent_bonded_s_matrix();
q_w = ther_dis_fac * q_r; 
  a_w = 1; b_w= 1;

%% ===== Intial and final state index ======%%

%  initial state index
initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));

%  final state index
final_state=find(all(s~=0,2) & all(s~=1,2) & all(s~=2,2));



[Q]=first_covalent_bonded_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w,w);
Q_temp = Q;
Q_temp(final_state, :) = [];
Q_temp(:, final_state) = [];

transient_states = setdiff(1:size(Q,1), final_state);
R = Q(transient_states,final_state); 

F = -inv(Q_temp);

B = F * R; 

numerator_C = sum(F(initial_state, :) .* B(:, 1)');           
T_C = numerator_C / (B(initial_state, 1) + eps); 

numerator_W = sum(F(initial_state, :) .* B(:, 2)');       
T_W = numerator_W / (B(initial_state, 2) + eps);
end
