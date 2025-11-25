
function[]=stalling_vs_delta_G_C()

Type='C';

switch Type
    case 'K'
        p=linspace(1,2.5,1000);
    case 'T'
        p=linspace(1,50,1000);
    case 'C'
        p=linspace(1,50,1000);
end


q_r=1;                                                         % dissocation rate od incorrect base pair
n=5;                                                            % Number of base pairs
E_position=3;                                             % Error position                                                            
ther_dis_fac=100;                                      % thermodynamic discrimiantion factor
a_r=5000; b_r=1/5000;                 % Kinetic modulating factors fot the correct base pair

w=10^6;

%% ==== conditional mfpt =====%%
ratio=zeros(1,length(p));

for i=1:length(p)
     [T_C,T_W]=mfpt_ratio(n,E_position,ther_dis_fac,a_r,b_r,p(i),q_r,w,Type);
     ratio(i)=T_W/T_C;
end

figure;
plot(log(p/q_r),ratio,'b-','linewidth',3);
axis tight;
xlabel('\DeltaG_{C}(in units of k_{B}T)\rightarrow');
switch Type
    case 'K'
        ylabel('\tau_{K}\rightarrow');
        set(gca, 'YLim', [26.5 31.5]);
    case 'T'
       ylabel('\tau_{T}\rightarrow');
    otherwise
       ylabel('\tau_{K+T}\rightarrow');
       set(gca, 'YLim', [118 250]);
end
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');
end


%% function to calculate the condtional mfpt ===%%
function [T_C,T_W]=mfpt_ratio(n,E_position,ther_dis_fac,a_r,b_r,p,q_r,w,Type)
%% ===== STATE MATRIX ======%%%
s=covalent_s_matrix(n,E_position);

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

%% ===== Intial and final state index ======%%

%  initial state index
initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));

%  final state index
final_state=find(all(s~=0,2) & all(s~=1,2) & all(s~=2,2));


%% =======PROBABILITY EVOLUTI ======%%%
Q = covalent_t_matrix(s, p, q_r, q_w, a_r, b_r, a_w, b_w,w);
Q_temp = Q;
Q_temp(final_state, :) = [];
Q_temp(:, final_state) = [];


% ===== Extract transient states =====
transient_states = setdiff(1:size(Q,1), final_state);
R = Q(transient_states,final_state); 
% ===== Fundamental matrix =====
F = -inv(Q_temp);

B = F * R; 

% % Step 5: Conditional absorption time into correct state
numerator_C = sum(F(initial_state, :) .* B(:, 1)');           % dot product
T_C = numerator_C / (B(initial_state, 1) + eps); 

% % Step 6: Conditional absorption time into wrong state
numerator_W = sum(F(initial_state, :) .* B(:, 2)');        % dot product
T_W = numerator_W / (B(initial_state, 2) + eps);
end
