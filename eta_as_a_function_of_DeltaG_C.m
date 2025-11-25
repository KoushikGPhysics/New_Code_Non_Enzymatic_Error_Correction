
function[]=eta_as_a_function_of_DeltaG_C()
%%% function to generate the figure 8

Type='T';                                      % Discrimination type('K','T','C')
n=5;                                               % Number of base pairs
E_position=3;                               % Error position
q_r = 1;                                          % dissocaition rate of the right base pair
ther_dis_fac=100;                         % thermodynamic discrimiantion factor
a_r=5000; b_r=1/5000 ;    % Kinetic modulating factors fot the correct base pair

%% ===== STATE MATRIX ======%%%
s=covalent_s_matrix(n,E_position);

%% ==== THERMODYNAMIC DISCRIMINATION ===== %%%

switch Type
    case 'K'  
        q_w = q_r;
         p = linspace(q_r,4,1000);
         a_w = 1;
         b_w= 1;
    case 'T'  
        q_w = ther_dis_fac * q_r;  
        p = linspace(q_r,50,1000); 
         a_w = a_r;  
        b_w = b_r;
    case 'C'  
        q_w = ther_dis_fac * q_r; 
         p = linspace(q_r,500,10^4);
         a_w = 1;
        b_w= 1;
    otherwise
        error('Unknown Type. Use ''K'', ''T'', or ''C''.');
end

w=10^6;
 
%% ==== INTIAL STATE / FINAL STATE INDEX =====%%

initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));

final_state=find(all(s~=0 & s~=1 & s~=2,2));


%%  ====== PROBABILITY RATIO ======%%%
prob_ratio=zeros(1,length(p));

for i=1:length(p)
    Q=covalent_t_matrix(s,p(i),q_r,q_w,a_r,b_r,a_w,b_w, w);
    Q_temp=Q;
    Q_temp(final_state,:)=[ ];
    Q_temp(:,final_state)=[ ];
    F=-inv(Q_temp);
    t_all=sum(F,2);
    t_ab=t_all(initial_state);
    [~, ~, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);
    prob_ratio(i)=(steady_state( end-1) ./ steady_state( end));
end

%% ====== VISUALIZATION ======%%%
figure;
plot(log(p/q_r),prob_ratio,'b-','lineWidth',3);
axis tight;
xlabel('\DeltaG_{C} (in units of k_{B}T)\rightarrow');

ax = gca;

switch Type
    case 'K'
        ax.YAxis.Exponent = 2;  
        ylabel('\eta_{K}\rightarrow');
        ylim([575 710]);
    case 'T'
        ax.YAxis.Exponent = 0;  
        ylabel('\eta_{T}\rightarrow');
    otherwise
        ax.YAxis.Exponent = 3;  
        ylabel('\eta_{K+T}\rightarrow');
        ylim([1800 4000]);
end
set(gca, 'Box', 'on', 'LineWidth',1, 'FontSize', 12, 'FontName', 'Times New Roman');  
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'XColor', 'k', 'YColor', 'k', 'FontSize', 14, 'FontName', 'Times New Roman');   
set(gcf, 'PaperUnits', 'inches');

end
