clear; clc;

Type='K';                                   % Discrimination type('K','T','C')
n=5;                                            % Number of base pairs
E_position=3;                            %  Error position
q_r = 1;                                       % dissocaition rate of the right base pair
ther_dis_fac=100;                      % thermodynamic discrimiantion factor
a_r=5000; b_r=1/5000 ;           % Kinetic modulating factors fot the correct base pair


switch Type
    case 'K'
        q_w = q_r;
        p = linspace(1,5,5000);   
        a_w = 1;
        b_w= 1;
    case 'T'
        q_w = ther_dis_fac * q_r;  
        p = linspace(1,100,500);   
        a_w = a_r;  
        b_w = b_r;
    case 'C'
        q_w = ther_dis_fac * q_r;  
        p = linspace(1,100,5000);
         a_w = 1;
        b_w= 1;
end
 
w=10^5; 

s = covalent_s_matrix(n, E_position);  

S_dot=zeros(1,length(p));

for i=1:length(p)
    Q = covalent_t_matrix(s, p(i), q_r, q_w, a_r, b_r, a_w, b_w, w); 
    [~,S_dot(i),~,~]=free_energy_dissipation_per_monomer_incorporation(s,Q);
end
initial_state = find(s(:,1) == 1 & all(s(:, 2:end) == 0, 2));

final_state=find(all(s~=0 & s~=1 & s~=2,2));

prob_ratio=zeros(1,length(p));

for i=1:length(p)
    Q=covalent_t_matrix(s,p(i),q_r,q_w,a_r,b_r,a_w,b_w, w);
    Q_temp=Q;
    Q_temp(final_state,:)=[ ];
    Q_temp(:,final_state)=[ ];
    F=-inv(Q_temp);
    t_all=sum(F,2);
    t_ab=t_all(initial_state);
    [prob_all, prob_ter, steady_state] = steady_state_probability(t_ab, initial_state, final_state, Q);
    prob_ratio(i)=(steady_state(end-1) ./ steady_state(end));
end

figure(1); clf;
loglog(S_dot, prob_ratio,'b-', 'LineWidth', 3);  
xlabel('EPR (sec^{-1})\rightarrow');
switch Type
    case 'K'
       ylabel('\eta_{K}\rightarrow');
    case 'T'
         ylabel('\eta_{T}\rightarrow');
    otherwise
        ylabel('\eta_{K+T}\rightarrow');
end
axis tight;

set(gca, 'LineWidth', 1.3, 'FontSize', 14, 'FontName', 'Times', ...
    'Box', 'on', 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on');
ax=gca;
ax.YAxis.Exponent = 3;  
switch Type
    case 'K'
       ax.YAxis.Exponent = 2; 
        ylim([570 720]);
    case 'T'
      ax.YAxis.Exponent = 0;  yticks([2   10  50 ]) 
    case 'C'
       ax.YAxis.Exponent = 3; 
% % %        ylim([1840 4000]);
end

function [G, Sdot, Jtot, flux_details] = free_energy_dissipation_per_monomer_incorporation(s, Q)
% ======================================================================
% Computes the free-energy dissipation per monomer incorporation (DeltaG)
% Based on: Shu et al., J. Phys.: Condens. Matter 27, 235105 (2015)
%
% Inputs:
%   s   : State matrix (each row = system microstate)
%   Q   : Rate matrix (generator of continuous-time Markov process)
%
% Outputs:
%   G             : Free-energy dissipation per incorporation (?G)
%   Sdot          : Total entropy production rate (steady state)
%   Jtot          : Total productive incorporation flux (to absorbing states)
%   flux_details  : Structure containing detailed flux components
%
% ----------------------------------------------------------------------
%   where A = absorbing (product) states
% ======================================================================

%% --- 1. Identify absorbing and transient states ---%
abs_idx  = find(all(s ~= 0 & s ~= 1 & s ~= 2, 2));    % absorbing (product) states
trans_idx = setdiff(1:size(s,1), abs_idx);                       % all other (transient) states

Q_TT = Q(trans_idx, trans_idx);                                % transitions among transient states
Q_TA = Q(trans_idx, abs_idx);                                 % transitions from transient -> absorbing

%% --- 2. Solve for steady-state distribution in transient subnetwork ---%

A = [Q_TT'; ones(1, size(Q_TT,1))];
b = [zeros(size(Q_TT,1),1); 1];
Pss = A\b;
Pss = max(Pss, 0);               
Pss = Pss / sum(Pss);            
Pss = Pss(:);                    

%% --- 3. Compute microscopic fluxes ---%
piQ = bsxfun(@times, Pss, Q_TT);
Jmat = piQ - piQ';               

%% --- 4. Compute entropy production rate (Sdot) ---%
Sdot = 0;
for i = 1:size(Jmat,1)
    for j = i+1:size(Jmat,1)
        pij = piQ(i,j);
        pji = piQ(j,i);
        if pij > 0 && pji > 0
            Sdot = Sdot + (pij - pji) * log(pij / pji);
        end
    end
end
% Units: k_B T per unit time (dimensionless if rates are scaled)

%% --- 5. Compute fluxes ---%%
% (a) Productive incorporation flux (transient -> absorbing)
Jtot = sum(Pss' * Q_TA);

% (b) Total internal cyclic flux within transient states
J_internal = 0.5 * sum(sum(abs(Jmat)));

% (c) Combined flux including all transitions
J_total_all_steps = Jtot + J_internal;

%% --- 6. Free energy dissipation per incorporation ---%%
G = Sdot / J_total_all_steps;

%% --- 7. Output detailed results ---%%
flux_details = struct();
flux_details.Jmat = Jmat;
flux_details.Pss = Pss;
flux_details.abs_idx = abs_idx;
flux_details.trans_idx = trans_idx;
flux_details.J_internal = J_internal;
flux_details.Jtot = Jtot;
flux_details.J_total_all_steps = J_total_all_steps;
flux_details.Sdot = Sdot;

end




function[Q]=covalent_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w, w)


%%% Function for transition matrix for reversible Markov chain consist of
%%% covalent bonding

%%% ==== Input =====%%%

% s= state matrix (s=covalent_s_matrix(n,E_positon)
%)
% p=base pair formation rate  (Thermodynamic)
% q_r = correct base pair dissociation rate (Thermodynamic)
% q_w= incorrect base pair dissociation rate (Thermodynamic) 
% a_r= right modulating factor for correct base  pair (Kinetic)
% b_r =left modulating factor for correct base pair (Kinetic)
% a_w = right modulating factor for incorrect base pair (Kinetic)
% b_w= left modulating factor for incorrect base pair (Kinetic)
% w = covlent bond formation rate

%%% ==== Output ====%%%

% Q =transition rate matrix

%%% ==== before absrobing final state index ===%%%

before_abs_final_state_idx=find(all(s~=0 & s~=5 & s~=6,2));


% ======= Absorbing state index========%
abs_state_idx=find(all(s~=0 & s~=1 & s~=2,2));


% ====== Making the rate matrix absorbing ========%

Q=zeros(size(s,1),size(s,1));

%%%  thermodynamic change
for i=1:size(s,1)
    for j=1:size(s,1)
        c=s(i,:)~=s(j,:);
        c_sum=sum(c(:,:));
        if c_sum==1
            k=find(c);
            if (s(i,k)==0 && s(j,k)==1) || (s(i,k)==0 && s(j,k)==2)
                Q(i,j)=p;
            end
            if (s(i,k)==1 && s(j,k)==0)
                Q(i,j)=q_r;
            end
            if (s(i,k)==2 && s(j,k)==0)
                Q(i,j)=q_w;
            end
        end
    end
end


%%% kinetic change
for i=1:size(s,1)
    for j=1:size(s,1)
        c=s(i,:)~=s(j,:);
        c_sum=sum(c(:,:));
        if c_sum==1
            k=find(c);
            if (s(i,k)==0 && s(j,k)==1) || (s(i,k)==1 && s(j,k)==0) ...
                    || (s(i,k)==0 && s(j,k)==2) || (s(i,k)==2 && s(j,k)==0)
                if k==1
                    if s(i,2)==1
                        Q(i,j)=Q(i,j)*b_r;
                    end
                     if s(i,2)==2
                        Q(i,j)=Q(i,j)*b_w;
                     end
                elseif k==size(s,2)
                    if s(i,size(s,2)-1)==1
                        Q(i,j)=Q(i,j)*a_r;
                    end
                     if s(i,size(s,2)-1)==2
                        Q(i,j)=Q(i,j)*a_w;
                     end
                else
                    if s(i,k-1)==1
                        Q(i,j)=Q(i,j)*a_r;
                    end
                    if s(i,k+1)==1
                        Q(i,j)=Q(i,j)*b_r;
                    end
                    if s(i,k-1)==2
                        Q(i,j)=Q(i,j)*a_w;
                    end
                    if s(i,k+1)==2
                        Q(i,j)=Q(i,j)*b_w;
                    end
                end
            end
        end
    end
end

Q(before_abs_final_state_idx(1),abs_state_idx(1))=w;
Q(before_abs_final_state_idx(2),abs_state_idx(2))=w;
Q(abs_state_idx(1),before_abs_final_state_idx(1))=10^2;
Q(abs_state_idx(2),before_abs_final_state_idx(2))=10^2;

%%% diagonal element
Q=Q-diag(sum(Q,2));


end





