

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

% % % % % Q(abs_state_idx(1),before_abs_final_state_idx(1))=w;
% % % % % Q(abs_state_idx(2),before_abs_final_state_idx(2))=w;
% % % % % 

%%% diagonal element
Q=Q-diag(sum(Q,2));


end


