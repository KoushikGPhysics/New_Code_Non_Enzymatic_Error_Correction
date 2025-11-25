
%%% Function for transition rate matrix of absorbing Markov chain

%%% ==== Input =====%%%

% s= state matrix
% p=base pair formation rate  (Thermodynamic)
% q_r = correct base pair dissociation rate (Thermodynamic)
% q_w= incorrect base pair dissociation rate (Thermodynamic) 
% a_r= right modulating factor for correct base  pair (Kinetic)
% b_r =left modulating factor for correct base pair (Kinetic)
% a_w = right modulating factor for incorrect base pair (Kinetic)
% b_w= left modulating factor for incorrect base pair (Kinetic)


%%% ==== Output ====%%%

% Q =transition rate matrix

function[Q]=absorbing_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w)

% ======= Final state index========%
final_state_idx=all(s~=0,2);

Q=zeros(size(s,1),size(s,1));

%%%   thermodynamic change
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

%%%  kinetic change
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

%%% Absobing state element
Q(final_state_idx,:)=0;


%%% diagonal element
Q=Q-diag(sum(Q,2));

end


