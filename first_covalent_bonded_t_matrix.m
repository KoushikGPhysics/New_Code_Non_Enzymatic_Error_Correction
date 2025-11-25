function[Q]=first_covalent_bonded_t_matrix(s,p,q_r,q_w,a_r,b_r,a_w,b_w,w)


N=size(s,1)-2;

Q=zeros(size(s,1),size(s,1));

for i=1:N
    for j=1:N
        state_i=s(i,1:5);
        state_j=s(j,1:5);
        c=state_i~=state_j;
         c_sum=sum(c(:,:));
         if c_sum==1
              k=find(c);
              if state_i(k)==0 && state_j(k)==1 || state_i(k)==0 && state_j(k)==2
                    Q(i,j)=p;
               end
               if state_i(k)==1 && state_j(k)==0
                    Q(i,j)=q_r;
               end
               if state_i(k)==2 && state_j(k)==0
                    Q(i,j)=q_w;
               end
         end
    end
end


for i=1:N
    for j=1:N
         c=s(i,:)~=s(j,:);
        c_sum=sum(c(:,:));
        if c_sum==1
            k=find(c);
            if k==6
               if all(s(i,1:3)~=0) && s(i,6)==0
                   Q(i,j)=w;
               end
            end
        end
    end
end

for i=1:N
    for j=1:N
        state_i=s(i,1:5);
        state_j=s(j,1:5);
        c=state_i~=state_j;
        c_sum=sum(c(:,:));
        if c_sum==1
            k=find(c);
            if (state_i(k)==0 && state_j(k)==1) || (state_i(k)==1 && state_j(k)==0) ...
                    || (state_i(k)==0 && state_j(k)==2) || (state_i(k)==2 && state_j(k)==0)
                if k==1
                    if state_i(2)==1
                        Q(i,j)=Q(i,j)*b_r;
                    end
                     if state_i(2)==2
                        Q(i,j)=Q(i,j)*b_w;
                     end
                elseif k==5
                    if s(i,5-1)==1
                        Q(i,j)=Q(i,j)*a_r;
                    end
                     if s(i,5-1)==2
                        Q(i,j)=Q(i,j)*a_w;
                     end
                else
                    if state_i(k-1)==1
                        Q(i,j)=Q(i,j)*a_r;
                    end
                    if state_i(k+1)==1
                        Q(i,j)=Q(i,j)*b_r;
                    end
                    if state_i(k-1)==2
                        Q(i,j)=Q(i,j)*a_w;
                    end
                    if state_i(k+1)==2
                        Q(i,j)=Q(i,j)*b_w;
                    end
                end
            end
        end
    end
end


% === MANUAL STATE LINKS ===%

Q(40,49) = w;
Q(48,50) = w;

% === DIAGONAL CORRECTION ===%
Q = Q - diag(sum(Q,2));

end


