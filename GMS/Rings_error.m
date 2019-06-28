%% Define error function

function lsq=Rings_error(p,Qexp,I_av,ind)
alpha=p(1);
qmax=p(2);

[Q_111,Q_1m1m1,Q_m11m1,Q_m1m11]= Model_Modulated_rings(alpha,qmax);
%Potential function along the rings
U_Q=@(Q,n) (Q(:,1)*n(1)+Q(:,2)*n(2)+Q(:,3)*n(3)).^2+alpha*(Q(:,1).^4+Q(:,2).^4+Q(:,3).^4);
U_111=U_Q(Q_111,[1,1,1]);
U_1m1m1=U_Q(Q_1m1m1,[1,-1,-1]);
U_m1m11=U_Q(Q_m1m11,[-1,-1,1]);
U_m11m1=U_Q(Q_m11m1,[-1,1,-1]);


%plot3(Q_111(:,1),Q_111(:,2),Q_111(:,3),'linewidth',2)
%plot3(Q_1m1m1(:,1),Q_1m1m1(:,2),Q_1m1m1(:,3),'linewidth',2)
%plot3(Q_m11m1(:,1),Q_m11m1(:,2),Q_m11m1(:,3),'linewidth',2)
%plot3(Q_m1m11(:,1),Q_m1m11(:,2),Q_m1m11(:,3),'linewidth',2)

Qmod=[Q_111;Q_1m1m1;Q_m11m1;Q_m1m11]; %model Q vectors organized as [Q_111(:,1),Q_111(:,2),Q_111(:,2); 
                                                                    %[Q_1m1m1(:,1),Q_1m1m1(:,2),Q_1m1m1(:,2)]...
Umod=[U_111;U_1m1m1;U_m1m11;U_m11m1]; %model potentials

lsq=0;
for i=1:size(Qexp,1)
    lsq=lsq+I_av(ind(i))*min(sum((bsxfun(@minus,Qexp(i,:),Qmod)).^2,2)); %Find q-vectors Qexp with minimal distance to model q-vectors, Qmod. Distance is weighed by the SANS intensity (I_av)
end

end