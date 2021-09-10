function [V,L,m,gama,d1,freq,theta,input,lambda,sys_d,A,C]=data_generation_nonreduced(cases,input,tfinal,t_sample,kron_reduction_flag,homogeneous_flag)
% Inputs:
% cases: benchmark case file
% input: power disturbance on generators (N*T: N:# buses,T:#time instances)
% t_final: duration of the oscillations in seconds
% t_sample: sampling intervals in seconds
% kron_reduction_flag = {0:no kron-reduction,1:kron reduction}
% homogeneous_flag = {1:homo,0:non-homo}


% Outputs:
% theta:rotor angles(rad)(N*T)% m: an N-sized vector containing the generators' moments of inertia
% d1: an N-sized vector containing the damping coefficient
% L: the kron reduced Laplacian matrix (N*N)
% V: eigenvectors of the scaled Laplacian (N*N)
% gama: homogenity coefficient (d1=gama*m)

trange  = 0:t_sample:tfinal;
Nsample = length(trange);
mpc     =loadcase(cases);

included = mpc.gen(:,1);

P_load = mpc.bus(:,3)/mpc.baseMVA;
Q_load = mpc.bus(:,4)/mpc.baseMVA;

excluded=setdiff(1:length(P_load),included);

if (kron_reduction_flag==1)
    %% Kron Reduction:
    input = input(included,:);
    
    m=mpc.M(:,2);
    N=length(included);
    
    Y = full(makeYbus(mpc));
    
    J = full(makeJac(mpc,'fullJac'));
    J = J(1:length(P_load),1:length(P_load));
    
%     Z = inv(Y);
%     Z_11=Z(included,included);
%     Z_22=Z(excluded,excluded);
%     Z_12=Z(included,excluded);
%     Z = Z_11-Z_12*inv(Z_22)*transpose(Z_12);
%     
    yg=-1i*1./(mpc.Ls(1:N,2));
    
    
    L = J(included,included)-J(included,excluded)*inv(J(excluded,excluded))*J(excluded,included);
    
    Y=Y-diag(P_load-Q_load*1i);
    
    Y_GL=Y(included,excluded);
    Y_LL=Y(excluded,excluded);
    Y_G_tilda=Y(included,included)+diag(yg);
    
    YA=Y_G_tilda-Y_GL*(eye(size(Y_LL,1))/Y_LL)*Y_GL.';
    YA=(eye(length(YA))/YA)*diag(yg);
    Y=conj(-diag(yg)*YA);
    
    BB=real(Y);
    
    L2=zeros(size(BB));
    for i=1:size(BB,1)
        for j=1:size(BB,2)
            if i==j
                L2(i,j)=-sum(BB(i,:))+BB(i,j);
            else
                L2(i,j)=BB(i,j);
            end
        end
    end
L=L2;

elseif (kron_reduction_flag==0)
    
    m           = 0.01*ones(length(P_load),1);
    m(included) = mpc.M(:,2);

    N           = length(m);
    
    Y=full(makeYbus(mpc));
    Z = inv(Y);
    L=-imag(Y);
    
end
%%


if (homogeneous_flag==1)
    gama=0.0057/mean(m);
    d1=gama*m;
elseif (homogeneous_flag==0)
    d1=0.0057*ones(N,1);
    gama = mean(d1)/mean(m);
end



LM=diag(1./sqrt(m))*L*diag(1./sqrt(m));

[V,lambda]=eig(LM);


A = [zeros(N,N)  eye(N);
    -(diag(m)\L)    -diag(m.\d1)];

B = [zeros(N,N);
    diag(1./m)];

C = [zeros(N,N) eye(N)];

D = zeros(N,N);

sys = ss(A,B,C,D);

sys_d = c2d(sys,t_sample);

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

ic = zeros(2*N,1);

freq = zeros(N,Nsample-1);

x = zeros(2*N,Nsample);
x(:,1) = ic;

for cnt = 1:Nsample-1
    x(:,cnt+1)  = Ad*x(:,cnt) + Bd*input(:,cnt);
    freq(:,cnt) = Cd*x(:,cnt) + Dd*input(:,cnt);
end

theta=x(1:N,:);



end