clear;

load('kron_wgn.mat');
% load('kron_impulse.mat');
% load('twf_metered.mat');

load('metered.mat');

buses = 1:N; 
 
freq_lin = freq_lin(:,1:50000);

metered=randsample(69,60);

non_metered = setdiff(buses,metered);
noise=randn(size(freq_lin(metered,:)));

lambda2=lambda2(2:6);
U2=U2(:,2:6);
V2=V2(:,2:6);

%% Preprocessing (filtering and downsampling):

k_downsample=20;
noise_var= (0.01)^2;

d1  = designfilt('bandpassiir','FilterOrder',8, ...
         'HalfPowerFrequency1',0.2,'HalfPowerFrequency2',0.8, ...
         'SampleRate',1000);

% d1 = designfilt('lowpassiir','FilterOrder',8, ...
% 'HalfPowerFrequency',2,'DesignMethod','butter','SampleRate',1000);

freq_filtered = filtfilt(d1,freq_lin');
freq_filtered = freq_filtered';

% fvtool(d1)

%%
freq_test  = freq_lin(metered,:) + sqrt(noise_var)*noise;

scale=100;
freq_test = filtfilt(d1,freq_test');
freq_test = scale*freq_test';
noise_var = noise_var*(scale^2);
Sq=sqrt((gama^2)-4*lambda2);
% 
a = 0.5 - gama./Sq/2;
b = 0.5 + gama./Sq/2;

% a = 1./Sq;
% b = -a;

c = (-gama + Sq)/2;
d = (-gama - Sq)/2; 

T_frame = 20000:28000;
freq_lin   = freq_lin(:,T_frame);
freq_test  = freq_test(:,T_frame);
freq_filtered  = freq_filtered(:,T_frame);

t = t(T_frame);
t = t(1:k_downsample:end);
T = length(t);

freq_test = freq_test(:,1:k_downsample:end);

[A,a_ij,b_ij]=finding_corr_constants(freq_test,U2(metered,:),a,b,c,d,noise_var,lambda2);
 
a_ji = transpose(a_ij);
b_ji = transpose(b_ij);

%% Construction of eigenstate covariances:
Delta = toeplitz(t-t(1));

K_ij = cell(size(U2,2),1);
K_ji = cell(size(U2,2),1);
K = cell(size(U2,2),1);

for i=1:1:size(U2,2)
    for j=1:1:size(U2,2)
   
        K_ij{i,j} = A(i,j)*(a_ij(i,j)*exp(c(i)*Delta) + b_ij(i,j)*exp(d(i)*Delta));
        K_ij{i,j} = tril(K_ij{i,j});
        K_ji{i,j} = A(i,j)*(a_ji(i,j)*exp(c(j)*Delta) + b_ji(i,j)*exp(d(j)*Delta));
        K_ji{i,j} = triu(K_ji{i,j});
        K{i,j} = K_ij{i,j}+K_ji{i,j}-diag(diag(K_ji{i,j}));
        
    end
end
K_bd = cell2mat(K);

%% System state covariances

IT  = eye(T);

US  = kron(U2(metered,:),IT);
Un  = kron(U2(non_metered,:),IT);

E_omega  = US * K_bd * US' + noise_var*eye(size(US,1));
E_n      = Un * K_bd * US';
E_nn     = Un * K_bd * Un';

%% Bayesian estimation:

E_no = E_n/E_omega;
me_omega = E_no*reshape(freq_test',[],1);

cov  = E_nn - E_no*E_n';
sigma = sqrt(abs(diag(cov)));

me_omega = reshape(me_omega,length(t),length(non_metered));
me_omega = me_omega'/scale;

sigma = reshape(sigma,length(t),length(non_metered));
sigma = sigma';

sigmal = me_omega - sigma;
sigmau = me_omega + sigma;

clear E_n E_theta_omega E_omega_inv E_nn 

%%
n1=1;
str_n = num2str(non_metered(n1)); 
figure;

hold on;box on;grid on;

shade(t,sigmal(n1,:)',t,sigmau(n1,:)','FillType',[1 2;2 1],'FillColor','b');
plot(t,me_omega(n1,:),'-ob');
plot(t,freq_filtered(non_metered(n1),1:k_downsample:end),'k');

error= abs(me_omega - freq_filtered(non_metered,1:k_downsample:length(freq_lin)));
mean(mean(error))
max(max(error))
