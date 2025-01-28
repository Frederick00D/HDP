% Sample size
n = 1000;
% Repitation times
S = 5000;

% contamination_percentage_list
contamination_percentage_list = [0,0.05,0.1,0.2,0.3];

% Real distribution, Normal
mu = 5;
sigma = 2;

% Bandwidth
h = 0.448;
% X = normrnd(mu,sigma,n,1);
% h = density_kernel.silverman(X);

% Sensitivity convergence rate, p_n \in (1,2)
p_n = 1.7;

% learning rate
eta = 0.5;

% Start point
theta_start = [1,1];

%% Generate data and resample
data = normrnd(mu,sigma,n,S);
K = @(x) OPT.epanechnikov(x);
m_n = 5000;
data_resample = zeros(m_n,S);
parfor i=1:S
    X = data(:,i);
    data_resample(:,i) = OPT.rnd_f_k_epa(m_n, X, h);
end

%% Generate contaminated data
data_contain = zeros(n,S,length(contamination_percentage_list));
m_n = 5000;
data_contain_resample = zeros(m_n,S,length(contamination_percentage_list));
q_985 = norminv(0.985, mu, sigma);
q_995 = norminv(0.995, mu, sigma);    
for j = 1:length(contamination_percentage_list)
    cp = contamination_percentage_list(j);
    num_to_replace = round(cp * n);
    uniform_values = (q_995 - q_985) * rand(num_to_replace, 1) + q_985;
    replace_indices = randperm(n, num_to_replace);
    data_contain(:,:,j) = data;
    for i=1:S
        data_contain(replace_indices,i,j) = uniform_values;
    end
    parfor i=1:S
        X = data_contain(:,i,j);
        data_contain_resample(:,i,j) = OPT.rnd_f_k_epa(m_n, X, h);
    end
end

%% Save data
dataname = append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat');
save(dataname);

%% Calculate PMHDE
clearvars -except n S;

% HDP data PNR
load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
ep_list = [2,0.6,0.2];
iteration_number = 30; % Algorithm maximum iteration number
noise_K = 5; % Iteration number for noise
method = 2; % Method, method=1 for Gradient descent, method=2 for Newton-Raphson
[collection,collection_CI,collection_CI_u,TRAJ] = simulation.collection_normal(data,data_resample,method,S,h,p_n,noise_K,iteration_number,ep_list,eta,theta_start);
name_data = append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat');
save(name_data, 'collection', 'collection_CI', 'collection_CI_u', 'TRAJ','ep_list','iteration_number','noise_K','method');
clearvars -except n S 

% HDP data PGD 
load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
ep_list = [2,0.6,0.2];
iteration_number = 80;
noise_K = 50;
method = 1;
[collection,collection_CI,collection_CI_u,TRAJ] = simulation.collection_normal(data,data_resample,method,S,h,p_n,noise_K,iteration_number,ep_list,eta,theta_start);
name_data = append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat');
save(name_data, 'collection', 'collection_CI', 'collection_CI_u', 'TRAJ','ep_list','iteration_number','noise_K','method');
clearvars -except n S 

% DPD data PNR
load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
lam = 1;
ep_list = [2, 0.6*2, 0.2*2];
iteration_number = 30; % Algorithm maximum iteration number
noise_K = 5; % Iteration number for noise
method = 2; % Method, method=1 for Gradient descent, method=2 for Newton-Raphson
[collection,collection_CI,collection_CI_u,TRAJ] = simulation_pdp.collection_normal(data,data_resample,method,S,h,p_n,noise_K,iteration_number,lam,ep_list,eta,theta_start);
name_data = append('PDP_n',num2str(n), ' ', 'lam_', num2str(lam), ' ', 'Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat');
save(name_data, 'collection', 'collection_CI', 'collection_CI_u', 'TRAJ','ep_list','iteration_number','noise_K','method');
clearvars -except n S  

% DPD data PGD
load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
lam = 1;
ep_list = [2, 0.6*2, 0.2*2];
iteration_number = 80;
noise_K = 50;
method = 1;
[collection,collection_CI,collection_CI_u,TRAJ] = simulation_pdp.collection_normal(data,data_resample,method,S,h,p_n,noise_K,iteration_number,lam,ep_list,eta,theta_start);
name_data = append('PDP_n',num2str(n), ' ', 'lam_', num2str(lam), ' ', 'Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat');
save(name_data, 'collection', 'collection_CI', 'collection_CI_u', 'TRAJ','ep_list','iteration_number','noise_K','method');
clearvars -except n S 

% HDP data contamination PNR 
load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
ep_list = [2,0.6,0.2];
iteration_number = 30; % Algorithm maximum iteration number
noise_K = 5; % Iteration number for noise
method = 2;
[collection,collection_CI,collection_CI_u] = simulation.collection_normal_contaminatin(data_contain,data_contain_resample,method,S,h,p_n,noise_K,iteration_number,ep_list,eta,theta_start, contamination_percentage_list);
name_data = append('Contamination',' ','Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat');
save(name_data, 'collection', 'collection_CI', 'collection_CI_u','ep_list','iteration_number','noise_K','method');
clearvars -except n S 

% HDP data contamination PGD
load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
ep_list = [2,0.6,0.2];
iteration_number = 80;
noise_K = 50;
method = 1;
[collection,collection_CI,collection_CI_u] = simulation.collection_normal_contaminatin(data_contain,data_contain_resample,method,S,h,p_n,noise_K,iteration_number,ep_list,eta,theta_start, contamination_percentage_list);
name_data = append('Contamination',' ','Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat');
save(name_data, 'collection', 'collection_CI', 'collection_CI_u','ep_list','iteration_number','noise_K','method');
clearvars -except n S 






