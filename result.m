% Sample size
n = 1000;
% Repitation times
S = 5000;

mu = 5;
sigma = 2;
contamination_percentage_list = [0,0.05,0.1,0.2,0.3];

%% Table for HDP (table 1-2)

% Method = 1 for PGD and Method = 2 for PNR
method = 1;
noise_K = 50;

% load(append('Data_', 'S_',num2str(S),' ', 'n_', num2str(n),'.mat'));
load(append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));

% Use this one for n=1000
simulation.table_normal(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

% % Use this one for n=200,300,500
% simulation.table_normal_threshold(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

method = 2;
noise_K = 5;
load(append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));

% Use this one for n=1000
simulation.table_normal(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

% % Use this one for n=200,300,500
% simulation.table_normal_threshold(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

%% Table for PDP (table 3-4)
lam = 1;

method = 1;
noise_K = 50;
load(append('PDP_n',num2str(n), ' ', 'lam_', num2str(lam), ' ', 'Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));

% Use this one for n=1000
simulation.table_normal(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

% % Use this one for n=200,300,500
% simulation.table_normal_threshold(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

method = 2;
noise_K = 5;
load(append('PDP_n',num2str(n), ' ', 'lam_', num2str(lam), ' ', 'Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));

% Use this one for n=1000
simulation.table_normal(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

% % Use this one for n=200,300,500
% simulation.table_normal_threshold(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)

%% Table for contamination (table 5-6)

method = 1;
noise_K = 50;
load(append('Contamination',' ','Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));

% Use this one for n=1000
simulation.table_normal_contamination(collection,collection_CI,collection_CI_u,n,ep_list,contamination_percentage_list,mu,sigma)

% % Use this one for n=200,300,500
% simulation.table_normal_contamination_threshold(collection,collection_CI,collection_CI_u,n,ep_list,contamination_percentage_list,mu,sigma)

method = 2;
noise_K = 5;
load(append('Contamination',' ','Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));

% Use this one for n=1000
simulation.table_normal_contamination(collection,collection_CI,collection_CI_u,n,ep_list,contamination_percentage_list,mu,sigma)

% % Use this one for n=200,300,500
% simulation.table_normal_contamination_threshold(collection,collection_CI,collection_CI_u,n,ep_list,contamination_percentage_list,mu,sigma)


%% Plots for PMHDE (Figure 1-4)

method = 1;
noise_K = 50;
load(append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));
S_short = 20;
simulation.plot_normal(collection,mu,sigma,noise_K,iteration_number,ep_list,S_short)
simulation.plot_normal_traj(TRAJ,noise_K,iteration_number,ep_list,S_short)

method = 2;
noise_K = 5;
load(append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', num2str(n),' ','K_', num2str(noise_K),'.mat'));
S_short = 20;
simulation.plot_normal(collection,mu,sigma,noise_K,iteration_number,ep_list,S_short)
simulation.plot_normal_traj(TRAJ,noise_K,iteration_number,ep_list,S_short)

%% Plot CI coverage
S = 5000;
method = 1;
noise_K = 50;
mu = 5;
sigma = 2;
n_list = [50,100,200,300,500,750,1000];
% This use the first 2 elements in the ep_list, ep_list(1) is non-private
CI_coverage = ones(length(n_list),3);
for i = 1:length(n_list)
    n = n_list(i);
    name_data = append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', ...
    num2str(n),' ','K_', num2str(noise_K),'.mat');
    load(name_data);

    % non-private coverage
    CI_coverage(i,1) = round(mean(squeeze(collection_CI(noise_K,1,1,:,1)) < mu & squeeze(collection_CI(noise_K,1,2,:,1)) > mu), 3);

    % private uncorrected coverage
    CI_coverage(i,2) = round(mean(squeeze(collection_CI_u(noise_K,1,1,:,2)) < mu & squeeze(collection_CI_u(noise_K,1,2,:,2)) > mu), 3);

    % private corrected coverage
    CI_coverage(i,3) = round(mean(squeeze(collection_CI(noise_K,1,1,:,2)) < mu & squeeze(collection_CI(noise_K,1,2,:,2)) > mu), 3);
end

% Prepare the plot
figure;
hold on;

% Define markers and colors for each method
markers = {'o', 's', 'd'}; % Circle, square, diamond
colors = {'r', 'g', 'b'};  % Red, green, blue
method_names = {'Non-private coverage', 'Private uncorrected coverage', 'Private corrected coverage'};

% Plot each method
for j = 1:3
    x = n_list;    % Apply offset to x-values
    y = CI_coverage(:, j);      % Corresponding coverage rates
    scatter(x, y, 50, colors{j}, markers{j}, 'filled', 'DisplayName', method_names{j});
end

% Add horizontal line at y = 0.95
yline(0.95, 'k--', 'LineWidth', 1,'HandleVisibility', 'off');

% Customize the plot
xlabel('Sample size');
ylabel('Coverage rate');
title('Coverage rate vs sample size');

ylim([0, 1]);  % Set y-axis range from 0 to 1
% Calculate x-axis limits with padding
x_padding = 20;  % Adjust this value to increase or decrease padding
x_min = min(n_list)  - x_padding;
x_max = max(n_list)  + x_padding;
xlim([x_min, x_max]);  % Set x-axis limits

legend('Location', 'best');
grid on;
hold off;


method = 2;
noise_K = 5;
mu = 5;
sigma = 2;
n_list = [50,100,200,300,500,750,1000];
% This use the first 2 elements in the ep_list, ep_list(1) is non-private
CI_coverage = ones(length(n_list),3);
for i = 1:length(n_list)
    n = n_list(i);
    name_data = append('Method_',num2str(method), ' ' ,'S_',num2str(S),' ', 'n_', ...
    num2str(n),' ','K_', num2str(noise_K),'.mat');
    load(name_data);

    % non-private coverage
    CI_coverage(i,1) = round(mean(squeeze(collection_CI(noise_K,1,1,:,1)) < mu & squeeze(collection_CI(noise_K,1,2,:,1)) > mu), 3);

    % private uncorrected coverage
    CI_coverage(i,2) = round(mean(squeeze(collection_CI_u(noise_K,1,1,:,2)) < mu & squeeze(collection_CI_u(noise_K,1,2,:,2)) > mu), 3);

    % private corrected coverage
    CI_coverage(i,3) = round(mean(squeeze(collection_CI(noise_K,1,1,:,2)) < mu & squeeze(collection_CI(noise_K,1,2,:,2)) > mu), 3);
end

% Prepare the plot
figure;
hold on;

% Define markers and colors for each method
markers = {'o', 's', 'd'}; % Circle, square, diamond
colors = {'r', 'g', 'b'};  % Red, green, blue
method_names = {'Non-private coverage', 'Private uncorrected coverage', 'Private corrected coverage'};

% Plot each method
for j = 1:3
    x = n_list;    % Apply offset to x-values
    y = CI_coverage(:, j);      % Corresponding coverage rates
    scatter(x, y, 50, colors{j}, markers{j}, 'filled', 'DisplayName', method_names{j});
end

% Add horizontal line at y = 0.95
yline(0.95, 'k--', 'LineWidth', 1,'HandleVisibility', 'off');

% Customize the plot
xlabel('Sample size');
ylabel('Coverage rate');
title('Coverage rate vs sample size');

ylim([0, 1]);  % Set y-axis range from 0 to 1
% Calculate x-axis limits with padding
x_padding = 20;  % Adjust this value to increase or decrease padding
x_min = min(n_list)  - x_padding;
x_max = max(n_list)  + x_padding;
xlim([x_min, x_max]);  % Set x-axis limits

legend('Location', 'best');
grid on;
hold off;


