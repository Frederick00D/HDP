classdef simulation
    methods(Static)

% Input: 
% Sample size: n
% Iteration number for noise: noise_K
% Repitation times: S
% Bandwidth: h
% Algorithm maximum iteration number: iteration_number
% Privacy level: ep
% Learning rate: eta
% Start point: theta_start
% Sensitivity convergence rate: p_n

function table_normal(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)
    T = ones(10,length(ep_list));
    RowNames = {'Sample size';'K';'mu: Estimator';'mu: std';'mu: Coverage'; 'mu: Uncorrected Coverage';
        'sigma: Estimator';'sigma: std';'sigma: Coverage'; 'sigma: Uncorrected Coverage'};
    VariableNames = cell(1,length(ep_list));
    for k=1:length(ep_list)
        ep = ep_list(k);
        VariableNames{k} = sprintf('epsilon = %.2f', ep);
        T(1,k) = n;
        T(2,k) = noise_K;

        data = collection(noise_K,1,:,k);
        data_short = data;
        
        estimator_mu = mean(data_short);
        estimator_mu_std = std(data_short);

        data = collection(noise_K,2,:,k);
        data_short = data;

        estimator_sigma = mean(data_short);
        estimator_sigma_std = std(data_short);

        T(3,k) = round(estimator_mu, 3);
        T(4,k) = round(estimator_mu_std, 3);
        T(5,k) = round(mean(collection_CI(noise_K,1,1,:,k)<mu & collection_CI(noise_K,1,2,:,k)>mu), 3); 
        T(6,k) = round(mean(collection_CI_u(noise_K,1,1,:,k)<mu & collection_CI_u(noise_K,1,2,:,k)>mu), 3); 
        
        T(7,k) = round(estimator_sigma, 3);
        T(8,k) = round(estimator_sigma_std, 3);
        T(9,k) = round(mean(collection_CI(noise_K,2,1,:,k)<sigma & collection_CI(noise_K,2,2,:,k)>sigma), 3);
        T(10,k) = round(mean(collection_CI_u(noise_K,2,1,:,k)<sigma & collection_CI_u(noise_K,2,2,:,k)>sigma), 3);
    end
    T = array2table(T,'RowNames',RowNames,'VariableNames',VariableNames);
    f = uifigure;
    uitable(f, 'Data', T);
    disp(T);
end

function table_normal_threshold(collection,collection_CI,collection_CI_u,mu,sigma,n,noise_K,ep_list)
    T = ones(10,length(ep_list));
    RowNames = {'Sample size';'K';'mu: Estimator';'mu: std';'mu: Coverage'; 'mu: Uncorrected Coverage';
        'sigma: Estimator';'sigma: std';'sigma: Coverage'; 'sigma: Uncorrected Coverage'};
    VariableNames = cell(1,length(ep_list));
    for k=1:length(ep_list)
        ep = ep_list(k);
        VariableNames{k} = sprintf('epsilon = %.2f', ep);
        T(1,k) = n;
        T(2,k) = noise_K;

        if k==1
            data_mu = collection(noise_K,1,:,k);
            data_sigma = collection(noise_K,2,:,k);
            data_short_mu = data_mu;
            data_short_sigma = data_sigma;
        else
            data_mu = collection(noise_K,1,:,k);
            data_sigma = collection(noise_K,2,:,k);
            mu_hat  = mean(collection(noise_K,1,:,1));
            sigma_hat = mean(collection(noise_K,2,:,1));
            q1 = norminv(0.007,mu_hat,sigma_hat);
            q9 = norminv(0.995,mu_hat,sigma_hat);
            data_short_mu = data_mu(data_mu >= q1 & data_mu <= q9);
            data_short_sigma = data_sigma(data_mu >= q1 & data_mu <= q9);
        end
        
        estimator_mu = mean(data_short_mu);
        estimator_mu_std = std(data_short_mu);

        estimator_sigma = mean(data_short_sigma);
        estimator_sigma_std = std(data_short_sigma);

        T(3,k) = round(estimator_mu, 3);
        T(4,k) = round(estimator_mu_std, 3);
        T(5,k) = round(mean(collection_CI(noise_K,1,1,:,k)<mu & collection_CI(noise_K,1,2,:,k)>mu), 3); 
        T(6,k) = round(mean(collection_CI_u(noise_K,1,1,:,k)<mu & collection_CI_u(noise_K,1,2,:,k)>mu), 3); 
        
        T(7,k) = round(estimator_sigma, 3);
        T(8,k) = round(estimator_sigma_std, 3);
        T(9,k) = round(mean(collection_CI(noise_K,2,1,:,k)<sigma & collection_CI(noise_K,2,2,:,k)>sigma), 3);
        T(10,k) = round(mean(collection_CI_u(noise_K,2,1,:,k)<sigma & collection_CI_u(noise_K,2,2,:,k)>sigma), 3);
    end
    T = array2table(T,'RowNames',RowNames,'VariableNames',VariableNames);
    f = uifigure;
    uitable(f, 'Data', T);
    disp(T);
end

function table_normal_contamination(collection, collection_CI, collection_CI_u, n, ep_list, contamination_percentage_list, mu, sigma)
    numEps = length(ep_list);
    numRows = 8 * numEps + 9; % Adjusted for all rows
    numCols = length(contamination_percentage_list);
    T = ones(numRows, numCols);
    RowNames = cell(numRows, 1); % Cell array for row names

    % Initial Row Names
    RowNames{1} = 'Sample size';
    RowNames{2} = 'MLE_mu';
    RowNames{3} = 'MLE_mu 95 percent coverage';
    RowNames{4} = 'MLE_mu uncorrected 95 percent coverage';
    RowNames{5} = 'MLE_sigma';
    RowNames{6} = 'MLE_sigma 95 percent coverage';
    RowNames{7} = 'MLE_sigma uncorrected 95 percent coverage';
    RowNames{8} = 'Std of MLE_mu';
    RowNames{9} = 'Std of MLE_sigma';

    VariableNames = cell(1, numCols); % Cell array for variable (column) names
    for j = 1:numCols
        VariableNames{j} = sprintf('Contamination = %.2f', contamination_percentage_list(j));
    end

    % Populate initial rows
    T(1, :) = n;

    % Compute MLE_mu and MLE_sigma with trimming and standard deviations
    data_mle = squeeze(collection(:,:,end,:));    
    for j = 1:numCols
        % Data for MLE_mu and MLE_sigma
        data_mu = squeeze(data_mle(1,:,j));
        data_sigma = squeeze(data_mle(2,:,j));

        data_mu_trimmed = data_mu;

        T(2, j) = round(mean(data_mu_trimmed), 3); % Mean of trimmed MLE_mu
        T(8, j) = round(std(data_mu_trimmed), 3);  % Std of trimmed MLE_mu

        data_sigma_trimmed = data_sigma;

        T(5, j) = round(mean(data_sigma_trimmed), 3); % Mean of trimmed MLE_sigma
        T(9, j) = round(std(data_sigma_trimmed), 3);  % Std of trimmed MLE_sigma

        % Corrected Coverage
        data_mle_CI = squeeze(collection_CI(:,:,:,end,j));
        T(3, j) = round(mean(data_mle_CI(1,1,:) < mu & data_mle_CI(1,2,:) > mu), 3);
        T(6, j) = round(mean(data_mle_CI(2,1,:) < sigma & data_mle_CI(2,2,:) > sigma), 3);

        % Uncorrected Coverage
        data_mle_CI_u = squeeze(collection_CI_u(:,:,:,end,j));
        T(4, j) = round(mean(data_mle_CI_u(1,1,:) < mu & data_mle_CI_u(1,2,:) > mu), 3);
        T(7, j) = round(mean(data_mle_CI_u(2,1,:) < sigma & data_mle_CI_u(2,2,:) > sigma), 3);
    end

    % Compute for each epsilon in ep_list
    for k = 1:numEps
        ep = ep_list(k);
        data = squeeze(collection(:,:,k,:));
        for j = 1:numCols
            % Data for mu and sigma
            data_mu = squeeze(data(1,:,j));
            data_sigma = squeeze(data(2,:,j));
 
            data_mu_trimmed = data_mu;

            T(9 + 8*k-6, j) = round(mean(data_mu_trimmed), 3); % Mean of trimmed mu
            T(9 + 8*k-5, j) = round(std(data_mu_trimmed), 3);  % Std of mu
 
            data_sigma_trimmed = data_sigma;
           
            T(9 + 8*k-3, j) = round(mean(data_sigma_trimmed), 3); % Mean of trimmed sigma
            T(9 + 8*k-2, j) = round(std(data_sigma_trimmed), 3);  % Std of sigma

            % Coverage (corrected and uncorrected)
            data_CI = squeeze(collection_CI(:,:,:,k,j)); 
            data_CI_u = squeeze(collection_CI_u(:,:,:,k,j));
            T(9 + 8*k-4, j) = round(mean(data_CI(1,1,:) < mu & data_CI(1,2,:) > mu), 3); % Corrected mu coverage
            T(9 + 8*k-1, j) = round(mean(data_CI(2,1,:) < sigma & data_CI(2,2,:) > sigma), 3); % Corrected sigma coverage
            T(9 + 8*k-7, j) = round(mean(data_CI_u(1,1,:) < mu & data_CI_u(1,2,:) > mu), 3); % Uncorrected mu coverage
            T(9 + 8*k-0, j) = round(mean(data_CI_u(2,1,:) < sigma & data_CI_u(2,2,:) > sigma), 3); % Uncorrected sigma coverage
        end

        % Update RowNames
        RowNames{9 + 8*k-7} = sprintf('ep = %.2f : mu uncorrected 95 percent coverage', ep);
        RowNames{9 + 8*k-6} = sprintf('ep = %.2f : mu', ep);
        RowNames{9 + 8*k-5} = sprintf('ep = %.2f : Std of mu', ep);
        RowNames{9 + 8*k-4} = sprintf('ep = %.2f : mu 95 percent coverage', ep);
        RowNames{9 + 8*k-1} = sprintf('ep = %.2f : sigma uncorrected 95 percent coverage', ep);
        RowNames{9 + 8*k-3} = sprintf('ep = %.2f : sigma', ep);
        RowNames{9 + 8*k-2} = sprintf('ep = %.2f : Std of sigma', ep);
        RowNames{9 + 8*k-0} = sprintf('ep = %.2f : sigma 95 percent coverage', ep);
    end

    % Create table
    T = array2table(T, 'RowNames', RowNames, 'VariableNames', VariableNames);
    f = uifigure;
    uitable(f, 'Data', T);
    disp(T);
end

function table_normal_contamination_threshold(collection, collection_CI, collection_CI_u, n, ep_list, contamination_percentage_list, mu, sigma)
    numEps = length(ep_list);
    numRows = 8 * numEps + 9; % Adjusted for all rows
    numCols = length(contamination_percentage_list);
    T = ones(numRows, numCols);
    RowNames = cell(numRows, 1); % Cell array for row names

    % Initial Row Names
    RowNames{1} = 'Sample size';
    RowNames{2} = 'MLE_mu';
    RowNames{3} = 'MLE_mu 95 percent coverage';
    RowNames{4} = 'MLE_mu uncorrected 95 percent coverage';
    RowNames{5} = 'MLE_sigma';
    RowNames{6} = 'MLE_sigma 95 percent coverage';
    RowNames{7} = 'MLE_sigma uncorrected 95 percent coverage';
    RowNames{8} = 'Std of MLE_mu';
    RowNames{9} = 'Std of MLE_sigma';

    VariableNames = cell(1, numCols); % Cell array for variable (column) names
    for j = 1:numCols
        VariableNames{j} = sprintf('Contamination = %.2f', contamination_percentage_list(j));
    end

    % Populate initial rows
    T(1, :) = n;

    % Compute MLE_mu and MLE_sigma with trimming and standard deviations
    data_mle = squeeze(collection(:,:,end,:));    
    for j = 1:numCols
        % Data for MLE_mu and MLE_sigma
        data_mu = squeeze(data_mle(1,:,j));
        data_sigma = squeeze(data_mle(2,:,j));

        data_mu_trimmed = data_mu;

        T(2, j) = round(mean(data_mu_trimmed), 3); % Mean of trimmed MLE_mu
        T(8, j) = round(std(data_mu_trimmed), 3);  % Std of trimmed MLE_mu

        data_sigma_trimmed = data_sigma;

        T(5, j) = round(mean(data_sigma_trimmed), 3); % Mean of trimmed MLE_sigma
        T(9, j) = round(std(data_sigma_trimmed), 3);  % Std of trimmed MLE_sigma

        % Corrected Coverage
        data_mle_CI = squeeze(collection_CI(:,:,:,end,j));
        T(3, j) = round(mean(data_mle_CI(1,1,:) < mu & data_mle_CI(1,2,:) > mu), 3);
        T(6, j) = round(mean(data_mle_CI(2,1,:) < sigma & data_mle_CI(2,2,:) > sigma), 3);

        % Uncorrected Coverage
        data_mle_CI_u = squeeze(collection_CI_u(:,:,:,end,j));
        T(4, j) = round(mean(data_mle_CI_u(1,1,:) < mu & data_mle_CI_u(1,2,:) > mu), 3);
        T(7, j) = round(mean(data_mle_CI_u(2,1,:) < sigma & data_mle_CI_u(2,2,:) > sigma), 3);
    end

    % Compute for each epsilon in ep_list
    for k = 1:numEps
        ep = ep_list(k);    
        for j = 1:numCols
            % Data for mu and sigma
            data = squeeze(collection(:,:,k,:));
            data_mu = squeeze(data(1,:,j));
            data_sigma = squeeze(data(2,:,j));

            if k==1
                data_mu_trimmed = data_mu;
                data_sigma_trimmed = data_sigma;
            else
                data = squeeze(collection(:,:,1,:));
                mu_hat = mean(squeeze(data(1,:,j)));
                sigma_hat = mean(squeeze(data(2,:,j)));
                q1_mu = norminv(0.007,mu_hat,sigma_hat);
                q9_mu = norminv(0.995,mu_hat,sigma_hat);
                data_mu_trimmed = data_mu(data_mu >= q1_mu & data_mu <= q9_mu);
                data_sigma_trimmed = data_sigma(data_mu >= q1_mu & data_mu <= q9_mu);
            end
            
            T(9 + 8*k-6, j) = round(mean(data_mu_trimmed), 3); % Mean of trimmed mu
            T(9 + 8*k-5, j) = round(std(data_mu_trimmed), 3);  % Std of mu

            T(9 + 8*k-3, j) = round(mean(data_sigma_trimmed), 3); % Mean of trimmed sigma
            T(9 + 8*k-2, j) = round(std(data_sigma_trimmed), 3);  % Std of sigma

            % Coverage (corrected and uncorrected)
            data_CI = squeeze(collection_CI(:,:,:,k,j)); 
            data_CI_u = squeeze(collection_CI_u(:,:,:,k,j));
            T(9 + 8*k-4, j) = round(mean(data_CI(1,1,:) < mu & data_CI(1,2,:) > mu), 3); % Corrected mu coverage
            T(9 + 8*k-1, j) = round(mean(data_CI(2,1,:) < sigma & data_CI(2,2,:) > sigma), 3); % Corrected sigma coverage
            T(9 + 8*k-7, j) = round(mean(data_CI_u(1,1,:) < mu & data_CI_u(1,2,:) > mu), 3); % Uncorrected mu coverage
            T(9 + 8*k-0, j) = round(mean(data_CI_u(2,1,:) < sigma & data_CI_u(2,2,:) > sigma), 3); % Uncorrected sigma coverage
        end

        % Update RowNames
        RowNames{9 + 8*k-7} = sprintf('ep = %.2f : mu uncorrected 95 percent coverage', ep);
        RowNames{9 + 8*k-6} = sprintf('ep = %.2f : mu', ep);
        RowNames{9 + 8*k-5} = sprintf('ep = %.2f : Std of mu', ep);
        RowNames{9 + 8*k-4} = sprintf('ep = %.2f : mu 95 percent coverage', ep);
        RowNames{9 + 8*k-1} = sprintf('ep = %.2f : sigma uncorrected 95 percent coverage', ep);
        RowNames{9 + 8*k-3} = sprintf('ep = %.2f : sigma', ep);
        RowNames{9 + 8*k-2} = sprintf('ep = %.2f : Std of sigma', ep);
        RowNames{9 + 8*k-0} = sprintf('ep = %.2f : sigma 95 percent coverage', ep);
    end

    % Create table
    T = array2table(T, 'RowNames', RowNames, 'VariableNames', VariableNames);
    f = uifigure;
    uitable(f, 'Data', T);
    disp(T);
end

% function table_normal_contamination(collection,collection_CI,collection_CI_u,n,ep_list,contamination_percentage_list,mu,sigma)
%     numRows = 4 * length(ep_list) + 5;
%     numCols = length(contamination_percentage_list);
%     T = ones(numRows, numCols);
%     RowNames = cell(numRows, 1); % Cell array for row names
%     RowNames{1} = 'Sample size';
%     RowNames{2} = 'MLE_mu';
%     RowNames{3} = 'MLE_mu 95 percent coverage';
%     RowNames{4} = 'MLE_sigma';
%     RowNames{5} = 'MLE_sigma 95 percent coverage';
%     VariableNames = cell(1, numCols); % Cell array for variable (column) names
%     for j = 1:numCols
%         VariableNames{j} = sprintf('Contamination = %.2f', contamination_percentage_list(j));
%         data_mle_CI = squeeze(collection_CI(:,:,:,end,j));
%         T(3,j) = round(mean(data_mle_CI(1,1,:)<mu & data_mle_CI(1,2,:)>mu), 3);
%         T(5,j) = round(mean(data_mle_CI(2,1,:)<sigma & data_mle_CI(2,2,:)>sigma), 3);
%     end
%     T(1,:) = n;
% 
%     data_mle = squeeze(collection(:,:,end,:));    
%     T([2,4],:) = round(squeeze(mean(data_mle,2)),3);
%     for k = 1 : length(ep_list)
%         ep = ep_list(k);
%         data = squeeze(collection(:,:,k,:));
%         T([5+6*k-5, 5+6*k-2],:)=round(squeeze(mean(data,2)),3);
%         for j = 1:numCols
%             data_CI = squeeze(collection_CI(:,:,:,k,j)); 
%             T(5+6*k-4,j) = round(mean(data_CI(1,1,:)<mu & data_CI(1,2,:)>mu), 3);
%             T(5+6*k-1,j) = round(mean(data_CI(2,1,:)<sigma & data_CI(2,2,:)>sigma), 3);
%             data_CI_u = squeeze(collection_CI_u(:,:,:,k,j)); 
%             T(5+6*k-3,j) = round(mean(data_CI_u(1,1,:)<mu & data_CI_u(1,2,:)>mu), 3);
%             T(5+6*k,j) = round(mean(data_CI_u(2,1,:)<sigma & data_CI_u(2,2,:)>sigma), 3);
%         end
%         RowNames{5+6*k-5} = sprintf('ep = %.2f : mu', ep);
%         RowNames{5+6*k-4} = sprintf('ep = %.2f : mu 95 percent coverage', ep);
%         RowNames{5+6*k-3} = sprintf('ep = %.2f : mu uncorrected 95 percent coverage', ep);
%         RowNames{5+6*k-2} = sprintf('ep = %.2f : sigma', ep);
%         RowNames{5+6*k-1} = sprintf('ep = %.2f : sigma 95 percent coverage', ep);
%         RowNames{5+6*k} = sprintf('ep = %.2f : sigma uncorrected 95 percent coverage', ep);
%     end
%     T = array2table(T,'RowNames',RowNames,'VariableNames',VariableNames);
%     f = uifigure;
%     uitable(f, 'Data', T);
%     disp(T);
% end

function plot_normal(collection,mu,sigma,noise_K,iteration_number,ep_list,S_short)
    figure_mu = figure; % Create figure for mu
    hold on
    figure_sigma = figure; % Create figure for sigma
    hold on
    colors = lines(length(ep_list)+1);  % Predefined set of colors 
    colors(3, :) = [];
    idx = 1;  % Color index
    for k=1:length(ep_list)
        ep = ep_list(k);
        data = collection(:,:,1:S_short,k);

        ave_theta = median(data,3);
        ave_theta_75 = quantile(data,0.75,3);
        ave_theta_25 = quantile(data,0.25,3);  
        color = colors(idx, :);

        figure(figure_mu);
        plot(1:iteration_number, ave_theta(:, 1), 'Color', color, 'LineWidth', 1.5,'DisplayName', sprintf('epsilon = %.2f', ep));
        % plot(1:iteration_number, ave_theta_75(:, 1), 'Color', color, 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        % plot(1:iteration_number, ave_theta_25(:, 1), 'Color', color, 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        
        figure(figure_sigma);
        plot(1:iteration_number, ave_theta(:, 2), 'Color', color, 'LineWidth', 1.5,'DisplayName', sprintf('epsilon = %.2f', ep));
        % plot(1:iteration_number, ave_theta_75(:, 2), 'Color', color, 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        % plot(1:iteration_number, ave_theta_25(:, 2), 'Color', color, 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        
        idx = idx + 1;
    end

    figure(figure_mu);
    yline(mu,'HandleVisibility', 'off')
    xline(noise_K,'HandleVisibility', 'off')
    title('\mu')
    xlabel('iteration times')
    legend();
    hold off;

    figure(figure_sigma);
    yline(sigma,'HandleVisibility', 'off')
    xline(noise_K,'HandleVisibility', 'off')
    title('\sigma')
    xlabel('Iteration times')
    legend();
    hold off;
end

function plot_normal_traj(TRAJ,noise_K,iteration_number,ep_list,S_short)
    figure;    
    colors = lines(length(ep_list))*0.8;  % Predefined set of colors
    idx = 1;  % Color index
    for k=1:length(ep_list)
        ep = ep_list(k);
        data = log(TRAJ(:,1:S_short,k));

        ave_traj = mean(data,2);
        color = colors(idx, :);
        hold on
        plot(1:iteration_number, ave_traj, 'Color', color, 'LineWidth', 1.5,'DisplayName', sprintf('epsilon = %.2f', ep));
        ylim([-5,1])

        idx = idx + 1;
    end
    xline(noise_K,'HandleVisibility', 'off')
    title('Trajectory')
    xlabel('Iteration times')
    ylabel('log(||gradient||_2)')
    legend();
end

function [collection,Collection_CI,Collection_CI_u,TRAJ] = collection_normal(data,data_resample,method,S,h,p_n,noise_K,iteration_number,ep_list,eta,theta_start)
    collection = ones(iteration_number,2,S,length(ep_list));
    Collection_CI = ones(iteration_number,2,2,S,length(ep_list));
    Collection_CI_u = ones(iteration_number,2,2,S,length(ep_list));
    TRAJ = ones(iteration_number,S,length(ep_list));
    for k = 1:length(ep_list)
        ep = ep_list(k);
        if method == 1
            parfor i=1:S
                X = data(:,i);
                X_resample=data_resample(:,i);
                [hat_theta,CI,CI_u,traj] = simulation.private_gradient_normal(X,X_resample,h,p_n,noise_K,iteration_number,ep,eta,theta_start);
                collection(:,:,i,k) = hat_theta;
                TRAJ(:,i,k) = traj;
                Collection_CI(:,:,:,i,k) = CI;
                Collection_CI_u(:,:,:,i,k) = CI_u;
            end
        elseif method ==2
            parfor i=1:S
                X = data(:,i);
                X_resample=data_resample(:,i);
                [hat_theta,CI,CI_u,traj] = simulation.private_Nt_normal(X,X_resample,h,p_n,noise_K,iteration_number,ep,eta,theta_start);
                collection(:,:,i,k) = hat_theta;
                TRAJ(:,i,k) = traj;
                Collection_CI(:,:,:,i,k) = CI;    
                Collection_CI_u(:,:,:,i,k) = CI_u;
            end
        else
            sprintf('Choose correct method, method=1 for gradient descent, method=2 for Newton\n')
        end
    end

end

function [collection,collection_CI,collection_CI_u] = collection_normal_contaminatin(data_contain,data_contain_resample,method,S,h,p_n,noise_K,iteration_number,ep_list,eta,theta_start, contamination_percentage_list)
    % Last one is MLE in 'length(ep_list)+1'.
    collection = ones(2, S, length(ep_list)+1,length(contamination_percentage_list));
    collection_CI = ones(2,2, S, length(ep_list)+1,length(contamination_percentage_list));
    collection_CI_u = ones(2,2, S, length(ep_list),length(contamination_percentage_list));
    for j = 1:length(contamination_percentage_list)
        for k = 1:length(ep_list)   
            ep = ep_list(k); 
            if method == 1
                parfor i=1:S
                    X = data_contain(:,i,j);
                    X_resample = data_contain_resample(:,i,j);
                    [hat_theta,CI,CI_u,~] = simulation.private_gradient_normal(X,X_resample,h,p_n,noise_K,iteration_number,ep,eta,theta_start);
                    collection(:,i,k,j) = hat_theta(noise_K,:);
                    collection_CI(:,:,i,k,j) = CI(noise_K,:,:); 
                    collection_CI_u(:,:,i,k,j) = CI_u(noise_K,:,:);
                end
            elseif method ==2
                parfor i=1:S
                    X = data_contain(:,i,j);
                    X_resample = data_contain_resample(:,i,j);
                    [hat_theta,CI,CI_u,~] = simulation.private_Nt_normal(X,X_resample,h,p_n,noise_K,iteration_number,ep,eta,theta_start);
                    collection(:,i,k,j) = hat_theta(noise_K,:);
                    collection_CI(:,:,i,k,j) = CI(noise_K,:,:);
                    collection_CI_u(:,:,i,k,j) = CI_u(noise_K,:,:);
                end
            else
                sprintf('Choose correct method, method=1 for gradient descent, method=2 for Newton\n')
            end
        end
        % Claculate MLE
        parfor i = 1:S
            X = data_contain(:,i,j);
            [MLE,MLE_CI] = simulation.MLE_norm(X);
            collection(:, i, end, j) = reshape(MLE,1,2);
            collection_CI(:,:, i, end, j) = MLE_CI;
            % mle_mu = mean(X);
            % mle_sigma = sqrt(mean((X - mle_mu).^2));
            % collection(:, i, end, j) = [mle_mu,mle_sigma];
        end
    end
    
end

%-----------------------------------------------------------------------

function [hat_theta,CI,CI_u,traj] = private_gradient_normal(X,X_resample,h,p_n,noise_K,iteration_number,ep,eta,theta_start)
    % X is input data, shape (n,1)
    n = length(X);
    CI = ones(iteration_number,length(theta_start),2);
    CI_u = ones(iteration_number,length(theta_start),2);
    traj = ones(1,iteration_number);
    if ep>0 && ep<2
        f = @(e) OPT.privacy_composition_HD(e,noise_K)-ep;
        e = fzero(f,0.3);
        noise_sd = 2*sqrt(6)/sqrt(-8*log(1-0.5*e)) / (n^(1/p_n));
    else
        noise_sd = 0;
    end
    noise = normrnd(0,noise_sd,iteration_number,length(theta_start));

    K = @(x) OPT.epanechnikov(x);
    f_k = @(x) OPT.f_fb(x,X,K,h);
    f_p = @(x,theta) normpdf(x,theta(1),abs(theta(2)));

    %--------------------------------------------------------------------
    % For CI
    if ep>0 && ep<2
        f_CI = @(e) OPT.privacy_composition_HD(e,2)-ep;
        e = fzero(f_CI,0.3);
        noise_sd_W1 = sqrt(118) / (n^(1/p_n)) /sqrt(-8*log(1-0.5*e));
        noise_sd_W2 = ( 2*sqrt(6)/ (n^(1/p_n)) )^2 /sqrt(-8*log(1-0.5*e)) ;
    else
        noise_sd_W1 = 0;
        noise_sd_W2 = 0;
    end

    %--------------------------------------------------
  
    argmin = ones(iteration_number,length(theta_start));
    theta = theta_start;
    for i=1:iteration_number
        grad_theta = 2 * reshape(OPT.Grad_HD_norm(X_resample,f_k,f_p,theta),1,length(theta));
        traj(i) = norm(grad_theta);
        theta = theta - eta * (grad_theta + noise(i,:)/theta(2));
        argmin(i,:) = theta;
        % grad_theta = NOPT.grad(HD,theta);
        % theta = theta - eta * (grad_theta + noise(i,:)/theta(2));        
        % argmin(i,:) = theta;

        CI_theta = reshape(theta,2,1);
        % Generate noise_W1, noise_W2    
        % Create an upper triangular matrix with i.i.d. normal elements
        R_1 = triu(normrnd(0, noise_sd_W1/(theta(2)^2), length(theta_start), length(theta_start)));  
        R_2 = triu(normrnd(0, noise_sd_W2/(theta(2)^2), length(theta_start), length(theta_start)));
        % Symmetrize the matrix
        W1 = R_1 + R_1' - diag(diag(R_1));
        W2 = R_2 + R_2' - diag(diag(R_1));
        M = 2 * OPT.Hess_HD_norm(X_resample,f_k,f_p,theta);
        % M = NOPT.hess(HD,theta);
        M_nois = M + W1;
        M_nois = OPT.regularize_hessian(M_nois,0.1); 
        M_nois_inv = inv(M_nois);
        Q = 4*OPT.Grad2_HD_norm(X_resample,f_k,f_p,theta);
        Q_nois = (Q + W2)/4;
        Q_nois = OPT.regularize_hessian(Q_nois,0);
        V_nois = diag(M_nois_inv * Q_nois * M_nois_inv);
        CI_u(i,:,1) = CI_theta - norminv(0.975)*sqrt( V_nois/n );
        CI_u(i,:,2) = CI_theta + norminv(0.975)*sqrt( V_nois/n );
        CI(i,:,1) = CI_theta - norminv(0.975)*sqrt( V_nois/n + 2*(eta*noise_sd/theta(2))^2 );
        CI(i,:,2) = CI_theta + norminv(0.975)*sqrt( V_nois/n + 2*(eta*noise_sd/theta(2))^2 );
    end

    hat_theta = argmin;    

end

% function argmin = NGD_normal(obj,theta_start,iteration_number,step_size,noise)
%     argmin = ones(iteration_number,length(theta_start));
%     theta = theta_start;
%     K = iteration_number;
%     eta = step_size;
%     for i=1:K
%         grad_theta = NOPT.grad(obj,theta);
%         theta = theta - eta * (grad_theta + noise(i,:)/theta(2));        
%         argmin(i,:) = theta;
%    end
% end

% --------------------------------------------------------------------

function [hat_theta,CI,CI_u,traj] = private_Nt_normal(X,X_resample,h,p_n,noise_K,iteration_number,ep,eta,theta_start)
    % X is input data, shape (n,1)
    n = length(X);
    CI = ones(iteration_number,length(theta_start),2);
    CI_u = ones(iteration_number,length(theta_start),2);
    traj = ones(1,iteration_number);
    if ep>0 && ep<2
        f = @(e) OPT.privacy_composition_HD(e,noise_K)-ep;
        e = fzero(f,0.3);
        noise_sd_G = 2*sqrt(6)/sqrt(-8*log(1-0.25*e)) / (n^(1/p_n));
        noise_sd_H = sqrt(118)/sqrt(-8*log(1-0.25*e)) / (n^(1/p_n));
    else
        noise_sd_G = 0;
        noise_sd_H = 0;
    end
    noise_G = normrnd(0,noise_sd_G,iteration_number,length(theta_start));

    K = @(x) OPT.epanechnikov(x);
    f_k = @(x) OPT.f_fb(x,X,K,h);
    f_p = @(x,theta) normpdf(x,theta(1),abs(theta(2)));
    
    %--------------------------------------------------------------------
    % For CI
    if ep>0 && ep<2
        f_CI = @(e) OPT.privacy_composition_HD(e,2)-ep;
        e = fzero(f_CI,0.3);
        noise_sd_W1 = sqrt(118) / (n^(1/p_n)) /sqrt(-8*log(1-0.5*e));
        noise_sd_W2 = ( 2*sqrt(6)/ (n^(1/p_n)) )^2 /sqrt(-8*log(1-0.5*e));
    else
        noise_sd_W1 = 0;
        noise_sd_W2 = 0;
    end

    %--------------------------------------------------
    % Preallocate noise_H
    noise_H = zeros(iteration_number, length(theta_start), length(theta_start));
    % Generate noise_H
    for k = 1:iteration_number
        % Create an upper triangular matrix with i.i.d. normal elements
        R_1 = triu(normrnd(0, noise_sd_H, length(theta_start), length(theta_start)));  
        % Symmetrize the matrix
        noise_H(k, :, :) = R_1 + R_1' - diag(diag(R_1));
    end

    argmin = ones(iteration_number,length(theta_start));
    theta = theta_start; % theta is row 
    for i=1:iteration_number
        grad_theta = 2 * reshape(OPT.Grad_HD_norm(X_resample,f_k,f_p,theta),1,length(theta));
        traj(i) = norm(grad_theta);
        % grad_theta = NOPT.grad(HD,theta);
        % if grad_theta == 0
        %     theta = theta_start;
        %     grad_theta = NOPT.grad(HD,theta);
        % end
        Hess = 2 * OPT.Hess_HD_norm(X_resample,f_k,f_p,theta);
        % Hess = NOPT.hess(HD,theta);
        Hess_nois = Hess + squeeze(noise_H(i,:,:))./(theta(2)^2);
        Hess_nois = OPT.regularize_hessian(Hess_nois,0.1);
        Hess_nois_T = transpose(Hess_nois);
        theta = theta - eta *(grad_theta + noise_G(i,:)/theta(2))/Hess_nois_T;
        argmin(i,:) = theta;

        CI_theta = transpose(theta);
        % Generate noise_W1, noise_W2    
        % Create an upper triangular matrix with i.i.d. normal elements
        R_1 = triu(normrnd(0, noise_sd_W1/(theta(2)^2), length(theta_start), length(theta_start)));  
        R_2 = triu(normrnd(0, noise_sd_W2/(theta(2)^2), length(theta_start), length(theta_start)));
        % Symmetrize the matrix
        W1 = R_1 + R_1' - diag(diag(R_1));
        W2 = R_2 + R_2' - diag(diag(R_1));
        M = Hess;
        M_nois = M + W1;
        M_nois = OPT.regularize_hessian(M_nois,0.1); 
        M_nois_inv = inv(M_nois);
        Q = 4 * OPT.Grad2_HD_norm(X_resample,f_k,f_p,theta);
        Q_nois = (Q + W2)/4;
        Q_nois = OPT.regularize_hessian(Q_nois,0);
        V_nois = diag(M_nois_inv * Q_nois * M_nois_inv);
        
        CI_u(i,:,1) = CI_theta - norminv(0.975)*sqrt( V_nois/n );
        CI_u(i,:,2) = CI_theta + norminv(0.975)*sqrt( V_nois/n );
        Hess_nois_inv = inv(Hess_nois);
        Sig = eta^2* (noise_sd_G)^2 * Hess_nois_inv*Hess_nois_inv;
        CI(i,:,1) = CI_theta - norminv(0.975)*sqrt( V_nois/n + diag(Sig) );
        CI(i,:,2) = CI_theta + norminv(0.975)*sqrt( V_nois/n + diag(Sig) );
    end

    hat_theta = argmin;
end

% function argmin = NNT_normal(obj,theta_start,iteration_number,step_size,noise_G,noise_H)
%             argmin = ones(iteration_number,length(theta_start));
%             theta = theta_start; % theta is row 
%             K = iteration_number;
%             eta = step_size;
%             for i=1:K
%                 grad_theta = NOPT.grad(obj,theta);
%                 if grad_theta == 0
%                     theta = theta_start;
%                     grad_theta = NOPT.grad(obj,theta);
%                 end
%                 Hess = NOPT.hess(obj,theta);
%                 Hess_nois = Hess + squeeze(noise_H(i,:,:))./(theta(2)^2);
%                 Hess_nois = NOPT.regularize_hessian(Hess_nois,0.5);
%                 Hess_nois_T = transpose(Hess_nois);
%                 theta = theta - eta *(grad_theta + noise_G(i,:)/theta(2))/Hess_nois_T;
%                 argmin(i,:) = theta;
%             end
% end

% Claculate sandwich variance for MLE
function [theta,CI] = MLE_norm(X)
    n = length(X);
    theta = ones(2,1);
    CI = ones(2,2);
    theta(1) = mean(X);
    theta(2) = std(X);

    DfDmu = @(x)  (x-theta(1)) ./ theta(2)^2;
    DfDsigma = @(x)  ( (x-theta(1)).^2-theta(2)^2 )/ theta(2)^3;
    u = @(x) [DfDmu(x) ; DfDsigma(x)];
    H_f = @(x) [((x-theta(1))^2-theta(2)^2) / theta(2)^4, (x-theta(1))*((x-theta(1))^2-3*theta(2)^2)/theta(2)^5;
        (x-theta(1))*((x-theta(1))^2-3*theta(2)^2)/theta(2)^5, ( (x-theta(1))^4-5*theta(2)^2*(x-theta(1))^2+2*theta(2)^4 )/theta(2)^6];
    s = zeros(2,2,n);
    for i=1:n
        xi = X(i);
        s(:,:,i) = ( H_f(xi)- u(xi) * u(xi)'  );
    end
    Hess = sum(s,3);
    Hess_inv = inv(Hess);

    s = zeros(2,2,n);
    for i=1:n
        xi = X(i);
        s(:,:,i) =  u(xi)  * u(xi)' ;
    end
    grad2 = sum(s,3);

    V = diag(Hess_inv*grad2*Hess_inv);
    CI(:,1) = theta - norminv(0.975)*sqrt( V );
    CI(:,2) = theta + norminv(0.975)*sqrt( V );
end

    end
end