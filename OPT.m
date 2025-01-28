classdef OPT
    methods(Static)

        % Epanechnikov kernel
        function K = epanechnikov(x)
            K = (3/4) * (1 - x.^2) .* (abs(x) <= 1);
        end

        % Calculate Silverman's bandwidth, H is the data set
        function c = silverman(H)
            n = length(H);
            q3 = quantile(H,0.75);
            q1 = quantile(H,0.25);
            IQR = q3-q1;
            c = 0.9 * min([std(H),IQR/1.34]) * n^(-1/5);
        end
        
        % Basic fixed bandwidth kernel density
        function f = f_fb(x,H,K,h)
            n = length(H);
            H = reshape(H,1,n);
            f = 0; 
            for i=1:n
                f = f + K((x-H(i))/h);
            end
            f = f/(n*h);
        end

        % Generate m_n samples from Epanechnikov kernel density
        function X_epan = rnd_epan(m_n)
            U = rand([3, m_n]); 
            X_epan = zeros(1, m_n);
            X_epan(U(1,:) > 0.5) = U(2,U(1,:) > 0.5).^(1/3) .* U(3,U(1,:) > 0.5);
            X_epan(U(1,:) <= 0.5) = -U(2,U(1,:) <= 0.5).^(1/3) .* U(3,U(1,:) <= 0.5);
        end

        % Generate from f_k (Epanechnikov kernel density estimator)
        % X is original data, h is bandwidth, m_n is resample size
        function X_resample = rnd_f_k_epa(m_n, X, h)
            % Inputs:
            % m_n: Number of samples to generate
            % X: Input vector of length n
            % h: Bandwidth parameter

            % Generate random values using the Epanechnikov kernel
            Epan = OPT.rnd_epan(m_n); 

            % Generate uniform random numbers
            U = rand([1, m_n]);

            % Initialize X_resample
            X_resample = zeros(1, m_n);

            % Length of input vector X
            n = length(X);

            % Loop through each random value in U
            for j = 1:m_n
                % Find i such that U(j) falls in ((i-1)/n, i/n)
                i = ceil(U(j) * n); % Convert U(j) to index (i)

                % Assign the resampled value
                X_resample(j) = Epan(j) * h + X(i);
            end
        end

        % HD = \int (\sqrt{g_n}-\sqrt{f_\theta})^2
        function HD = HD(X,f_k,f_p,theta)
            n = length(X);
            X = reshape(X, 1, n);
            s = sum(sqrt( f_p(X,theta)./f_k(X) ),2);
            HD = 2-2*s/n;
        end

        % This is for HD = \int (\sqrt{g_n}-\sqrt{f_\theta})^2
        function grad = Grad_HD_norm(X,f_k,f_p,theta)
            n = length(X);
            X = reshape(X, 1, n);
            DfDmu = @(x)  (x-theta(1)) ./ theta(2)^2;
            DfDsigma = @(x)  ( (x-theta(1)).^2-theta(2)^2 )/ theta(2)^3;
            u = @(x) [DfDmu(x) ; DfDsigma(x)];
            s = sum(sqrt( f_p(X,theta)./f_k(X) ) .* u(X),2);
            grad = -s/n;
        end
        % function grad = Grad_HD_norm(X,f_k,f_p,theta)
        %     n = length(X);
        %     X = reshape(X, 1, n);
        %     DfDmu = @(x)  (x-theta(1)) ./ theta(2)^2;
        %     DfDsigma = @(x)  ( (x-theta(1)).^2-theta(2)^2 )/ theta(2)^3;
        %     u = @(x) [DfDmu(x) ; DfDsigma(x)];
        %     s = [0;0];
        %     for i=1:n
        %         xi = X(i);
        %         s = s + sqrt( f_p(xi,theta)/f_k(xi) ) * u(xi);
        %     end
        %     grad = -s/n;
        % end

        % This is for HD = \int (\sqrt{g_n}-\sqrt{f_\theta})^2
        function grad2 = Grad2_HD_norm(X,f_k,f_p,theta)
            n = length(X);
            X = reshape(X, 1, n);
            DfDmu = @(x)  (x-theta(1)) ./ theta(2)^2;
            DfDsigma = @(x)  ( (x-theta(1)).^2-theta(2)^2 )/ theta(2)^3;
            u = @(x) [DfDmu(x) ; DfDsigma(x)];
            s = zeros(2,2,n);
            for i=1:n
                xi = X(i);
                s(:,:,i) = f_p(xi,theta)/f_k(xi) * u(xi)  * u(xi)' ;
            end
            grad2 = sum(s,3)/n;
        end
        
        % This is for HD = \int (\sqrt{g_n}-\sqrt{f_\theta})^2
        function Hess = Hess_HD_norm(X,f_k,f_p,theta)
            n = length(X);
            X = reshape(X, 1, n);
            DfDmu = @(x)  (x-theta(1)) ./ theta(2)^2;
            DfDsigma = @(x)  ( (x-theta(1)).^2-theta(2)^2 )/ theta(2)^3;
            u = @(x) [DfDmu(x) ; DfDsigma(x)];
            H_f = @(x) [((x-theta(1))^2-theta(2)^2) / theta(2)^4, (x-theta(1))*((x-theta(1))^2-3*theta(2)^2)/theta(2)^5;
                (x-theta(1))*((x-theta(1))^2-3*theta(2)^2)/theta(2)^5, ( (x-theta(1))^4-5*theta(2)^2*(x-theta(1))^2+2*theta(2)^4 )/theta(2)^6];
            s = zeros(2,2,n);
            for i=1:n
                xi = X(i);
                s(:,:,i) = sqrt( f_p(xi,theta)/f_k(xi) ) * ( 0.5* u(xi) * u(xi)' -H_f(xi) );
            end
            Hess = sum(s,3)/n;
        end
        % function Hess = Hess_HD_norm2(X,f_k,f_p,theta)
        %     n = length(X);
        %     X = reshape(X, 1, n);
        %     DfDmu = @(x)  (x-theta(1)) ./ theta(2)^2;
        %     DfDsigma = @(x)  ( (x-theta(1)).^2-theta(2)^2 )/ theta(2)^3;
        %     u = @(x) [DfDmu(x) ; DfDsigma(x)];
        %     H_f = @(x) [((x-theta(1)).^2-theta(2).^2) / theta(2)^4, (x-theta(1)).*((x-theta(1)).^2-3*theta(2)^2)/theta(2)^5;
        %         (x-theta(1)).*((x-theta(1)).^2-3*theta(2)^2)/theta(2)^5, ( (x-theta(1)).^4-5*theta(2)^2*(x-theta(1)).^2+2*theta(2)^4 )/theta(2)^6];
        %     s = sum(sqrt( f_p(X,theta)./f_k(X) ) .* ( 0.5* u(X) * u(X)' -H_f(X) ),2);
        %     Hess = s/n;
        % end

        function H_reg = regularize_hessian(H, c)
            % Symmetrize H
            H = (H + H') / 2;
            [V, D] = eig(H);
            D = diag(D);
            D(D < c) = c;
            H_reg = V * diag(D) * V';
        end

        % HD = \int (\sqrt{g_n}-\sqrt{f_\theta})^2
        function argmin = GD(X,f_k,f_p,theta_start,iteration_number,step_size)
            argmin = ones(iteration_number,length(theta_start));
            theta = theta_start;
            K = iteration_number;
            eta = step_size;
            for i=1:K
                grad_theta = reshape(OPT.Grad_HD_norm(X,f_k,f_p,theta),1,length(theta));
                theta = theta - eta * (grad_theta);
                argmin(i,:) = theta;
            end
        end

        % HD = \int (\sqrt{g_n}-\sqrt{f_\theta})^2
        function argmin = NT(X,f_k,f_p,theta_start,iteration_number,step_size)
            argmin = ones(iteration_number,length(theta_start));
            theta = theta_start; % theta is row 
            K = iteration_number;
            eta = step_size;
            for i=1:K
                grad_theta = reshape(OPT.Grad_HD_norm(X,f_k,f_p,theta),1,length(theta));
                if grad_theta == 0
                    theta = theta_start;
                end
                Hess = OPT.Hess_HD_norm(X,f_k,f_p,theta);
                Hess = OPT.regularize_hessian(Hess,0.1);
                Hess = transpose(Hess);
                theta = theta - eta * grad_theta / Hess;
                argmin(i,:) = theta;
            end
        end

         % Calculate HD composition privacy
        function e_composition = privacy_composition_HD(e,n)
            if n==1
                e_composition = e;
            else 
                e_pre = OPT.privacy_composition_HD(e,n-1);
                e_composition = e + e_pre - 0.5*e*e_pre;
            end
        end

        % Calculate composition privacy
        function e_composition = privacy_composition_PDP(lam,e,n)
            if n==1
                e_composition = e;
            else 
                e_pre = OPT.privacy_composition_PDP(lam,e,n-1);
                e_composition = e + e_pre + lam*(lam+1)*e*e_pre;
            end
        end

    end
end