function [cdf_lifetime_distribution, pdf_lifecycle] = Gamma_Process(R_0, threshold, theta, a, b, end_of_time_series, timesteps_per_year)

% Input
x_t_max = end_of_time_series;             % Maximum deterioration to focus on
time = 1:1/timesteps_per_year:end_of_time_series; % year
x_t = 1:x_t_max/length(time):x_t_max; % Deterioration at time t
N_random_values_percentiles = 10000; % 100000 values for straight line

% Scale parameter beta > 1
beta = theta;

% Limit state
limit_state = R_0 - threshold; % Limit state: Value = maximal accepted deterioration

%%
for year = 1:1:length(time)
    t = time(1,year);
    
% Shape parameter alfa > 1
alfa = (a.*power(t,b))/theta;
    
% Gamma distribution
p(year,:) = gampdf(x_t,alfa,beta); %#ok<AGROW>

% Confidence intervals
Dataset (year,:)=gamrnd(alfa,beta,[1 N_random_values_percentiles]); %#ok<AGROW>
Con_95(year,:) = prctile(Dataset (year,:),95);%#ok<AGROW>
Con_5(year,:) = prctile(Dataset (year,:),5);%#ok<AGROW>

% Determining percentiles
Expectation_Deterioration_at_time_t (year,:) = alfa.*beta;%#ok<AGROW>
R_time_expectation(year,:) = R_0 - Expectation_Deterioration_at_time_t (year,:);%#ok<AGROW>
R_time_95(year,:) = R_0 - Con_95(year,:);%#ok<AGROW>
R_time_5(year,:) = R_0 - Con_5(year,:);%#ok<AGROW>

% Lifetime Distribution
cdf_lifetime_distribution(year,:) = gammainc(limit_state/theta,alfa,'upper');%#ok<AGROW>
pdf_lifecycle(year,:) = gampdf(limit_state,alfa,beta);%#ok<AGROW>
end

R_time_expectation = [0; R_time_expectation];
R_time_95 = [0; R_time_95];
R_time_5 = [0; R_time_5];
time = [0 time];


%%
figure
xlabel('Depth scour hole (m)')
ylabel('Pr(X(t) = x)')
xlim([0 end_of_time_series])
for yearnumbers = 1:1:end_of_time_series
col(yearnumbers,:) = find(time == yearnumbers);%#ok<AGROW>
for row = 1:1:length(col)
hold on
plot(x_t,p(row,:),'Color','black')
end
end
hold off

figure
plot(time,R_time_expectation(:,1))
xlim([0 50])
hold on
plot(time,R_time_95(:,1),':black')
hold on
plot(time,R_time_5(:,1),':black')
hold on
line([0,end_of_time_series],[threshold threshold],'Color', 'r')
hold off

figure
plot(time,[0;cdf_lifetime_distribution])
xlim([0 end_of_time_series])
xlabel('Year')
ylabel('Pr(T \leq t) = Pr(X(t) \geq y)')

figure
plot(time,[0;pdf_lifecycle])
xlim([0 end_of_time_series])
xlabel('Year')
ylabel('Pr(T = t)')