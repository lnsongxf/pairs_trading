%% Summary of Statistics - Daily Portfolio Returns
% Lets first compare the performance of the different types of Pairs Trading strategies being
% analysed in this study.

clear all
clc

%% Defining the subset we will study:
% Selecting whichever portfolios we are going to study:
port_set = ["DM_5", "DM_10", "DM_20", "DM_40"];
% Selection whatever titles we want for the figures:
titles = ["Distance Method - 5 Pairs", "Distance Method - 10 Pairs", ...
    "Distance Method - 20 Pairs", "Distance Method - 40 Pairs"];

%% Loading the Data

fprintf('Reading data... ')
for portfolio = port_set
    filename = strcat('Daily_r_',portfolio);
    load(filename)   
    r.(portfolio) = daily_returns;
    cum_r.(portfolio) = cum_returns;
end
fprintf('Done!\n')

%% Graph of Daily Portfolio Returns:

d = num2str(r.DM_5(:, 1)); 
date = datetime(d,'InputFormat','yyyyMMdd','Format','dd/MM/yyyy');

figure()
i=0;
for portfolio = port_set
    i=i+1;
    date_id = size(r.(portfolio), 1) - size(r.DM_5, 1) +1;
    subplot(2,2,i);
    plot(date, r.(portfolio)(date_id:end,2))
    title(titles(i))
    ylim([-0.05, 0.05])
end

%% Rolling Sample Standard Deviation - Preliminary Analysis of Variance

% Set the size of the sample in which we will compute the standard deviation:
n = 30;

% Calculating and plotting the rolling standard deviation:
i = 0;
x = date(n:end);
figure()
for portfolio = port_set
    sds_r.(portfolio) = NaN(numel(r.(portfolio)(:,2))-n+1,1);
    for j=1:numel(sds_r.(portfolio))
        sds_r.(portfolio)(j) = std(r.(portfolio)(j:(j+n-1),2));
    end
    % Plotting the series:
    i = i +1; 
    date_id = size(sds_r.(portfolio), 1) - size(sds_r.DM_5, 1) +1;
    subplot(2,2,i);
    plot(x, sds_r.(portfolio)(date_id:end))
    ylim([0, 0.025])
    title(titles(i))
end

%% Checking the Performance of the Different Pairs Trading Strategies

% Loading SP Data
load SP_Index
1

%% Returns per year:

% Calculating and plotting returns per year:
i = 0;
figure()
for portfolio = port_set
    % Creading the location index for the first day of the year:
    delta_y = [0; diff(floor(r.(portfolio)(:, 1)./10000))];
    d_year_first =[1];
    for j = 1:numel(r.(portfolio)(:, 1))
        if delta_y(j) > 0
            d_year_first = [d_year_first; j];
        end 
    end
    d_year_first = [d_year_first; numel(r.(portfolio)(:, 1))+1];
    % Calculating returns per year:
    year_returns = nan(numel(d_year_first)-1,2);
    for j=1:(numel(d_year_first)-1)
        year_returns(j,2) = (prod(1+r.(portfolio)(d_year_first(j):(d_year_first(j+1)-1),2))-1);
        year_returns(j,1) = floor(r.(portfolio)(d_year_first(j), 1)/10000);
    end
    % Plotting 'year_returns':
    i=i+1;
    subplot(2,2,i);
    bar(year_returns(:,1), year_returns(:,2))
    ylim([-0.1, 0.7])
    title(titles(i))
end

%% Histogram of Daily Portfolio Returns:

i=0;
figure()
for portfolio = port_set
    i=i+1;
    subplot(2,2,i);
    histogram(r.(portfolio)(:,2))
    xlim([-0.02, 0.02])
    title(titles(i))    
end 
