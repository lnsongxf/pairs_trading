%% Pre-Processing Before Analysis
% Before analysing the relationship between daily market returns and daiy Pairs Trading returns, we
% need to construct the series of daily Pairs Trading returns for each of the strategies being
% considered in this study.

clear all 
clc

%% Uploading Data for Portfolio

% Uploading data for whichever portfolio we want to analyse
fprintf( 'Reading data... ' );
load Results_DM_40
fprintf('Done!\n');

%% Processing Returns for this Portfolio

% First, joining all profits in tha same matrix:
all_profits = [];
for i=1:numel(Profits)
   Portfolio = Profits{i};
   for j=1:numel(Portfolio)
       x = Portfolio{j};
       if ~isempty(x)
           all_profits = [all_profits; x];
       end  
   end
end
clear i j x Portfolio

% Now using this matrix with all profits to get profits by day:
AP_Dates = all_profits(:,3);
AP_Returns = all_profits(:,1);
%Lets loop for all days where we have profit info:
days = unique(AP_Dates);
daily_returns = nan(numel(days),2); 
for i =1:numel(days)
   x = AP_Returns(AP_Dates == days(i));
   day_profit = mean(x);
   daily_returns(i, 1) = dates(days(i));
   daily_returns(i, 2) = day_profit;
end
clear x i days AP_Returns AP_Dates

% Now getting cumulative returns:
cum_returns = nan(size(daily_returns));
cum_returns(:,1) = daily_returns(:,1);
cum_returns(:,2) = cumprod(1+daily_returns(:,2));

% Saving daily returns:
save('Daily_r_DM_40', 'daily_returns', 'cum_returns')