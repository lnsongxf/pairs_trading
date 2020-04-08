%% General Dataset Adjustments
% The 'Full_Data.mat' has been imported using CRSP from WRDS database. This dataset includes the
% variables: NCUSIP COMNAM TICKER HEXCD PRC VOL CFACPR for daily observations between 29Dec1961 - 
% 31Dec2018. The only conditional statement imposed in the dataset is hexcd = 1 OR hexcd = 2 OR 
% hexcd = 3; which is intended to limit the analysis of stocks traded on AMEX, NYSE and NASDAQ.

% The data set 'Full_Data.mat' can be found at 'www.talespadilha.com/reseach'

clear all
clc

% Loading the data for the full period:
fprintf( 'Reading data... ' );
load('Full_Data.mat')
fprintf('Done!\n');

% Brief correction of prices (CRSP data files have negative prices to indicate mean between bid and
% ask when the are no trades):
Data.PRC = abs(Data.PRC);

% Cleaning any missing value:
Data = Data(~isnan(Data.PRC), :);

% Calculating Stock USD Volume:
USD_VOL = Data.PRC .* Data.VOL;
Data = [Data, table(USD_VOL)];
clear USD_VOL

% Adjusting prices for events:
ADJ_PRC = Data.PRC ./ Data.CFACPR;
Data = [Data, table(ADJ_PRC)];
clear ADJ_PRC


%% Building the large matrices with all the USD_VOLs and Prices
% To begin our analysis, we will build a large matrix with all the USD_VOLs and another one with all
% the ADJ_PRC, where the rows will indicate the dates and the columns will indicate the PERMNOs.

% We will use these vectors with PERMNOs and Dates to guide us through the large matrices that will 
% include the given data being analysed for each stock/date:
stocks = unique(Data.PERMNO);
dates = unique(Data.date);

% Creating an index for the dates (will be used to assign stocks info to their dates):
minDate = min(dates);
shiftedDates = dates - minDate + 1;
dates_LT = zeros(max(shiftedDates), 1);
dates_LT(shiftedDates) = 1 : numel(shiftedDates);

% Creating the big matrix of Volumes and Prices to be filled:
volumes = nan(numel(dates), numel(stocks));
prices = nan(numel(dates), numel(stocks));

% Filling the 'volumes' and 'prices' matrices (looping by stock):
x = numel(stocks);
for s_id = 1 : x
    if mod(s_id, 100) == 0
    fprintf('Building %d of %d...\n', s_id, x);
    end
    % Building an index for the stock being analysed in each step for each row in the data:
    r_id = table2array(Data(:,1)) == stocks(s_id);
    % Getting the date for these rows:
    r_dates = table2array(Data(r_id,2));
    % Using the index that we have built to get a vector of indices mapping to each row in the 
    % matrix volumes:
    data_id = dates_LT(r_dates - minDate + 1) ;
    % Filling the volume matrix:
    volumes(data_id,s_id) = table2array(Data(r_id,9));
    % And finally the prices matrix:
    prices(data_id,s_id) = table2array(Data(r_id,10));
end
clear x  data_id r_dates r_id s_id shiftedDates Data
% Note that the table data won't be saved. This is fine, we won't use it anymore. 
save('First_Step_Full', '-v7.3')
fprintf('Done!\n');


%% Actually cleaning the data
% The aim in this section is to fill the matrix with the Stock IDs and PERMNOs for all stocks to be 
% considered in the pairs trading formation, for each of the formation periods.
fprintf('Cleaning dataset...');

% Getting the indices for t when the months change:
delta_t = [0; diff(mod(dates,100))];
t_month_first = [];
for i = 1:numel(dates)
    if delta_t(i) < 0
        t_month_first = [t_month_first; i];
    end 
end   
clear i delta_t

% Creating the "vectors" where we will input the PERMNO of the stocks and their location indices
% to be considered for pairs formation in each of the formation periods:
PERMNOs = [];
Stocks_id = [];

% Looping for all pair formation sets:
for m = 13:(numel(t_month_first)-5)
    % Restraining the training set:
    t_training_months = t_month_first(m-12:m);
    training_date_low = dates(t_training_months(1)-1); 
    training_date_high = dates(t_training_months(end)-1);
    
    % Using our large matrix to get the data on volume:
    training_date_low_id = dates_LT(training_date_low - minDate + 1);
    training_date_high_id = dates_LT(training_date_high - minDate + 1);
    training_volume = volumes(training_date_low_id:training_date_high_id, :);
    training_prices = prices(training_date_low_id:training_date_high_id, :);
    
    % Leaving only the most traded stocks:
    total_volume = sum(sum(training_volume), 'omitnan');
    [ord_vol, stock_id] = sort(sum(training_volume)', 'descend','MissingPlacement', 'last');
    cap = 0;
    i = 1;
    training_stocks_id =[];
    while cap < 0.9 
        % Not considering stocks with adjusted price lower than one dollar:
        if min(training_prices(:,stock_id(i)))>1
            cap = cap + (ord_vol(i) / total_volume);
            training_stocks_id = [training_stocks_id; stock_id(i)];
        end
        i = i + 1;
    end
   %training_stocks_id = stock_id(1:i);
   clear i cap ord_vol total_volume stock_id
   
   %Using 'training_stocks_id' to get the PERMOs:
   Stocks_id = [Stocks_id; {training_stocks_id}];
   PERMNOs = [PERMNOs; {stocks(training_stocks_id)}];
end
clear training_volume training_prices training_stocks_id training_date_high_id training_date_low_id training_date_high training_date_low t_training_months m
save('Clean_Data', '-v7.3')
fprintf('Done!\n');

