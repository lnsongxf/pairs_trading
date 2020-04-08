%% Pairs Trading - Co-Integration Approach
% We will use the clean data (see Data_Cleaning.m for info on how to clean the data) to identify and
% implement Pairs Trading following Gatev (2006).
clear all
clc

%% Setting the size of the portfolio
% Before we proceed to identify the pairs, we need to define how many pairs our portfolio of pairs
% will have. This will be defined by the variable 'P'. For example, setting P=20:
P=5;

%% Matching Pairs
fprintf( 'Reading data... ' ) ;
load('Clean_Data.mat')
fprintf('Done!\n');

% Lets now build a matrix with the normalized prices for each training set and get the index of the 
% pairs:
Pairs_id = [];
Pairs_Betas_j = [];
Pairs_Spreads = [];
for m = 13:(numel(t_month_first)-5)
    if mod(m-12, 10) == 0
    fprintf('Processing %d of %d training sets...\n', m-12, numel(t_month_first)-17);
    end
    % Restraining the training set:
    t_training_months = t_month_first(m-12:m);
    training_date_low = dates(t_training_months(1)-1); 
    training_date_high = dates(t_training_months(end)-1);
  
    % Getting the indices for the dates of the training set:
    training_date_low_id = dates_LT(training_date_low - minDate + 1);
    training_date_high_id = dates_LT(training_date_high - minDate + 1);
    
    % Now getting the indices of the stocks:
    s_id = Stocks_id{m-12};
    
    % Building the matrix of normalized prices for this sample:
    norm_prices = prices(training_date_low_id:training_date_high_id, s_id)./prices(training_date_low_id, s_id);

    % Using these values to build the matrix of SSDs:
    SSDs = nan(numel(s_id) , numel(s_id));
    for i = 1:numel(s_id)  
        Dif = norm_prices - norm_prices(:, i);
        SSDs(i,:) = sum(Dif.^2);
    end
    SSDs = triu(SSDs);
    SSDs(SSDs==0)= NaN;    
    
    % Getting the info for the top P pairs:
    i_id = nan(P, 1);
    j_id = nan(P, 1);
    betas = nan(P, 1); %Betas from the co-integration regression
    spreads = nan(P, 2); %Spreads mean and standard deviation
    k = 1;
    while k < (P+1) 
       [A, i] = min(SSDs);
       [b, j] = min(A);
       SSDs(:, i(j)) = NaN;
       SSDs(i(j), :) = NaN;
       SSDs(:, j) = NaN;
       SSDs(j, :) = NaN;
       
       % Testing for co-integration between i(j) and j (no constant/no trend):
       Y = [norm_prices(:,i(j)), norm_prices(:,j)];
       [h,~,~,~, beta] = egcitest(Y, 'creg', 'nc');
       % If the series do co-integrate:
       if h == 1
          i_id(k) = i(j);
          j_id(k) = j;
          betas(k) = beta.coeff;
          spread = norm_prices(:,i(j)) - beta.coeff .* norm_prices(:,j); %check
          spreads(k,:) = [mean(spread), std(spread)];          
          k = k +1;
       end
    end   
       
    Pairs_id = [Pairs_id;{[i_id,j_id]}];
    Pairs_Betas_j = [Pairs_Betas_j; {betas}];
    Pairs_Spreads = [Pairs_Spreads; {spreads}];  
end
fprintf('Done!\n');
clear A b beta betas Dif i i_id i_id_ssd j j_id j_id_ssd k m h norm_prices s_id spread spreads SSDs t_training_months training_date_low training_date_low_id training_date_high training_date_high_id Y

%% Trading Pairs

% Quick adjustment for the index forward:
t_month_first = [t_month_first; numel(dates)+1];

Profits = [];
for m = 13:(numel(t_month_first)-6)
    if mod(m-12, 10) == 0
    fprintf('Processing %d of %d test sets...\n', m-12, numel(t_month_first)-18);
    end
    
    % First getting the index of the stocks that will be considered in this test set:
    s_id = Stocks_id{m-12};
    i_id = Pairs_id{m-12}(:,1);
    j_id = Pairs_id{m-12}(:,2);
    % Getting the info of the pairs:
    betas = Pairs_Betas_j{m-12};
    spreads = Pairs_Spreads{m-12};
    
    % Restraining the test set:
    t_test_months = t_month_first(m:m+6);
    test_date_low = dates(t_test_months(1)-1); 
    test_date_high = dates(t_test_months(end)-1);
    % Getting the indices for the dates of the training set:
    test_date_low_id = dates_LT(test_date_low - minDate + 1);
    test_date_high_id = dates_LT(test_date_high - minDate + 1);
    
    % Getting the prices from the test set (required for trading):
    test_prices = prices(test_date_low_id:test_date_high_id, s_id);
    
    % Building the matrix of normalized prices for this test set:
    norm_prices = prices(test_date_low_id:test_date_high_id, s_id)./prices(test_date_low_id, s_id);
    
    % Let's trade these Pairs!!!
    PR = []; 
    for pair=1:P
        % Setting initial positions to zero:       
        position_i = 0;
        position_j = 0;
        pr_ij = [];
        for t=1:size(norm_prices,1)
            % First calculating the spread (normalized) between the two stocks of the pair:
            spread = norm_prices(t,i_id(pair)) - betas(pair) * norm_prices(t,j_id(pair));
            spread_norm = (spread - spreads(pair, 1))/spreads(pair, 2);
            % If we are not trading the pair in 't':
            if position_i == 0
                % Following the rule of Rad et al. (2016):
                if spread_norm < -2 || spread_norm > 2  
                    if spread > 0 
                        position_i = -1;
                        position_j = betas(pair);
                    else
                        position_i = 1;
                        position_j = -betas(pair);                    
                    end
                    pair_t0 = t;
                    spread_0 = spread;
                end
           % If we are currently trading the pair:
            else
                if isnan(test_prices(t,i_id(pair))) || isnan(test_prices(t,j_id(pair)))
                    position_i = 0;
                    position_j = 0;
                else    
                    profit_ij = position_i*((test_prices(t,i_id(pair)) - test_prices(t-1,i_id(pair)))/test_prices(t-1,i_id(pair))) + position_j*((test_prices(t,j_id(pair)) - test_prices(t-1,j_id(pair)))/test_prices(t-1,j_id(pair)));
                    pr_ij =[pr_ij; profit_ij, t-pair_t0, t-1+test_date_low_id];
                    if sign(spread) ~= sign(spread_0)
                        position_i = 0;
                        position_j = 0;
                    end
                end
            end
        end
        
        PR = [PR; {pr_ij}];

    end
    
    Profits = [Profits; {PR}];

end
clear i_id j_id norm_prices pair pair_t0 position_i position_j PR pr_ij profit_ij s_id Sigma_ij spread spread_norm spread_0 t t_test_months test_date_high test_date_high_id test_date_low test_date_low_id test_prices
save('Results_Cointegration_5', '-v7.3')
fprintf('Done!\n');


