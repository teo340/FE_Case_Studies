function NPV = compute_NPV(discounts, dates, Initial_Flow, growth)
% INPUT:
%   discounts    : Vector of discount factors obtained from bootstrap
%   dates        : Complete set of dates used in the bootstrap
%   Initial_Flow : Initial monthly cash flow value
%   growth       : Annual Average Growth Rate (AAGR)
%
% OUTPUT:
%   NPV          : Net Present Value computed over a 20-year period

% Set the start date for the cash flows 25-Set-2026
startDate=datetime(2026,9,25);

% Generate monthly cash flow dates for 20 years (20*12 months)
cashFlowDates = startDate+calmonths(0:20*12-1); 

% The assignment did not specify according to which convention we were
% supposed to handle cashFlowDates that are not business days.
% If we were to use the same convention as in the complete set of swaps
% rate, the code would be the following:
% for i=1:length(cashFlowDates)
%     if weekday(cashFlowDates(i)) == 1
%         cashFlowDates(i) = cashFlowDates(i) + 1;
%     elseif weekday(cashFlowDates(i)) == 7
%         cashFlowDates(i) = cashFlowDates(i) + 2;
%     end
% end

% Initialize a vector for yearly cash flows over 20 years; set the first year's flow to Initial_Flow
Flows=ones(1,20);
Flows(1)=Initial_Flow;

% Calculate cash flows for each subsequent year by applying AAGR
for i=2:20
    Flows(i)=Flows(i-1)*(1+growth);
end

% For each cash flow date, interpolate the corresponding discount factor using the bootstrapped curve
B=ones(1,length(cashFlowDates));
for i=1:length(cashFlowDates)
    B(i) = interpolation(discounts, dates, cashFlowDates(i));
end

% Convert the yearly cash flows into a monthly series by repeating each yearly value 12 times
flow_months = repelem(Flows, 12); 

% Compute the Net Present Value
NPV = sum(flow_months.*B);     
end
