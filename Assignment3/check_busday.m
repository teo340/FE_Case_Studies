function dates_bus = check_busday(dates)

% INPUT: dates to be checked
% OUTPUT: business dates, using a backward convention

dates_bus=dates;

for i=1:length(dates)
    if weekday(dates(i)) == 1
        dates_bus(i) = dates(i)-2;
     elseif weekday(dates(i)) == 7
         dates_bus(i) = dates(i)-1;
     end
end

end
