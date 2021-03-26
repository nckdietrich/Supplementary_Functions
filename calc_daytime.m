function daytime_s = calc_daytime(time_cdfepoch)
% function to calculate from the cdfepoch time to the second of the current
% day
    daytime_s = NaN(length(time_cdfepoch), 1);
    for i = 1:length(time_cdfepoch)
        time_datenum = todatenum(time_cdfepoch{i});
        time_datetime = datetime(time_datenum, 'ConvertFrom', 'datenum');
        daytime_s(i) = 3600*hour(time_datetime) + 60*minute(time_datetime) + second(time_datetime);
    end
end