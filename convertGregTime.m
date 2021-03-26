function dt = convertGregTime(time_int)
%{
Function to convert observation time to a time to a datetime that can
actually be used.

Input: Time of observation in 'days since 1601-1-1' in gregorian calender
    > Input can be a scalar or array
Output: datetime
    > scalar or array

Author: Nick Dietrich
Version: 3.18.2021
%}

dt = datetime(time_int*24*3600, 'ConvertFrom', 'epochtime', 'Epoch', '1601-1-1');

end