function solar_time =  find_solartime(doy, lla, UTC_hr)
%{ 
This function will input lonlat and UTC time to find the corresponding 
solar time

lla: [lon, lat, alt] - units:[deg, deg, km]

https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
%}

% change 365 to 366 for leap years (I don't care too much)
gamma = (2*pi/365)*(doy - 1 + (UTC_hr - 12)/24);

eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma)...
             - 0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma));
% decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) ...
%            - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) ...
%            - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma);
lon = lla(1);
% 360deg/24hr = 15 deg per timezone
if lon > 180
    lon =  lon - 360;
end
timezone = round(lon/15);
% Local time at location
loc_time = UTC_hr+timezone;
if loc_time < 0
    loc_time = loc_time + 24;
end
% Time offset in minutes
time_offset = eqtime + 4*lon - 60*timezone;
% True solar time (in hours)
solar_time = ((loc_time)*60 + time_offset)/60; % + mn (0-59) + sc/60 (0-59);

end