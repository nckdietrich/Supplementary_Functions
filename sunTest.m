%% Sun.m class tests
clc; clear; close all;

% dts = datetime(2021, 3, 24, 0:24, 0, 0); 
% dts = datetime(2020, 1, 1, 0:24, 0, 0);
dts = datetime(2008, 6, 15, 0:24, 0, 0);
jds = juliandate(dts);

% Position of the sun in inertial frame (does not depend on earth's rotation)
[ra_rad,dec_rad] = Sun.solar_position_almanac(jds);

% The local solar time and solar zenith angle are representations of solar position that are relative to a particlar location
% glat = 40.01;
glat = 0;
glon = 0;
% glon = -105.27 + 360;
% glon = -105.27; %Boulder, Colorado (geospacepy-lite uses west-longitudes-are-negative convention)
% glon = 100;

lsts = Sun.local_mean_solar_time(jds,glon);
szas_rad = Sun.solar_zenith_angle(jds,glat,glon);

% Geospacepy-lite functions return angles in radians (unless angle is a latitude or longitude [degrees],
% or a time, that can also be understood as an angle, e.g. local solar time [hours])
figure()
subplot(4,1,1);
plot(dts, rad2deg(ra_rad), 'b.','DisplayName','Solar Right Ascension [deg]');
legend() 
subplot(4,1,2);
plot(dts, rad2deg(dec_rad), 'r.','DisplayName', 'Solar Declination [deg]');
legend()
subplot(4,1,3);
plot(dts, lsts/pi*12,'g.', 'DisplayName','Local Solar Time [hr]');
legend()
subplot(4,1,4);
plot(dts, rad2deg(szas_rad),'k.','DisplayName','Solar Zenith Angle [deg]');
hold on;
yline(90);
legend()


% theta_GST = Sun.greenwich_mean_siderial_time(jds);
lhas = Sun.local_hour_angle(jds, glon);
% figure()
% subplot(2,1,1);
% plot(dts, rad2deg(theta_GST),'b.','DisplayName','Greenwhich mean Sideral Time');
% legend
% subplot(2,1,2);
figure()
plot(dts, lhas*(12/pi) ,'r.','DisplayName','Local Hour Angle');
legend

% glons = -180:2.5:177.5;
% lhas_lons1 = NaN(length(glons),1);
% for i = 1:length(glons)
%     lhas_lons1(i) = Sun.local_hour_angle(jds(2), glons(i));
% end

% lhas_lons2 = Sun.local_hour_angle(jds(2),glons);
% figure()
% plot(glons, lhas_lons2, 'b.', 'DisplayName','Longitudinal Local Hour Angle (loop)');
% hold on;
% plot(glons, lhas_lons1, 'r.', 'DisplayName','Longitudinal Local Hour Angle (one go)');
% plot(glon, lhas(2), 'k.', 'DisplayName', 'One point');
% legend()