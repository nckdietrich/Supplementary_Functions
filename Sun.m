<<<<<<< HEAD
classdef Sun
    %{ 
    (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group 
    Written by Liam M. Kilcommons
    Translated to matlab code by Nick Dietrich
    lil test edit
    %}
    methods (Static)
        function [jd] = datetime2jd(dt)
            %{
            Converts between Python datetime and Julian Date
            (days since 12:00 PM on January 1, 4713 B.C.)
            Implementation is valid from 1900-2199
            
            Parameters
            ----------
            dt : datetime.datetime
            
            Returns
            -------
            jd : float
            %}
            if year(dt) < 1900
                raise ValueError('Year must be 4 digit year')
            end
            t1 = 367.*year(dt);
            t2 = floor(7.*(year(dt)+floor((month(dt)+9.)/12.))/4.);
            t3 = floor(275.*month(dt)./9);
            t4 = day(dt) + 1721013.5;
            t5 = ((second(dt)./60 + minute(dt))./60 + hour(dt))./24;
            jd = t1-t2+t3+t4+t5;
        end
        
        function [alpha_r,delta_r] = solar_position_almanac(jds)
            %{
            Finds the apparent solar right ascension and
            declination for any number of julian dates.
            
            This algorithm is entitled 'Low precision formulas for the Sun'
            and can be found in section C of the Naval Research Laboratory
            Astronomical Almanac (2019). This formula should also be in
            other editions of the Almanac.
            
            The Almanac describes this formula as yeilding a precision
            better than 1' (1/60 degrees) for the years 1950 to 2050
            
            Parameters
            ----------
            
            jds : array or scalar
                Julian dates for which to calculate solar positions
             
            Returns
            -------
            
            alpha : array or scalar (matches input)
                Solar apparent right ascension (angle in equatorial
                plane measured clockwise from the vernal equinox direction),
                in radians
            
            delta : array or scalar (matches input)
                Solar declination (equiv. to subsolar latitude), in radians
            %}
                        
            dt_j2000 = datetime(2000,1,1,12,0,0);
            jd_j2000_epoch = Sun.datetime2jd(dt_j2000);
            jd2000 = jds - jd_j2000_epoch; %J2000 epoch = 2451545.0
            
            %Solar mean longitude (degrees)
            L = 280.460 + .9856474.*jd2000;
            L = mod(L,360.);
            
            %Solar mean anomaly (degrees)
            g = 357.528+0.9856003.*jd2000;
            g = mod(g,360.);
            
            %Solar ecliptic longitude (degrees)
            g_rad = deg2rad(g);
            lam = L + 1.915.*sin(g_rad)+.020.*sin(2*g_rad);
            
            %Solar ecliptic latitude (degrees)
            beta = 0.;
            
            %Obliquity of the ecliptic (degrees)
            epsilon = 23.439 - .0000004.*jd2000;
            
            %Right Ascension (there are 2 formulae...prefer the second b/c arctan)
            %Similar to subsolar longitude but measured counterclockwise from
            %vernal equinox direction instead of counterclockwise from
            %prime meridian
            epsilon_r = deg2rad(epsilon);
            lam_r = deg2rad(lam);
            t = tan(epsilon_r/2).^2;
            f = 180/pi;
            lam_r = deg2rad(lam);
            alpha = lam - f*t.*sin(2.*lam_r) + (f/2.)*t.^2.*sin(4*lam_r);
            %alpha = np.arctan2(np.cos(epsilon_r)*np.sin(lam_r),np.cos(lam_r))
            alpha_r = deg2rad(alpha);
            
            %Declination (equiv. to subsolar latitude, always < 90.)
            delta = rad2deg(asin(sin(epsilon_r).*sin(lam_r)));
            delta_r = deg2rad(delta);
        end
        
        function [gst,sdec,sransn] = solar_position_russell(dt)
            %{
            This function is DEPRECATED use solar_position_almanac instead.
            There is reliable documentation for the solar_position_almanac
            algorithm, whereas the Russell algorithm only references
            'private communication'

            The following is the fortran code from which this code was translated:
            From C.T. Russell, (1971) "Geophysical Coordinate Transformations",
            Cosmic. Electrodyn. 2, 184-196
            ...
            G.D. Mead (private communication) has written a simple subroutine to\
            calculate the position of the Sun in GEI coordinates. It is accurate
            for years 1901 through 2099, to within 0.006 deg. The input is the
            year, day of year and seconds of the day in UT. The output is
            Greenwich Mean Sideral Time in degrees, the ecliptic longitude,
            apparent right ascension and declination of the Sun in degrees.
            The listing of this program follows. We note that the cartesian
            coordinates of the vector from the Earth to the Sun are:
              X = cos(SRASN) cos(SDEC)
              Y = sin(SRASN) cos(SDEC)
              Z = sin(SDEC)
              SUBROUTINE SUN(IYR, IDAY, SECS, GST, SLONG, SRASN, SDEC)
            C PROGRAM TO CALCULATE SIDEREAL TIME AND POSITION OF THE SUN.
            C GOOD FOR YEARS 1901 THROUGH 2099. ACCURACY 0.006 DEGREE.
            C INPUT IS IYR, IDAY (INTEGERS), AND SECS, DEFINING UN. TIME.
            C OUTPUT IS GREENWICH MEAN SIDEREAL TIME (GST) IN DEGREES,
            C LONGITUDE ALONG ECLIPTIC (SLONG), AND APPARENT RIGHT ASCENSION
            C AND DECLINATION (SRASN, SDEC) OF THE SUN, ALL IN DEGREES
              DATA RAD /57.29578/
              DOUBLE PRECISION DJ, FDAY
              IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
              FDAY = SECS/86400
              DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5D0
              T = DJ / 36525
              VL = DMOD (279.696678 + 0.9856473354*DJ, 360.D0)
              GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.D0)
              G = DMOD (358.475845 + 0.985600267*DJ, 360.D0) / RAD
              SLONG = VL + (1.91946 -0.004789*T)*SIN(G) + 0.020094*SIN (2.*G)
              OBLIQ = (23.45229 -0.0130125*T) / RAD
              SLP = (SLONG -0.005686) / RAD
              SIND = SIN (OBLIQ)*SIN (SLP)
              COSD = SQRT(1.-SIND**2)
              SDEC = RAD * ATAN (SIND/COSD)
              SRASN = 180. -RAD*ATAN2
              (COTAN (OBLIQ)*SIND/COSD, -COS (SLP)/COSD)
              RETURN
              END
            
            %}
            iyear = year(dt);
            iday = day(dt, 'dayofyear');
            secs = hour(dt).*3600+ minute(dt).*60+ second(dt);
            fday = secs/86400;
            dj = 365.*(iyear-1900)+(iyear-1901)./4 + iday + fday - 0.5;
            t = dj./36525;
            vl = mod(279.696678 + 0.9856473354.*dj, 360);
            gst = mod(279.690983 + 0.9856473354.*dj + 360.*fday + 180., 360.);
            g = mod(358.475845 + 0.985600267.*dj, 360.) * pi/180;
            slong = vl + (1.91946 -0.004789.*t).*sin(g) + 0.020094.*sin(2.*g);
            obliq = (23.45229 -0.0130125.*t) * pi/180.;

            slp = (slong - 0.005686) .* pi/180;
            sin_d = sin(obliq).*sin(slp);
            cos_d = sqrt(1-sin_d.^2);
            sdec = atan(sin_d./cos_d);
            sransn = pi - atan2(1./tan(obliq).*sin_d./cos_d, -1.*cos(slp)./cos_d);
            %GST is in degrees    
            gst = deg2rad(gst);
        end
        
        function [theta_GST] = greenwich_mean_siderial_time(jds)
            %{
            Calculate the angle in the plane of the equator
            between the vernal equinox direction and the prime meridian (the 
            line of longitude through Greenwich, England).

            Parameters
            ----------
            jds : array or scalar
                The julian date(s) of the times for which the GMST should
                be calculated
            
            Returns
            -------
            theta_GST : array or scalar
                The Greenwich Mean Siderial Time in radians
            
            Notes
            -----
            Because this calculation depends on the actual exact number of earth
            rotations since the J2000 epoch, the time (julian date) strictly speaking
            should be in the UT1 system (the time system determined from observations
            of distant stars), because this system takes into account the small changes
            in earth's rotation speed.

            Generally though, UTC is available instead of UT1. UTC is determined
            from atomic clocks, and is kept within +- 1 second of UT1 
            by the periodic insertion of leap seconds.
            %}
            
            dt_j2000 = datetime(2000,1,1,12,0,0);
            jd_j2000 = Sun.datetime2jd(dt_j2000);
            
            t_ut1 = (jds-jd_j2000)./36525; %Get Julian centuries since the j2000.0 epoch
            %Note that this formula can be broken up into a two part (hours and seconds) version using a two part
            %T_UT1. Where 876600 is multiplied by 3600., and in the exponentiation, the accuracy can be increased
            %by breaking up the T_UT1
            theta_GST_s = 67310.54841+(876600*3600+8640184.812866).*t_ut1+0.093104.*t_ut1.^2-6.2e-6.*t_ut1.^3;

            % NOTE: In Python (and Numpy), the output of modulus is 
            % always the same sign as the divisor (360. in the case of angles)
            % so we can treat negative out of bounds the same
            % as positive. Same convention in Matlab.

            %Make sure abs(theta_GST) <= 86400 seconds
            theta_GST_s = mod(theta_GST_s,86400);

            %Convert theta_GST to degrees from seconds
            theta_GST = theta_GST_s./240;

            % Ensure in 0 to 360.
            theta_GST = mod(theta_GST,360);

            % Radians
            theta_GST = theta_GST .* pi / 180;
        end
        
        function [lhas] = local_hour_angle(jds,glons)
            %{
            Finds local hour angle in radians. The sign convention is that of astronomy 
            (positive to the west, meaning angle increases opposite the 
            rotation direction of earth).
            Under this sign convention the hour angle (in units of hours) is the
            time it will be before for the sun is directly overhead at the specified
            location and time.
            
            Parameters
            ----------
            jds : array or scalar
                Time (as julian date)
            glons : array or scalar
                Geographic longitude
            
            Returns
            -------
            lhas : array or scalar
                Local hour angles for locations at specified times (in radians)

            .. note::
                See also Vallado Figure 3-9 (pp.157)
            %}
            
            [sra,~] = Sun.solar_position_almanac(jds);
            gmst = Sun.greenwich_mean_siderial_time(jds);

            %Greenwich Mean Sideral Time, Right Ascension, and longitude are
            %all measured with positive angles counterclockwise about the north pole
            %(with the earth's rotation direction).

            %But the hour angle is defined as positive opposite the earth's rotation
            phi = deg2rad(glons);
            lhas = (gmst+phi) - sra;
            lhas = wrapToPi(lhas);
%             lhas = mod(lhas + pi, 2*pi) - pi;
        end
        
        function [lmsts] = local_mean_solar_time(jds,glons)
            %{
            Find the local solar time (using the mean equinox)
            
            Parameters
            ----------
            jds : Geographic longitude
                Time (as julian date)
            glons : Geographic longitude
                Geographic longitude
            
            Returns
            -------
            lmsts : Geographic longitude
                    Local mean solar time at locations at specified times (in radians)
            
            Notes
            -----
            As an angle, the solar local time increases in the same direction
            as the ISO 6709 longitude convention (eastward positive).
            The differences between hour angle and local time is that hour angle
            is:
            1.) Measured with positive in the westward direction
            2.) Zero when the sun is directly overhead (instead of 12 for solar time)
            See also Vallado pp. 184
            %}
            
            lhas = Sun.local_hour_angle(jds,glons);
            % lhas = -1*lhas; %Convert to positive in the eastward direction
            lmsts = lhas + pi;  %Equiv to + 12 in hours units
%             lmsts( lmsts < 0 ) = lmsts( lmsts < 0) + 2*pi; % Wrap around to be 0 - 24, 12 directly overhead
        end
        
        function [szas] = solar_zenith_angle(jds,glats,glons)
            %{
            Finds solar zenith angle using Astronomical Almanac low-accuracy
            solar position.
            sza < 90 - dayside
            sza > 90 - nightside
            Parameters
            ----------
            jds : array or scalar
                Time (as julian date)
            glats : array or scalar
                Geographic (geocentric-spherical) latitude of the location
            glons : array or scalar
                Geographic longitude of the location
            Returns
            -------
            szas : array or scalar
                Solar Zenith angles for the time/location combinations
                specified (in radians)
            %}
            lam = deg2rad(glats);
            phi = deg2rad(glons);

            [sra,sdec] = Sun.solar_position_almanac(jds);
            sha = Sun.local_hour_angle(jds,glons);

            cossza = sin(lam).*sin(sdec) + cos(lam).*cos(sdec).*cos(sha);
            szas = acos(cossza);
        end
    end
end
=======
classdef Sun
    %{ 
    (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group 
    Written by Liam M. Kilcommons
    Translated to matlab code by Nick Dietrich
    %}
    methods (Static)
        function [jd] = datetime2jd(dt)
            %{
            Converts between Python datetime and Julian Date
            (days since 12:00 PM on January 1, 4713 B.C.)
            Implementation is valid from 1900-2199
            
            Parameters
            ----------
            dt : datetime.datetime
            
            Returns
            -------
            jd : float
            %}
            if year(dt) < 1900
                raise ValueError('Year must be 4 digit year')
            end
            t1 = 367.*year(dt);
            t2 = floor(7.*(year(dt)+floor((month(dt)+9.)/12.))/4.);
            t3 = floor(275.*month(dt)./9);
            t4 = day(dt) + 1721013.5;
            t5 = ((second(dt)./60 + minute(dt))./60 + hour(dt))./24;
            jd = t1-t2+t3+t4+t5;
        end
        
        function [alpha_r,delta_r] = solar_position_almanac(jds)
            %{
            Finds the apparent solar right ascension and
            declination for any number of julian dates.
            
            This algorithm is entitled 'Low precision formulas for the Sun'
            and can be found in section C of the Naval Research Laboratory
            Astronomical Almanac (2019). This formula should also be in
            other editions of the Almanac.
            
            The Almanac describes this formula as yeilding a precision
            better than 1' (1/60 degrees) for the years 1950 to 2050
            
            Parameters
            ----------
            
            jds : array or scalar
                Julian dates for which to calculate solar positions
             
            Returns
            -------
            
            alpha : array or scalar (matches input)
                Solar apparent right ascension (angle in equatorial
                plane measured clockwise from the vernal equinox direction),
                in radians
            
            delta : array or scalar (matches input)
                Solar declination (equiv. to subsolar latitude), in radians
            %}
                        
            dt_j2000 = datetime(2000,1,1,12,0,0);
            jd_j2000_epoch = Sun.datetime2jd(dt_j2000);
            jd2000 = jds - jd_j2000_epoch; %J2000 epoch = 2451545.0
            
            %Solar mean longitude (degrees)
            L = 280.460 + .9856474.*jd2000;
            L = mod(L,360.);
            
            %Solar mean anomaly (degrees)
            g = 357.528+0.9856003.*jd2000;
            g = mod(g,360.);
            
            %Solar ecliptic longitude (degrees)
            g_rad = deg2rad(g);
            lam = L + 1.915.*sin(g_rad)+.020.*sin(2*g_rad);
            
            %Solar ecliptic latitude (degrees)
            beta = 0.;
            
            %Obliquity of the ecliptic (degrees)
            epsilon = 23.439 - .0000004.*jd2000;
            
            %Right Ascension (there are 2 formulae...prefer the second b/c arctan)
            %Similar to subsolar longitude but measured counterclockwise from
            %vernal equinox direction instead of counterclockwise from
            %prime meridian
            epsilon_r = deg2rad(epsilon);
            lam_r = deg2rad(lam);
            t = tan(epsilon_r/2).^2;
            f = 180/pi;
            lam_r = deg2rad(lam);
            alpha = lam - f*t.*sin(2.*lam_r) + (f/2.)*t.^2.*sin(4*lam_r);
            %alpha = np.arctan2(np.cos(epsilon_r)*np.sin(lam_r),np.cos(lam_r))
            alpha_r = deg2rad(alpha);
            
            %Declination (equiv. to subsolar latitude, always < 90.)
            delta = rad2deg(asin(sin(epsilon_r).*sin(lam_r)));
            delta_r = deg2rad(delta);
        end
        
        function [gst,sdec,sransn] = solar_position_russell(dt)
            %{
            This function is DEPRECATED use solar_position_almanac instead.
            There is reliable documentation for the solar_position_almanac
            algorithm, whereas the Russell algorithm only references
            'private communication'

            The following is the fortran code from which this code was translated:
            From C.T. Russell, (1971) "Geophysical Coordinate Transformations",
            Cosmic. Electrodyn. 2, 184-196
            ...
            G.D. Mead (private communication) has written a simple subroutine to\
            calculate the position of the Sun in GEI coordinates. It is accurate
            for years 1901 through 2099, to within 0.006 deg. The input is the
            year, day of year and seconds of the day in UT. The output is
            Greenwich Mean Sideral Time in degrees, the ecliptic longitude,
            apparent right ascension and declination of the Sun in degrees.
            The listing of this program follows. We note that the cartesian
            coordinates of the vector from the Earth to the Sun are:
              X = cos(SRASN) cos(SDEC)
              Y = sin(SRASN) cos(SDEC)
              Z = sin(SDEC)
              SUBROUTINE SUN(IYR, IDAY, SECS, GST, SLONG, SRASN, SDEC)
            C PROGRAM TO CALCULATE SIDEREAL TIME AND POSITION OF THE SUN.
            C GOOD FOR YEARS 1901 THROUGH 2099. ACCURACY 0.006 DEGREE.
            C INPUT IS IYR, IDAY (INTEGERS), AND SECS, DEFINING UN. TIME.
            C OUTPUT IS GREENWICH MEAN SIDEREAL TIME (GST) IN DEGREES,
            C LONGITUDE ALONG ECLIPTIC (SLONG), AND APPARENT RIGHT ASCENSION
            C AND DECLINATION (SRASN, SDEC) OF THE SUN, ALL IN DEGREES
              DATA RAD /57.29578/
              DOUBLE PRECISION DJ, FDAY
              IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
              FDAY = SECS/86400
              DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5D0
              T = DJ / 36525
              VL = DMOD (279.696678 + 0.9856473354*DJ, 360.D0)
              GST = DMOD (279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.D0)
              G = DMOD (358.475845 + 0.985600267*DJ, 360.D0) / RAD
              SLONG = VL + (1.91946 -0.004789*T)*SIN(G) + 0.020094*SIN (2.*G)
              OBLIQ = (23.45229 -0.0130125*T) / RAD
              SLP = (SLONG -0.005686) / RAD
              SIND = SIN (OBLIQ)*SIN (SLP)
              COSD = SQRT(1.-SIND**2)
              SDEC = RAD * ATAN (SIND/COSD)
              SRASN = 180. -RAD*ATAN2
              (COTAN (OBLIQ)*SIND/COSD, -COS (SLP)/COSD)
              RETURN
              END
            
            %}
            iyear = year(dt);
            iday = day(dt, 'dayofyear');
            secs = hour(dt).*3600+ minute(dt).*60+ second(dt);
            fday = secs/86400;
            dj = 365.*(iyear-1900)+(iyear-1901)./4 + iday + fday - 0.5;
            t = dj./36525;
            vl = mod(279.696678 + 0.9856473354.*dj, 360);
            gst = mod(279.690983 + 0.9856473354.*dj + 360.*fday + 180., 360.);
            g = mod(358.475845 + 0.985600267.*dj, 360.) * pi/180;
            slong = vl + (1.91946 -0.004789.*t).*sin(g) + 0.020094.*sin(2.*g);
            obliq = (23.45229 -0.0130125.*t) * pi/180.;

            slp = (slong - 0.005686) .* pi/180;
            sin_d = sin(obliq).*sin(slp);
            cos_d = sqrt(1-sin_d.^2);
            sdec = atan(sin_d./cos_d);
            sransn = pi - atan2(1./tan(obliq).*sin_d./cos_d, -1.*cos(slp)./cos_d);
            %GST is in degrees    
            gst = deg2rad(gst);
        end
        
        function [theta_GST] = greenwich_mean_siderial_time(jds)
            %{
            Calculate the angle in the plane of the equator
            between the vernal equinox direction and the prime meridian (the 
            line of longitude through Greenwich, England).

            Parameters
            ----------
            jds : array or scalar
                The julian date(s) of the times for which the GMST should
                be calculated
            
            Returns
            -------
            theta_GST : array or scalar
                The Greenwich Mean Siderial Time in radians
            
            Notes
            -----
            Because this calculation depends on the actual exact number of earth
            rotations since the J2000 epoch, the time (julian date) strictly speaking
            should be in the UT1 system (the time system determined from observations
            of distant stars), because this system takes into account the small changes
            in earth's rotation speed.

            Generally though, UTC is available instead of UT1. UTC is determined
            from atomic clocks, and is kept within +- 1 second of UT1 
            by the periodic insertion of leap seconds.
            %}
            
            dt_j2000 = datetime(2000,1,1,12,0,0);
            jd_j2000 = Sun.datetime2jd(dt_j2000);
            
            t_ut1 = (jds-jd_j2000)./36525; %Get Julian centuries since the j2000.0 epoch
            %Note that this formula can be broken up into a two part (hours and seconds) version using a two part
            %T_UT1. Where 876600 is multiplied by 3600., and in the exponentiation, the accuracy can be increased
            %by breaking up the T_UT1
            theta_GST_s = 67310.54841+(876600*3600+8640184.812866).*t_ut1+0.093104.*t_ut1.^2-6.2e-6.*t_ut1.^3;

            % NOTE: In Python (and Numpy), the output of modulus is 
            % always the same sign as the divisor (360. in the case of angles)
            % so we can treat negative out of bounds the same
            % as positive. Same convention in Matlab.

            %Make sure abs(theta_GST) <= 86400 seconds
            theta_GST_s = mod(theta_GST_s,86400);

            %Convert theta_GST to degrees from seconds
            theta_GST = theta_GST_s./240;

            % Ensure in 0 to 360.
            theta_GST = mod(theta_GST,360);

            % Radians
            theta_GST = theta_GST .* pi / 180;
        end
        
        function [lhas] = local_hour_angle(jds,glons)
            %{
            Finds local hour angle in radians. The sign convention is that of astronomy 
            (positive to the west, meaning angle increases opposite the 
            rotation direction of earth).
            Under this sign convention the hour angle (in units of hours) is the
            time it will be before for the sun is directly overhead at the specified
            location and time.
            
            Parameters
            ----------
            jds : array or scalar
                Time (as julian date)
            glons : array or scalar
                Geographic longitude
            
            Returns
            -------
            lhas : array or scalar
                Local hour angles for locations at specified times (in radians)

            .. note::
                See also Vallado Figure 3-9 (pp.157)
            %}
            
            [sra,~] = Sun.solar_position_almanac(jds);
            gmst = Sun.greenwich_mean_siderial_time(jds);

            %Greenwich Mean Sideral Time, Right Ascension, and longitude are
            %all measured with positive angles counterclockwise about the north pole
            %(with the earth's rotation direction).

            %But the hour angle is defined as positive opposite the earth's rotation
            phi = deg2rad(glons);
            lhas = (gmst+phi) - sra;
            lhas = wrapToPi(lhas);
%             lhas = mod(lhas + pi, 2*pi) - pi;
        end
        
        function [lmsts] = local_mean_solar_time(jds,glons)
            %{
            Find the local solar time (using the mean equinox)
            
            Parameters
            ----------
            jds : Geographic longitude
                Time (as julian date)
            glons : Geographic longitude
                Geographic longitude
            
            Returns
            -------
            lmsts : Geographic longitude
                    Local mean solar time at locations at specified times (in radians)
            
            Notes
            -----
            As an angle, the solar local time increases in the same direction
            as the ISO 6709 longitude convention (eastward positive).
            The differences between hour angle and local time is that hour angle
            is:
            1.) Measured with positive in the westward direction
            2.) Zero when the sun is directly overhead (instead of 12 for solar time)
            See also Vallado pp. 184
            %}
            
            lhas = Sun.local_hour_angle(jds,glons);
            % lhas = -1*lhas; %Convert to positive in the eastward direction
            lmsts = lhas + pi;  %Equiv to + 12 in hours units
%             lmsts( lmsts < 0 ) = lmsts( lmsts < 0) + 2*pi; % Wrap around to be 0 - 24, 12 directly overhead
        end
        
        function [szas] = solar_zenith_angle(jds,glats,glons)
            %{
            Finds solar zenith angle using Astronomical Almanac low-accuracy
            solar position.
            sza < 90 - dayside
            sza > 90 - nightside
            Parameters
            ----------
            jds : array or scalar
                Time (as julian date)
            glats : array or scalar
                Geographic (geocentric-spherical) latitude of the location
            glons : array or scalar
                Geographic longitude of the location
            Returns
            -------
            szas : array or scalar
                Solar Zenith angles for the time/location combinations
                specified (in radians)
            %}
            lam = deg2rad(glats);
            phi = deg2rad(glons);

            [sra,sdec] = Sun.solar_position_almanac(jds);
            sha = Sun.local_hour_angle(jds,glons);

            cossza = sin(lam).*sin(sdec) + cos(lam).*cos(sdec).*cos(sha);
            szas = acos(cossza);
        end
    end
end
>>>>>>> 3ac73b763697ebb99ac0156bbe7d658b0b27fbe6
