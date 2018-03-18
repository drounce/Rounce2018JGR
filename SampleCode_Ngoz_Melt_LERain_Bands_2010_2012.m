%% OVERVIEW: DEBRIS THICKNESS ESTIMATES BASED ON DEM DIFFERENCING BANDS
% Objective: derive debris thickness using an iterative approach that
% solves for when the modeled melt rate agrees with the DEM differencing

% Notes for running the model:
%   The model uses relative paths to the input and output.  This was done
%   to make it easier to share the model and not have to worry about
%   changing the filenames.  That said, the script needs to be run in
%   Matlab with the current folder open to the location of the script.

% Assumptions:
% - Slope and aspect of pixel are not changing through time
% - Surface lowering only considers flux divergence and surface mass
%   balance
% - Climatic mass balance between Oct 15 - May 15 is zero, i.e., any snow
%   that has accumulated has melted

% Limitations:
% - Does not account for shading from surrounding terrain
% - Does not account for snowfall thereby potentially overestimating melt
%   which would underestimate debris thickness
% - Does not account for changes in debris thickness over time 

% Files required to run model successfully:
% - DEM (raster)
% - Slope (raster)
% - Aspect (raster)
% - Elevation change (raster)
% - Elevation change (.csv each box with uncertainty)
% - Glacier outline (raster - 1/0 for glacier/non-glacier pixels)
% - Bands/boxes outlines (raster - each pixel is assigned value of box,
%                         0 represents pixel outside of box)
% - Emergence velocity (.csv each box with uncertainty)
% - Meteorological Data (.csv)
    
% Other important input data
% - Lat/Long starting point (center of upper left pixel)
% - Lat/Long pixel resolution in degrees
% - Number of Monte Carlo simulations
% - Time step (delta_t)
% - Number of years modeled
% - Number of iterations (will affect thickness resolution)

% Monte Carlo simulation allows uncertainty to be incorporated into model
% performance.  The debris properties that are considered are:
% (1) albedo, (2) surface roughness, and (3) thermal conductivity.
% Additional uncertainty regarding products are:
% (4) error associated with elevation change

% Turn on timer
    tic
    
%% INPUT DATA:
% Number of Monte Carlo Simulations
    File_Name = 'Ngoz_600mBands_1000MC_includeponds';
    MC_Simulations_Total = 1000;
    MC_extra = 0; % # of sims already done to ensure count stays correct
    MC_Simulations_emergence = 1000; % # of sims for emergence velocity
    Year_Start = 2010;
    Year_End = 2012;
    N_Years_ModeledMelt = 3; 
    %  Required for emergence velocity
    
% Option ponds and cliffs included (0) or excluded (1)
    option_pondscliffs_excluded = 0;
    %  0 - no, ponds/cliffs included (use column 2)
    %  1 - yes, ponds/cliffs excluded (use column 3)
% Surface lowering and uncertainty for each band
    DEM_Diff_w_uncertainty_Table = csvread([pwd,'/input/Ngoz_SurfaceLowering_Bands.csv']);
    %  col 1 - band number
    %  col 2 = mean surface lowering for each band (ponds/cliffs included)
    %  col 3 = mean surface lowering for each band (ponds/cliffs excluded)
    %  col 4 = error associated with mean surface lowering
    if option_pondscliffs_excluded == 0
        bands_DEM_Diff_raw = DEM_Diff_w_uncertainty_Table (:,2);
    elseif option_pondscliffs_excluded == 1
        bands_DEM_Diff_raw = DEM_Diff_w_uncertainty_Table (:,3);
    end
    bands_DEM_Diff_raw = -1*bands_DEM_Diff_raw;
    % change to positive value for surface lowering such that it agrees 
    % with the modeled melt rate for each band, which is considered to be 
    % positive
    bands_DEM_Diff_raw_uncertainty = DEM_Diff_w_uncertainty_Table(:,4);

% Important inputs (make sure spatial resolution, rows, and columns agree)
%   Note: NoData values are assumed to have a value of zero
    Elevation = imread([pwd,'/input/tif_files/Ngozumpa/Ngoz_2012DEM_10m.tif']);
        nrows = size(Elevation,1);
        ncols = size(Elevation,2);
    Slope_deg = imread([pwd,'/input/tif_files/Ngozumpa/Ngoz_2012DEM_10m_Slope.tif']);
        Slope_rad = Slope_deg.*(pi/180);
    Aspect_deg = imread([pwd,'/input/tif_files/Ngozumpa/Ngoz_2012DEM_10m_Aspect.tif']);
        Aspect_rad = Aspect_deg.*(pi/180);
    DEM_Diff_raw = imread([pwd,'/input/tif_files/Ngozumpa/Ngoz_2010_2012_dhdt_10m.tif']);
        DEM_Diff_raw = DEM_Diff_raw*-1;
        % Now positive indicates melt as opposed to change in elevation
    glacier = imread([pwd,'/input/tif_files/Ngozumpa/Ngoz_rgi60_raster_10m.tif']); 
    bands = imread([pwd,'/input/tif_files/Ngozumpa/Ngoz_Bands_600m_10m.tif']);
        n_bands = max(max(bands));
    % Filter out off-glacier points and NoData points for bands and 
    %   elevation change, which will be used for average computations of 
    %   elevation to calculate temperature differences due to lapse rates 
    %   in the model (NoData assumed to have a value of zero)
        DEM_Diff_raw(glacier==0) = 0;
        bands(DEM_Diff_raw==0) = 0;   
        
% % Enter Emergence Velocities
    emergence_velocity_sims = csvread([pwd,'/input/Ngoz_EmergenceVelocity_simulations.csv']);
    %  rows are each band (# bands)
    %  columns are a specific simulation for the band (1000 simulations)
%     emergence_velocity_sims = zeros(n_bands,MC_Simulations_Total);
%     % neglect emergence velocity by making all zeros
    
% Date of start of simulation
    date_start = datestr(now, 'YYYYmmdd');
%     date_start = '20180317';
    % If date_start causes an error, enter the date manually 'YYYYmmdd'
    
% Enter Latitude and Longitude of Top Left Pixel (top left corner of
% pixel as well) such that the lat and long of each point can be
% estimated.
    % Starting latitude, longitude, and pixel size in angular units
        lat_deg_start = 28.06005611;
        long_deg_start = 86.65907611;
        latlong_delta_deg = 0.00009600326515; % for 10 m pixel
    % Compute latitude and longitude of each pixel, which will be
    % needed for the solar corrections
        lat_deg = zeros(nrows,ncols);
        long_deg = zeros(nrows,ncols);
        for r = 1:nrows
            for c = 1:ncols
                lat_deg(r,c) = lat_deg_start - (r-0.5)*latlong_delta_deg;
                    % subtract because going from top to bottom
                long_deg(r,c) = long_deg_start + (c-0.5)*latlong_delta_deg;
                    % add because going from left to right
            end
        end
        lat_rad = lat_deg.*(pi/180); 
        long_rad = long_deg.*(pi/180);
        
%% Generate random variables for Monte Carlo Simulations
    % Albedo - uniform distribution from 0.1 - 0.4
    Albedo_random = 0.3*rand(MC_Simulations_Total,1)+0.1;
    % z0 - uniform distribution from 0.0035 - 0.0600 m
    z0_random = 0.0565*rand(MC_Simulations_Total,1)+0.0035;
    % Thermal Conductivity - uniform distribution from 0.47 - 1.62 W m-1 K-1    
    k_random = 1.15*rand(MC_Simulations_Total,1)+0.47;
    % Elevation change - normal distribution with mean 0 and stdev 1 that
    %   will be adjusted below for the error of each band
    DEM_Diff_random = normrnd(0,1,[MC_Simulations_Total,1]);
    % Sample from the 1000 simulations of emergence velocity
    emergence_random = round(rand(1000,1)*(MC_Simulations_emergence-1))+1;    
        
%% Define constants to be used with Melt Model work
    row_d = 2700;           % Density of debris (kg m-3)
    c_d = 750;              % Specific heat capacity of debris (J kg-1 K-1)
    I0 = 1368;              % Solar constant (W/m2)
    transmissivity=0.75;    % Vertical atmospheric clear-sky transmissivity (assumed)
    emissivity = 0.95;      % Emissivity of debris surface (Nicholson and Benn, 2006)
    P0 = 101325;            % Standard Atmospheric Pressure (Pa) 
    density_air_0 = 1.29;   % Density of air (kg/m3)
    density_water = 1000;   % Density of water (kg/m3)
    density_ice = 900;      % Density of ice (kg/m3) (Nicholson and Benn, 2006)
    LapseRate = 0.0065;     % Temperature Lapse Rate (K/m)
    Kvk = 0.41;             % Von Karman's Constant
    Lv = 2.49*10^6;         %- Latent Heat of Vaporation of water (J/kg)
    Lf = 3.34*10^5;         % Latent Heat of Fusion of Water (J/kg)
    cA = 1010;              % Specific Heat Capacity of Air (J kg-1 K-1)
    cW = 4.179*10^3;        % Specific Heat Capacity of Water (J kg-1 K-1)
    R = 461;                % Gas Constant for water vapor (J Kg-1 K-1)
    
%% Data about Pyramid Station
    Elev_PyrStat = 5035;  %m
    Slope_PyrStat_deg = 0; % assuming Pyramid Station flat on top of ridge or Sky View Factor = 1
    Slope_PyrStat_rad = Slope_PyrStat_deg.*(pi/180);
    Aspect_PyrStat_deg = 0; % deg CW from N
    Aspect_PyrStat_rad = Aspect_PyrStat_deg.*(pi/180);
    P_PyrStat = P0*exp(-0.0289644*9.81*Elev_PyrStat/(8.31447*288.15)); % Pressure at Pyramid Station
    za = 2;                 % m, height of air temperature instrument
    delta_t = 60*60; % s, time step of Pyramid Station data (hourly)

%% Meteorological Data for the entire study period (DEM Differencing)
    Meteorological_data = csvread([pwd,'/input/PyrStat_20100604_20121015_forNgozumpa.csv']);
        % 5 days for model spinup is included
        Year = Meteorological_data(:,1);
        Julian_Day_of_year = Meteorological_data(:,2);    
        Month = Meteorological_data(:,3);
        Day = Meteorological_data(:,4);
        Hour = Meteorological_data(:,5);
        Minute = Meteorological_data(:,6);
            Time = Hour + Minute/60;
            timezone = 5.75;
        % Column 7 is pressure, which is not needed in the model
        Tair_PyrStat = Meteorological_data(:,8) + 273.15; %celsius to kelvin
        RH_PyrStat = Meteorological_data(:,9) / 100; % percentage to decimal
        u_PyrStat_raw = Meteorological_data(:,10);
        Rain_PyrStat = Meteorological_data(:,11)/1000; % mm to meters
        Sin_PyrStat = Meteorological_data(:,12);
        Lin_PyrStat = Meteorological_data(:,13);

    % Table of missing meteorological data (days per month per year)
    %   (required to perform melt corrections)
        Table_MissingDaysPerMonthPerYear = [0,14,28,0,3,3;
                                            2,14,1,0,0,0;
                                            10,15,0,0,0,0];
            % Row 1 = 2010; Row 2 = 2011; Row 3 = 2012
            % Col 1 = May; Col 2 = June, ..., Col 6 = Oct
        Table_DaysPerMonthPerYear = [0, 22, 31, 31, 30, 15;
                                    17, 30, 31, 31, 30, 15;
                                    17, 30, 31, 31, 30, 15];

    % Enter the number of days for model start-up
        days_spinup = 5;
        days_spinup_nsteps = days_spinup*24*60*60/delta_t;
    % Compute the number of time steps
        NumberTimeSteps = length(Sin_PyrStat);
        
%% Average elevation of each band needed to adjust air temperatures
bands_Elev = zeros(n_bands,1);
for n = 1:n_bands
    ind = find(bands == n);
    bands_Elev(n,1) = sum(Elevation(ind))/(sum(sum(bands==n)));
end

%% Perform solar calculations to adjust for the position of the sun (see NOAA_Solar_Calcs_forDEMDiff script)
    Year_NOAA = zeros(NumberTimeSteps,1);
    JulianDay_NOAA = zeros(NumberTimeSteps,1);
    JulianCentury = zeros(NumberTimeSteps,1);
    GeomMeanLongSun_deg = zeros(NumberTimeSteps,1);
    GeomMeanLongSun_rad = zeros(NumberTimeSteps,1);
    GeomMeanAnomSun_deg = zeros(NumberTimeSteps,1);
    GeomMeanAnomSun_rad = zeros(NumberTimeSteps,1);
    EccentEarthOrbit = zeros(NumberTimeSteps,1);
    SunEqofCtr = zeros(NumberTimeSteps,1);
    SunTrueLong_deg = zeros(NumberTimeSteps,1);
    SunTrueAnom_deg = zeros(NumberTimeSteps,1);
    SunTrueAnom_rad = zeros(NumberTimeSteps,1);
    SunRadVector = zeros(NumberTimeSteps,1);
    SunAppLong_deg = zeros(NumberTimeSteps,1);
    SunAppLong_rad = zeros(NumberTimeSteps,1);
    MeanObliqEcliptic_deg = zeros(NumberTimeSteps,1);
    MeanObliqEcliptic_rad = zeros(NumberTimeSteps,1);
    ObliqCorr_deg = zeros(NumberTimeSteps,1);
    ObliqCorr_rad = zeros(NumberTimeSteps,1);
    SunRtAscen_deg = zeros(NumberTimeSteps,1);
    SunDeclin_deg = zeros(NumberTimeSteps,1);
    SunDeclin_rad = zeros(NumberTimeSteps,1);
    VarY = zeros(NumberTimeSteps,1);
    EqofTime = zeros(NumberTimeSteps,1);
    TrueSolarTime = zeros(NumberTimeSteps,1);
    HourAngle_deg = zeros(NumberTimeSteps,1);
    HourAngle_rad = zeros(NumberTimeSteps,1);
    SolarZenithAngle_deg = zeros(NumberTimeSteps,1);
    SolarZenithAngle_rad = zeros(NumberTimeSteps,1);
    SolarElevationAngle_deg = zeros(NumberTimeSteps,1);
    SolarElevationAngle_rad = zeros(NumberTimeSteps,1);
    ApproxAtmosRefrac_deg = zeros(NumberTimeSteps,1);
    SolarElevationAngleCorr_deg = zeros(NumberTimeSteps,1);
    SolarZenithAngleCorr_deg = zeros(NumberTimeSteps,1);
    SolarZenithAngleCorr_rad = zeros(NumberTimeSteps,1);
    SolarAzimuthAngle_deg = zeros(NumberTimeSteps,1);
    SolarAzimuthAngle_rad = zeros(NumberTimeSteps,1);
    rm_r2 = zeros(NumberTimeSteps,1);

    for i = 1:NumberTimeSteps
        % Julian Day 
            JulianDay_NOAA(i) = floor(365.25*(Year(i)-1900)+1) + Julian_Day_of_year(i) + (Time(i)-timezone)/24 + 2415018.5;
                % +1 accounts for the fact that Day 1 is January 1, 1900
                % 2415018.5 converts from 1900 to NOAA Julian Day of Year
        % Julian Century
            JulianCentury(i) = (JulianDay_NOAA(i)-2451545)/36525;
        % Geom Mean Long Sun (deg)
            GeomMeanLongSun_deg(i) = mod(280.46646+JulianCentury(i).*(36000.76983+JulianCentury(i)*0.0003032),360);
            GeomMeanLongSun_rad(i) = GeomMeanLongSun_deg(i).*(pi/180);
            % MOD fxn in excel returns the remainder after a number is divided by a divisor
        % Geom Mean Anom Sun (deg)
            GeomMeanAnomSun_deg(i) = 357.52911+JulianCentury(i).*(35999.05029-0.0001537*JulianCentury(i));
            GeomMeanAnomSun_rad(i) = GeomMeanAnomSun_deg(i).*(pi/180);
        % Eccent Earth Orbit
            EccentEarthOrbit(i) = 0.016708634-JulianCentury(i).*(0.000042037+0.0000001267*JulianCentury(i));
        % Sun Eq of Ctr
            SunEqofCtr(i) = sin(GeomMeanAnomSun_rad(i)).*(1.914602-JulianCentury(i).*(0.004817+0.000014*JulianCentury(i)))+sin(2*GeomMeanAnomSun_rad(i)).*(0.019993-0.000101*JulianCentury(i))+sin(3*GeomMeanAnomSun_rad(i))*0.000289;
        % Sun True Long (deg)
            SunTrueLong_deg(i) = GeomMeanLongSun_deg(i) + SunEqofCtr(i);
        % Sun True Anom (deg)
            SunTrueAnom_deg(i) = GeomMeanAnomSun_deg(i) + SunEqofCtr(i);
            SunTrueAnom_rad(i) = SunTrueAnom_deg(i).*(pi/180);
        % Sun Rad Vector (AUs)
            SunRadVector(i) = (1.000001018*(1-EccentEarthOrbit(i).*EccentEarthOrbit(i)))./(1+EccentEarthOrbit(i).*cos(SunTrueAnom_rad(i)));
        % Sun App Long (deg)
            SunAppLong_deg(i) = SunTrueLong_deg(i)-0.00569-0.00478*sin((125.04-1934.136*JulianCentury(i)).*(pi/180));
            SunAppLong_rad(i) = SunAppLong_deg(i).*(pi/180);
        % Mean Obliq Ecliptic (deg)
            MeanObliqEcliptic_deg(i) = 23+(26+((21.448-JulianCentury(i).*(46.815+JulianCentury(i).*(0.00059-JulianCentury(i)*0.001813))))/60)/60;
            MeanObliqEcliptic_rad(i) = MeanObliqEcliptic_deg(i).*(pi/180);
        % Obliq Corr (deg)
            ObliqCorr_deg(i) = MeanObliqEcliptic_deg(i)+0.00256*cos((125.04-1934.136*JulianCentury(i)).*(pi/180));
            ObliqCorr_rad(i) = ObliqCorr_deg(i).*(pi/180);
        % Sun Rt Ascen (deg)
            SunRtAscen_deg(i) = 180/pi*atan((cos(ObliqCorr_rad(i)).*sin(SunAppLong_rad(i)))./cos(SunAppLong_rad(i)));
        % Sun Declin (deg)
            SunDeclin_deg(i) = 180/pi*asin(sin(ObliqCorr_rad(i)).*sin(SunAppLong_rad(i)));
            SunDeclin_rad(i) = SunDeclin_deg(i).*(pi/180);
        % VarY
            VarY(i) = tan(ObliqCorr_deg(i)/2.*(pi/180)).*tan(ObliqCorr_deg(i)/2.*(pi/180));
        % Eq of Time (min)
            EqofTime(i) = 4*180/pi*(VarY(i).*sin(2*GeomMeanLongSun_rad(i))-2*EccentEarthOrbit(i).*sin(GeomMeanAnomSun_rad(i))+4*EccentEarthOrbit(i).*VarY(i).*sin(GeomMeanAnomSun_rad(i)).*cos(2*GeomMeanLongSun_rad(i))-0.5*VarY(i).*VarY(i).*sin(4*GeomMeanLongSun_rad(i))-1.25*EccentEarthOrbit(i).*EccentEarthOrbit(i).*sin(2*GeomMeanAnomSun_rad(i)));
        % True Solar Time (min)
            TrueSolarTime(i) = mod((Time(i)*60*1440+Time(i)*60+EqofTime(i)+4*long_deg(r,c)-60*timezone),1440);
        % Hour Angle (deg)
            if TrueSolarTime(i)/4 < 0
                HourAngle_deg(i) = TrueSolarTime(i)/4+180;
            else
                HourAngle_deg(i) = TrueSolarTime(i)/4-180;
            end
                HourAngle_rad(i) = HourAngle_deg(i).*(pi/180);
        % Solar Zenith Angle (deg)
            SolarZenithAngle_deg(i) = 180/pi*acos(sin(lat_rad(r,c)).*sin(SunDeclin_rad(i))+cos(lat_rad(r,c)).*cos(SunDeclin_rad(i)).*cos(HourAngle_rad(i)));
            SolarZenithAngle_rad(i) = SolarZenithAngle_deg(i).*(pi/180);
        % Solar Elevation Angle (deg)
            SolarElevationAngle_deg(i) = 90-SolarZenithAngle_deg(i);
            SolarElevationAngle_rad(i) = SolarElevationAngle_deg(i).*(pi/180);
        % Approx Atmospheric Refraction (deg)  
            if SolarElevationAngle_deg(i) > 85
                ApproxAtmosRefrac_deg(i) = 0;
            elseif SolarElevationAngle_deg(i) > 5
                ApproxAtmosRefrac_deg(i) = 58.1./tan(SolarElevationAngle_rad(i))-0.07./((tan(SolarElevationAngle_rad(i))).^3)+0.000086./((tan(SolarElevationAngle_rad(i))).^5);
            elseif SolarElevationAngle_deg > -0.575
                ApproxAtmosRefrac_deg(i) = 1735+SolarElevationAngle_deg(i).*(-518.2+SolarElevationAngle_deg(i).*(103.4+SolarElevationAngle_deg(i).*(-12.79+SolarElevationAngle_deg(i)*0.711)));
            else
                ApproxAtmosRefrac_deg(i) = -20.772./tan(SolarElevationAngle_rad(i));
            end
                ApproxAtmosRefrac_deg(i) = ApproxAtmosRefrac_deg(i)/3600;
        % Solar Elevation Correct for Atm Refraction (deg)
            SolarElevationAngleCorr_deg(i) = SolarElevationAngle_deg(i) + ApproxAtmosRefrac_deg(i);
        % Solar Zenith Angle Corrected for Atm Refraction (deg)
            SolarZenithAngleCorr_deg(i) = 90 - SolarElevationAngleCorr_deg(i);
            SolarZenithAngleCorr_rad(i) = SolarZenithAngleCorr_deg(i).*(pi/180);
        % Solar Azimuth Angle (deg CW from N)    
            if HourAngle_deg(i) > 0
                SolarAzimuthAngle_deg(i) = ((180/pi*(acos(((sin(lat_rad(r,c)).*cos(SolarZenithAngle_rad(i)))-sin(SunDeclin_rad(i)))./(cos(lat_rad(r,c)).*sin(SolarZenithAngle_rad(i)))))+180)/360-floor((180/pi*(acos(((sin(lat_rad(r,c)).*cos(SolarZenithAngle_rad(i)))-sin(SunDeclin_rad(i)))./(cos(lat_rad(r,c)).*sin(SolarZenithAngle_rad(i)))))+180)/360))*360;
            else
                SolarAzimuthAngle_deg(i) = ((540-180/pi*(acos(((sin(lat_rad(r,c)).*cos(SolarZenithAngle_rad(i)))-sin(SunDeclin_rad(i)))./(cos(lat_rad(r,c)).*sin(SolarZenithAngle_rad(i))))))/360-floor((540-180/pi*(acos(((sin(lat_rad(r,c)).*cos(SolarZenithAngle_rad(i)))-sin(SunDeclin_rad(i)))./(cos(lat_rad(r,c)).*sin(SolarZenithAngle_rad(i))))))/360))*360;
            end  
                SolarAzimuthAngle_rad(i) = SolarAzimuthAngle_deg(i).*(pi/180);
    % Distance from sun based on eccentricity of orbit (dbar/d)^2
        J = 0:364; 
        phi_d = (J.*2*pi)./365; %Julian Date radians (eq. A.5 Hartmann 1994)
        d0 = 1.000110.*cos(0*phi_d);
        d1 = 0.034221.*cos(1*phi_d) + 0.001280.*sin(1.*phi_d);
        d2 = 0.000719.*cos(2*phi_d) + 0.000077.*sin(2.*phi_d);
        d = d0 + d1 + d2; %(dbar/d)^2-one value for each day of year
        rm_r2(i) = d(Julian_Day_of_year(i)); % (rm/r)^2 for specific day of satellite image
    end % end solar calcs corrections
    
%% Compute Solar Corrections to adjust incoming radiation for the position of the sun 
%     Sin_timeseries = csvread([pwd,'/input/Ngoz_Sin_timeseries_bands.csv'])';
    %  Sin for each time step for each band.  This is the output of the
    %  code below, but saves the computational time of running it.

% Note: the following loop needs to be run the first time the model is run.
%       It computes the incoming solar radiation at every pixel for every
%       time step and then averages the incoming solar radiation for each
%       band.  Depending on the resolution of the data, the script can be
%       very slow; hence, it is useful to save the Sin_timeseries as a .csv
%       to skip this step in future runs.

     P = P0*exp(-0.0289644*9.81*Elevation/(8.31447*288.15)); % Pressure (barometric pressure formula)
     theta_PyrStat = zeros(1,NumberTimeSteps);
     I_PyrStat = zeros(1,NumberTimeSteps);
     theta = zeros(1,NumberTimeSteps);
     I = zeros(1,NumberTimeSteps);
     Sin = zeros(1,NumberTimeSteps);
     Sin_timeseries = zeros(n_bands,NumberTimeSteps);
     bands_Sin = zeros(n_bands,1);
     for i = 1:NumberTimeSteps

                % Angle of Incidence b/w normal to grid slope at PyrStat and solar beam 
                     theta_PyrStat(i) = acos(cos(Slope_PyrStat_rad).*cos(SolarZenithAngleCorr_rad(i))+sin(Slope_PyrStat_rad).*sin(SolarZenithAngleCorr_rad(i)).*cos(SolarAzimuthAngle_rad(i)-Aspect_PyrStat_rad));
                % Potential Clear-Sky Solar Radiation at Pyramid Station
                    I_PyrStat(i) = I0*rm_r2(i)*transmissivity.^(P_PyrStat./(P0*cos(SolarZenithAngleCorr_rad(i)))).*cos(theta_PyrStat(i));
                        % Impossible to have negative incoming shortwave radiation
                            if I_PyrStat(i) < 1 % Use less than 1 instead of 0 such that the correction for the effect of shading doesn't multiply to excessively large values
                                I_PyrStat(i) = 0;
                            end
                % Angle of Incidence b/w normal to grid slope and solar beam
                    theta_matrix = acos(cos(Slope_rad).*cos(SolarZenithAngleCorr_rad(i))+sin(Slope_rad).*sin(SolarZenithAngleCorr_rad(i)).*cos(SolarAzimuthAngle_rad(i)-Aspect_rad));
                % Potential Clear-Sky Solar Radiation
                    I_matrix = I0*rm_r2(i)*transmissivity.^(P./(P0*cos(SolarZenithAngleCorr_rad(i)))).*cos(theta_matrix);
                        if SolarZenithAngleCorr_rad(i) > (pi/2)
                            I_matrix = zeros(nrows,ncols);
                        end
                % Effect of Shading (from Hock and Noetzli, 1997)
                        % Not included - see old scripts
                        if I_PyrStat(i) == 0
                            Sin_matrix = zeros(nrows,ncols);
                        else
                            Sin_matrix = I_matrix*Sin_PyrStat(i)/I_PyrStat(i);
                        end  
                        Sin_matrix(Sin_matrix < 0) = 0;
                
            for n = 1:n_bands
                ind = find(bands == n);
                bands_Sin(n,1) = sum(Sin_matrix(ind))/sum(sum(bands==n));
            end
            
            Sin_timeseries(:,i) = bands_Sin;
     end % end loop that does solar corrections
    
%% Start the Parallel Computing Loop
% Note: if not using Parallels, then you can change "parfor" to "for"

% Output data to .txt file
    % Print the header as the other file will not contain this information
    Output_FileName_DebrisThickness = [pwd,'/output/',File_Name,'_',num2str(date_start),'_header.txt'];     
    fileID_header = fopen(Output_FileName_DebrisThickness,'a');
    fprintf(fileID_header,'%1s %1s %3s %5s %3s %4s %3s %6s %6s %6s %6s %6s %4s\r\n','n','MC','d','Melt','alb','z0','k','DEMerr','DEMDif','em_mc','em_vel','TotDif','n_it');
    fclose(fileID_header);
    
% Loop through each elevation band
parfor n = 1:36
    
    % Output data to .txt file
    Output_FileName_DebrisThickness = [pwd,'/output/',File_Name,'_',num2str(date_start),'_band_',num2str(n),'.txt'];     
    fileID = fopen(Output_FileName_DebrisThickness,'a');
        
    parfor MC = 1:MC_Simulations_Total 
%     for MC = 1:MC_Simulations_Total    
    % Note: If not using parallel computing, then change "parfor" to "for"
    % Note: If you get an error within the parfor loop, switch to
    % the "for" loop as this will identify the specific line where
    % the error occured (this doesn't happen in parfor loop).

        Elevation_pixel = bands_Elev(n,1);
        DEM_Diff_raw_pixel = bands_DEM_Diff_raw(n,1);
        DEM_Diff_raw_pixel_uncertainty = bands_DEM_Diff_raw_uncertainty(n,1);
        Sin = Sin_timeseries(n,:);
        
        % Set debris properties
            % Albedo
                Albedo = Albedo_random(MC,1);
            % Surface roughness (m)
                z0 = z0_random(MC,1);
            % Thermal conductivity (W m-1 K-1)
                k = k_random(MC,1);
            % Emergence velocity (m, total time period)
                emergence_velocity = emergence_velocity_sims(n,emergence_random(MC,1))*N_Years_ModeledMelt;
            % DEM difference (m, total time period)
                DEM_Diff_noem = DEM_Diff_raw_pixel + DEM_Diff_random(MC,1)*DEM_Diff_raw_pixel_uncertainty;
                DEM_Diff = DEM_Diff_raw_pixel + DEM_Diff_random(MC,1)*DEM_Diff_raw_pixel_uncertainty + emergence_velocity;
            % Additional properties
                % Dimensionless transfer coefficient for turbulent heat fluxes
                    A_Benn = Kvk^2/(log(za/z0))^2;
                % Adjust wind speed from 5 m to 2 m accounting for surface roughness
                    u_PyrStat = u_PyrStat_raw*(log(2/z0)/(log(5/z0)));

        % Set debris iterations
            debris_thickness_min = 0.02; % m
            debris_thickness_middle = 0.15; % m
            debris_thickness_max = 5; % m
            debris_thickness_round_interval = 0.01;
            N_iterations_Break_lower = 7;
            N_iterations_Break_upper = 13;
                % For 0.03-5 m (0.24 middle), 6 lower & 10 upper gives thickness to the nearest 5 cm
                % For 0.02-5 m (0.15 middle), 7 lower & 13 upper gives thickness to the nearest 1 cm (must change rounding to nearest 0.001)
            N_iterations_Break = max(N_iterations_Break_lower,N_iterations_Break_upper);
            Table_Melt = zeros(NumberTimeSteps,5);
            Julian_Day_of_year_adjustedMelt = Julian_Day_of_year;
            Melt_Modeled_pixel = 0;
            N_iterations = 0;
            
        % Set temporary variables to avoid error messages in parfor loop
            debris_thickness = 0;
            debris_thickness_i1 = 0;
            debris_thickness_i2 = 0;
            Melt_i1 = 0;
            Melt_i2 = 0;
            
        % Set while loop to vary the debris thickness for each pixel
%             while N_iterations <= N_iterations_Break
            while N_iterations <= N_iterations_Break & bands_DEM_Diff_raw(n,1) ~= 0
                if N_iterations == N_iterations_Break + 1; % if the number of iterations exceeds the break number, then you have your solution
                    break
                end
                N_iterations = N_iterations + 1;
                
                if N_iterations == 1 % run middle debris thickness
                    debris_thickness = debris_thickness_middle;
                elseif N_iterations == 2 % run upper or lower bound
                    if DEM_Diff < Melt_Modeled_pixel % debris thickness is too thin (increase debris thickness to get less melt)
                        debris_thickness_i1 = debris_thickness; % record debris thickness (initial lower bound)
                        Melt_i1 = Melt_Modeled_pixel; % record melt associated with this debris thickness
                        debris_thickness = debris_thickness_max; % next iteration run max debris thickness (initial upper bound)
                    else % debris thickness is too large (reduce debris thickness to get more melt)
                        debris_thickness_i2 = debris_thickness; % record debris thickness (initial upper bound)
                        Melt_i2 = Melt_Modeled_pixel; % record melt associated with this debris thickness
                        debris_thickness = debris_thickness_min; % next iteration run min debris thickness (initial lower bound)
                    end
                elseif N_iterations == 3 % iterate to solution or flag for being outside of upper and lower bounds
                    if debris_thickness == debris_thickness_max
                        if DEM_Diff < Melt_Modeled_pixel % if melt is less than max debris thickness, then thickness is greater than max
                            debris_thickness = debris_thickness_max;
                            break
                        end
                        debris_thickness_i2 = debris_thickness; % i1 is the middle; i2 is the max
                        Melt_i2 = Melt_Modeled_pixel;
                        N_iterations_Break = N_iterations_Break_upper;
                        debris_thickness = round((debris_thickness_i1 + debris_thickness_i2)/2*1000)/1000;
%                             % rounding to nearest centimeter is required to get the proper height of pixel layers
                    elseif debris_thickness == debris_thickness_min
                        if DEM_Diff > Melt_Modeled_pixel % if melt exceeds melt of min debris thickness, then thickness is less than min (ponds, ice, etc.)
                            debris_thickness = 0;
                            break
                        end
                        debris_thickness_i1 = debris_thickness; % i1 is the min; i2 is the middle
                        Melt_i1 = Melt_Modeled_pixel;
                        N_iterations_Break = N_iterations_Break_lower;
                        debris_thickness = round((debris_thickness_i1 + debris_thickness_i2)/2*1000)/1000;
%                             % rounding to nearest centimeter is required to get the proper height of pixel layers
                    end
                else % continue iterating until reach iteration break
                    debris_thickness_i3 = debris_thickness;
                    Melt_i3 = Melt_Modeled_pixel;
                    if DEM_Diff > Melt_i3 & DEM_Diff < Melt_i1 % iterate to move closer to target (here: debris thickness decreases)
                        debris_thickness_i2 = debris_thickness_i3; % Update debris thickness bounds for next iteration
                        Melt_i2 = Melt_i3; % Update melt bounds for next iteration
                        debris_thickness = round((debris_thickness_i3 + debris_thickness_i1)/2*1000)/1000;
                    elseif DEM_Diff < Melt_i3 & DEM_Diff > Melt_i2 % iterate to move closer to target (here: debris thickness increases)
                            debris_thickness_i1 = debris_thickness_i3; % Update debris thickness bounds for next iteration
                            Melt_i1 = Melt_i3;% Update melt bounds for next iteration
                            debris_thickness = round((debris_thickness_i3 + debris_thickness_i2)/2*1000)/1000;
                    end
                end % end N_iterations loop that varies debris thickness
                
                % Run the debris-covered energy balance model
                    % Compute height of each internal layer
                        h = debris_thickness / 10;
                    % Compute various information needed for 
                        P = P0*exp(-0.0289644*9.81*Elevation_pixel/(8.31447*288.15)); % Pressure (barometric pressure formula)  
                        C = k*delta_t/(2*row_d*c_d*h^2);    % constant defined by Reid and Brock (2010) for Crank-Nicholson Scheme
                        N = debris_thickness/h + 1;         % Number of internal calculation layers + 1 to include the surface layer
                        Tair = Tair_PyrStat - LapseRate*(Elevation_pixel-Elev_PyrStat);
                        eZ_Saturated = 611*exp(-Lv/R*(1./Tair-1/273.15));
                    
                    % "Crank Nicholson Newton Raphson" Method for LE Rain (No Snow)
                    % Compute Ts from surface energy balance model using Newton-Raphson Method
                    % at each time step and Td at all points in debris layer
                        Td = zeros(N, NumberTimeSteps);
                            % Td is the temperature in the debris with row 1 being Ts
                        a_Crank = zeros(N,NumberTimeSteps);
                        b_Crank = zeros(N,NumberTimeSteps);
                        c_Crank = zeros(N,NumberTimeSteps);
                        d_Crank = zeros(N,NumberTimeSteps);
                        A_Crank = zeros(N,NumberTimeSteps);
                        S_Crank = zeros(N,NumberTimeSteps);
                        n_iterations = zeros(1, NumberTimeSteps);
                        Ts_past = zeros(1, NumberTimeSteps);
                        eS_Saturated = zeros(1, NumberTimeSteps);
                        eS_dry = zeros(1, NumberTimeSteps);
                        eS = zeros(1, NumberTimeSteps);
                        eZ = zeros(1, NumberTimeSteps);
                        LE_Benn = zeros(1, NumberTimeSteps);
                        Rn = zeros(1, NumberTimeSteps);
                        H_Benn = zeros(1, NumberTimeSteps);
                        Qc = zeros(1, NumberTimeSteps);
                        P_Flux = zeros(1, NumberTimeSteps);
                        dLE_Benn = zeros(1, NumberTimeSteps);
                        dRn = zeros(1, NumberTimeSteps);
                        dH_Benn = zeros(1, NumberTimeSteps);
                        dQc = zeros(1, NumberTimeSteps);
                        dP_Flux = zeros(1, NumberTimeSteps);
                        F_Ts_Benn = zeros(1, NumberTimeSteps);
                        dF_Ts_Benn = zeros(1, NumberTimeSteps);
                        Qc_ice = zeros(1, NumberTimeSteps);
                        Melt = zeros(1, NumberTimeSteps);

                        for i = 1:NumberTimeSteps
                            n_iterations(i) = 0;
                            Ts_past(i) = 0;
                            Td(N,i) = 273.15;
                            % Initially assume Ts = Tair, for all other time steps assume it's equal to previous Ts
                                if i == 1
                                    Td(1,i) = Tair(i);
                                else
                                    Td(1,i) = Td(1,i-1);
                                end
                            % Calculate temperature profile in the debris
                                % For t = 0, which is i = 1, assume initial condition of linear temperature profile in the debris
                                if i == 1
                                    Td_gradient = (Td(1,1) - Td(N,1))/debris_thickness;
                                    for j = 2:1:(N-1)
                                        Td(j,1) = Td(1,1) - (j*h)*Td_gradient;
                                    end
                                else
                                    % Perform Crank-Nicholson Scheme
                                        for j = 2:1:(N-1)
                                            % Equations A8 in Reid and Brock (2010) 
                                                a_Crank(j,i) = C;
                                                b_Crank(j,i) = 2*C+1;
                                                c_Crank(j,i) = C;

                                            % Equations A9 in Reid and Brock (2010) 
                                                if j == 2
                                                    d_Crank(j,i) = C*Td(1,i) + C*Td(1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                elseif j < (N-1)
                                                    d_Crank(j,i) = C*Td(j-1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                elseif j == (N-1)
                                                    d_Crank(j,i) = 2*C*Td(N,i) + C*Td(N-2,i-1) + (1-2*C)*Td(N-1,i-1);
                                                end
                                                    % note notation:
                                                        % "i-1" refers to the past
                                                        % "j-1" refers to the cell above it
                                                        % "j+1" refers to the cell below it          

                                            % Equations A10 and A11 in Reid and Brock (2010)
                                                if j == 2
                                                    A_Crank(j,i) = b_Crank(j,i);
                                                    S_Crank(j,i) = d_Crank(j,i);
                                                else
                                                    A_Crank(j,i) = b_Crank(j,i) - a_Crank(j,i)/A_Crank(j-1,i)*c_Crank(j-1,i);
                                                    S_Crank(j,i) = d_Crank(j,i) + a_Crank(j,i)/A_Crank(j-1,i)*S_Crank(j-1,i);
                                                end
                                        end

                                    % Equations A12 in Reid and Brock (2010)
                                        for j = N-1:-1:2
                                            if j == (N-1)
                                                Td(j,i) = S_Crank(j,i)/A_Crank(j,i);
                                            else
                                                Td(j,i) = 1/A_Crank(j,i)*(S_Crank(j,i)+c_Crank(j,i)*Td(j+1,i));
                                            end
                                        end
                                end
                                % Assume snow-free surface and compute fluxes normally
                                    % Calculate Surface Energy Fluxes
                                        if Rain_PyrStat(i) > 0
                                            eS_Saturated(i) = 611*exp(-Lv/R*(1/Td(1,i)-1/273.15));
                                            eS(i) = eS_Saturated(i);
                                            eZ(i) = RH_PyrStat(i)*eZ_Saturated(i);
                                            LE_Benn(i) = 0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*(eZ(i)-eS(i));
                                        else
                                            eS_Saturated(i) = 0;
                                            LE_Benn(i) = 0;
                                        end
                                        Rn(i) = Sin(i)*(1-Albedo) + emissivity*(Lin_PyrStat(i) - (5.67*10^-8*Td(1,i)^4));
                                        H_Benn(i) = density_air_0*(P/P0)*cA*A_Benn*u_PyrStat(i)*(Tair(i)-Td(1,i));
                                        P_Flux(i) = density_water*cW*Rain_PyrStat(i)/(delta_t)*(Tair(i)-Td(1,i));
                                        Qc(i) = k*(Td(2,i) - Td(1,i))/h;
                                        F_Ts_Benn(i) = Rn(i) + H_Benn(i) + LE_Benn(i) + Qc(i) + P_Flux(i);
                                        dRn(i) = -4*emissivity*(5.67*10^-8)*Td(1,i)^3;
                                        dH_Benn(i) = -1*density_air_0*P/P0*cA*A_Benn*u_PyrStat(i);
                                        if Rain_PyrStat(i) > 0
                                            dLE_Benn(i) = -0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*611*exp(-Lv/R*(1/Td(1,i)-1/273.15))*(Lv/R*Td(1,i)^-2);
                                        else
                                            dLE_Benn(i) = 0;
                                        end
                                        dP_Flux(i) = -density_water*cW*Rain_PyrStat(i)/(delta_t);
                                        dQc(i) = -k/h;
                                        dF_Ts_Benn(i) = dRn(i) + dH_Benn(i) + dLE_Benn(i) + dQc(i) + dP_Flux(i);
                                % Newton-Raphson method to solve for surface temperature
                                        while abs(Td(1,i)-Ts_past(i)) > 0.01 & n_iterations < 100
                                            n_iterations(i) = n_iterations(i) + 1;
                                            Ts_past(i) = Td(1,i);
                                            % max step size is 1 degree C
                                                Td(1,i) = Ts_past(i) - F_Ts_Benn(i)/dF_Ts_Benn(i);
                                                    if (Td(1,i) - Ts_past(i)) > 1
                                                        Td(1,i) = Ts_past(i) + 1;
                                                    elseif (Td(1,i) - Ts_past(i)) < -1
                                                        Td(1,i) = Ts_past(i) - 1;
                                                    end
                                     % Calculate temperature profile in the debris
                                            % For t = 0, which is i = 1, assume initial condition of linear temperature profile in the debris
                                                if i == 1
                                                    Td_gradient = (Td(1,1) - Td(N,1))/debris_thickness;
                                                    for j = 2:1:(N-1)
                                                        Td(j,1) = Td(1,1) - (j*h)*Td_gradient;
                                                    end
                                                else
                                                    % Perform Crank-Nicholson Scheme
                                                        for j = 2:1:(N-1)
                                                            % Equations A8 in Reid and Brock (2010) 
                                                                    a_Crank(j,i) = C;
                                                                    b_Crank(j,i) = 2*C+1;
                                                                    c_Crank(j,i) = C;

                                                            % Equations A9 in Reid and Brock (2010) 
                                                                    if j == 2
                                                                        d_Crank(j,i) = C*Td(1,i) + C*Td(1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                                    elseif j < (N-1)
                                                                        d_Crank(j,i) = C*Td(j-1,i-1) + (1-2*C)*Td(j,i-1) + C*Td(j+1,i-1);
                                                                    elseif j == (N-1)
                                                                        d_Crank(j,i) = 2*C*Td(N,i) + C*Td(N-2,i-1) + (1-2*C)*Td(N-1,i-1);
                                                                    end
                                                                        % note notation:
                                                                            % "i-1" refers to the past
                                                                            % "j-1" refers to the cell above it
                                                                            % "j+1" refers to the cell below it          

                                                            % Equations A10 and A11 in Reid and Brock (2010)
                                                                    if j == 2
                                                                        A_Crank(j,i) = b_Crank(j,i);
                                                                        S_Crank(j,i) = d_Crank(j,i);
                                                                    else
                                                                        A_Crank(j,i) = b_Crank(j,i) - a_Crank(j,i)/A_Crank(j-1,i)*c_Crank(j-1,i);
                                                                        S_Crank(j,i) = d_Crank(j,i) + a_Crank(j,i)/A_Crank(j-1,i)*S_Crank(j-1,i);
                                                                    end
                                                        end

                                                            % Equations A12 in Reid and Brock (2010)
                                                                    for j = N-1:-1:2
                                                                        if j == (N-1)
                                                                            Td(j,i) = S_Crank(j,i)/A_Crank(j,i);
                                                                        else
                                                                            Td(j,i) = 1/A_Crank(j,i)*(S_Crank(j,i)+c_Crank(j,i)*Td(j+1,i));
                                                                        end
                                                                    end
                                                end
                                            % Assume snow-free surface and compute fluxes normally
                                                    % Calculate Surface Energy Fluxes
                                                    if Rain_PyrStat(i) > 0
                                                        eS_Saturated(i) = 611*exp(-Lv/R*(1/Td(1,i)-1/273.15));
                                                        eS(i) = eS_Saturated(i);
                                                        eZ(i) = RH_PyrStat(i)*eZ_Saturated(i);
                                                        LE_Benn(i) = 0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*(eZ(i)-eS(i));
                                                    else
                                                        eS_Saturated(i) = 0;
                                                        LE_Benn(i) = 0;
                                                    end
                                                    Rn(i) = Sin(i)*(1-Albedo) + emissivity*(Lin_PyrStat(i) - (5.67*10^-8*Td(1,i)^4));
                                                    H_Benn(i) = density_air_0*(P/P0)*cA*A_Benn*u_PyrStat(i)*(Tair(i)-Td(1,i));
                                                    P_Flux(i) = density_water*cW*Rain_PyrStat(i)/(delta_t)*(Tair(i)-Td(1,i));
                                                    Qc(i) = k*(Td(2,i) - Td(1,i))/h;
                                                    F_Ts_Benn(i) = Rn(i) + H_Benn(i) + LE_Benn(i) + Qc(i) + P_Flux(i);
                                                    dRn(i) = -4*emissivity*(5.67*10^-8)*Td(1,i)^3;
                                                    dH_Benn(i) = -1*density_air_0*P/P0*cA*A_Benn*u_PyrStat(i);
                                                    if Rain_PyrStat(i) > 0
                                                        dLE_Benn(i) = -0.622*density_air_0/P0*Lv*A_Benn*u_PyrStat(i)*611*exp(-Lv/R*(1/Td(1,i)-1/273.15))*(Lv/R*Td(1,i)^-2);
                                                    else
                                                        dLE_Benn(i) = 0;
                                                    end
                                                    dP_Flux(i) = -density_water*cW*Rain_PyrStat(i)/(delta_t);
                                                    dQc(i) = -k/h;

                                                    dF_Ts_Benn(i) = dRn(i) + dH_Benn(i) + dLE_Benn(i) + dQc(i) + dP_Flux(i);
                                                if n_iterations == 100
                                                    Td(1,i) = (Td(1,i) + Ts_past(i)) / 2;
                                                end
                                        end
                                        Qc_ice(i) = k*(Td(N-1,i) - Td(N,i))/h;
                                        if Qc_ice(i) < 0
                                            Qc_ice(i) = 0;
                                        end
                                        Melt(i) = Qc_ice(i)*delta_t / (density_ice*Lf);
                                            % meters of ice melt
                        end % end Crank-Nicholson Newton-Raphson method for loop
                        
                    % Melt Table to compute monthly melt per year:
                        Table_Melt(:,1) = Year;
                        Table_Melt(:,2) = Month;
                        Table_Melt(:,3) = Day;
                        Table_Melt(:,4) = Julian_Day_of_year;
                        Table_Melt(:,5) = Melt';
                        
                    %% Compute the monthly melt per year based on script "Compute_Melt_MonthlyPerYear" 
                            % Step 1: Figure out the raw melt for each month
                            % Step 2: Compute average melt for that month
                            % Step 3: Adjust total melt for missing days
                            % Step 4: Compute the average melt for each month for all years
                            % Step 5: If entire month is missing, use avg melt from all the years
                            % Adjust Julian_Day for leap years to facilitate easier sorting/filling data gaps
                                Julian_Day_of_year_adjustedMelt = Julian_Day_of_year;
                                for i = 1:NumberTimeSteps
                                    if mod(Year(i),4) == 0
                                        Julian_Day_of_year_adjustedMelt(i) = Julian_Day_of_year(i) - 1;
                                    end
                                end
                            % Step 1: Compute the raw melt for each month 
                                Melt_MonthlyPerYear_Raw = zeros(3,6);
                                Melt_MonthlyPerYear_Raw_Count = zeros(3,6);
                                Melt_MonthlyPerYear_Total = zeros(3,6);
                                Melt_MonthlyPerYear_Avg = zeros(3,6);
                            for i = (days_spinup_nsteps+1):NumberTimeSteps
                                if Year(i) == 2010
                                    if Julian_Day_of_year_adjustedMelt(i) < 152 & Julian_Day_of_year_adjustedMelt(i) >= 135 % May
                                        Melt_MonthlyPerYear_Raw(1,1) = Melt_MonthlyPerYear_Raw(1,1) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(1,1) = Melt_MonthlyPerYear_Raw_Count(1,1) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 182 % June
                                        Melt_MonthlyPerYear_Raw(1,2) = Melt_MonthlyPerYear_Raw(1,2) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(1,2) = Melt_MonthlyPerYear_Raw_Count(1,2) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 213 % July
                                        Melt_MonthlyPerYear_Raw(1,3) = Melt_MonthlyPerYear_Raw(1,3) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(1,3) = Melt_MonthlyPerYear_Raw_Count(1,3) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 244 % Aug
                                        Melt_MonthlyPerYear_Raw(1,4) = Melt_MonthlyPerYear_Raw(1,4) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(1,4) = Melt_MonthlyPerYear_Raw_Count(1,4) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 274 % Sept
                                        Melt_MonthlyPerYear_Raw(1,5) = Melt_MonthlyPerYear_Raw(1,5) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(1,5) = Melt_MonthlyPerYear_Raw_Count(1,5) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 289 % Oct
                                        Melt_MonthlyPerYear_Raw(1,6) = Melt_MonthlyPerYear_Raw(1,6) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(1,6) = Melt_MonthlyPerYear_Raw_Count(1,6) + 1;
                                    end
                                elseif Year(i) == 2011
                                    if Julian_Day_of_year_adjustedMelt(i) < 152 & Julian_Day_of_year_adjustedMelt(i) >= 135 % May
                                        Melt_MonthlyPerYear_Raw(2,1) = Melt_MonthlyPerYear_Raw(2,1) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(2,1) = Melt_MonthlyPerYear_Raw_Count(2,1) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 182 % June
                                        Melt_MonthlyPerYear_Raw(2,2) = Melt_MonthlyPerYear_Raw(2,2) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(2,2) = Melt_MonthlyPerYear_Raw_Count(2,2) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 213 % July
                                        Melt_MonthlyPerYear_Raw(2,3) = Melt_MonthlyPerYear_Raw(2,3) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(2,3) = Melt_MonthlyPerYear_Raw_Count(2,3) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 244 % Aug
                                        Melt_MonthlyPerYear_Raw(2,4) = Melt_MonthlyPerYear_Raw(2,4) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(2,4) = Melt_MonthlyPerYear_Raw_Count(2,4) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 274 % Sept
                                        Melt_MonthlyPerYear_Raw(2,5) = Melt_MonthlyPerYear_Raw(2,5) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(2,5) = Melt_MonthlyPerYear_Raw_Count(2,5) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 289 % Oct
                                        Melt_MonthlyPerYear_Raw(2,6) = Melt_MonthlyPerYear_Raw(2,6) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(2,6) = Melt_MonthlyPerYear_Raw_Count(2,6) + 1;
                                    end
                                elseif Year(i) == 2012
                                    if Julian_Day_of_year_adjustedMelt(i) < 152 & Julian_Day_of_year_adjustedMelt(i) >= 135 % May
                                        Melt_MonthlyPerYear_Raw(3,1) = Melt_MonthlyPerYear_Raw(3,1) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(3,1) = Melt_MonthlyPerYear_Raw_Count(3,1) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 182 % June
                                        Melt_MonthlyPerYear_Raw(3,2) = Melt_MonthlyPerYear_Raw(3,2) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(3,2) = Melt_MonthlyPerYear_Raw_Count(3,2) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 213 % July
                                        Melt_MonthlyPerYear_Raw(3,3) = Melt_MonthlyPerYear_Raw(3,3) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(3,3) = Melt_MonthlyPerYear_Raw_Count(3,3) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 244 % Aug
                                        Melt_MonthlyPerYear_Raw(3,4) = Melt_MonthlyPerYear_Raw(3,4) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(3,4) = Melt_MonthlyPerYear_Raw_Count(3,4) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 274 % Sept
                                        Melt_MonthlyPerYear_Raw(3,5) = Melt_MonthlyPerYear_Raw(3,5) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(3,5) = Melt_MonthlyPerYear_Raw_Count(3,5) + 1;
                                    elseif Julian_Day_of_year_adjustedMelt(i) < 289 % Oct
                                        Melt_MonthlyPerYear_Raw(3,6) = Melt_MonthlyPerYear_Raw(3,6) + Melt(i);
                                        Melt_MonthlyPerYear_Raw_Count(3,6) = Melt_MonthlyPerYear_Raw_Count(3,6) + 1;
                                    end
                                end
                            end
                            % Step 2: Compute average for each month per year
                                % Days per melt that were recorded
                                    Melt_MonthlyPerYear_Raw_Count_days = Melt_MonthlyPerYear_Raw_Count/24;
                                % If the entire month isn't missing, then compute the average
                                    for r_m = 1:3
                                        for c_m = 1:6
                                            if Melt_MonthlyPerYear_Raw_Count_days(r_m,c_m) > 0
                                                Melt_MonthlyPerYear_Avg(r_m,c_m) = Melt_MonthlyPerYear_Raw(r_m,c_m)/Melt_MonthlyPerYear_Raw_Count_days(r_m,c_m);
                                            end
                                        end
                                    end
                            % Step 3: Adjust total melt for missing days each month each year
                                    for r_m = 1:3
                                        for c_m = 1:6
                                            if Melt_MonthlyPerYear_Raw_Count_days(r_m,c_m) > 0
                                                Melt_MonthlyPerYear_Total(r_m,c_m) = Melt_MonthlyPerYear_Avg(r_m,c_m)*Table_DaysPerMonthPerYear(r_m,c_m);
                                            end
                                        end
                                    end
                            % Step 4: Compute average melt per month for all years that have some
                                Melt_MonthlyAllYears_Avg = sum(Melt_MonthlyPerYear_Total,1)./sum(Melt_MonthlyPerYear_Total~=0,1);
                                    % Note: sum(A~=0,1) counts the number of non-zero values in each column of matrix A
                                % Replace NaN values with zero in case there was no melt for a given month 
                                    Melt_MonthlyAllYears_Avg(isnan(Melt_MonthlyAllYears_Avg)==1)=0;
                            % Step 5: If entire month is missing, use avg monthly melt rate of all years
                                for r_m = 1:3
                                    for c_m = 1:6
                                        if Melt_MonthlyPerYear_Total(r_m,c_m) == 0 & Table_DaysPerMonthPerYear(r_m,c_m) > 0
                                            Melt_MonthlyPerYear_Total(r_m,c_m) = Melt_MonthlyAllYears_Avg(1,c_m);
                                        end
                                    end
                                end
                                
                        Melt_Modeled_pixel = sum(sum(Melt_MonthlyPerYear_Total));
                            % for May 15 - Oct 15 Melt season
                        
            end % end while loop that varies debris thickness
            
            % If band is not debris cover (Elevation change = 0), then
            % skips the while loop and put 0 for debris thickness and melt
            if bands_DEM_Diff_raw(n,1) == 0
                debris_thickness_rounded = 0;
                Melt_Modeled_pixel = 0;
            else
                debris_thickness_rounded = round(debris_thickness/debris_thickness_round_interval)*debris_thickness_round_interval;
            end
            
            % Write Results of MC Simulations
                Table_MC_Simulation_Output = [n, MC + MC_extra, debris_thickness_rounded, Melt_Modeled_pixel, Albedo, z0, k, DEM_Diff_random(MC,1), DEM_Diff_noem, emergence_random(MC,1), emergence_velocity, DEM_Diff, N_iterations-1];
                fprintf(fileID,'%1d %1d %3.2f %5.3f %3.2f %4.3f %3.2f %6.3f %6.3f %6.3f %6.3f %6.3f %4d\r\n',Table_MC_Simulation_Output);
                Display_Table = [n, MC];
                display(Display_Table)
                
        if MC == MC_Simulations_Total;
            fclose(fileID);
        end

    end % end parfor loop
end % end n = 1:n_bands

fclose('all');
% Turn off timer
    toc

% Notes: n_iterations = 2 indicates debris thickness is outside of bounds
%        Debris thickness set to 0 or maximum (not higher than 5 m according to Nicholson's dissertation) as a result