function [allTime, allTemp] = loadScrippsPier()
% This function loads multiple years of Scripps Data
% OUTPUTS:
%   allTime - Datetime array containing the date and time of each measurement
%   allTemp - Temperature measurements from all years

% Load all Scripps Pier temperature data (in a loop)
years= 2005:2021;

% Blanks
allTime= datetime.empty(0,1);
allTemp= [];

% Define the reference date for the NetCDF time values
d0= datetime(1970,1,1);

% Loop through each year
for i= 1:length(years)

    % Filename
    filename= sprintf("scripps_pier-%d.nc", years(i));

    % Read time and temperature from the NetCDF file
    time= ncread(filename, 'time');
    temp= ncread(filename, 'temperature');

    % Convert time from seconds since January 1, 1970 to datetime
    dateTime= d0 + seconds(double(time));

    % Append this year's data to the combined arrays
    allTime= [allTime; dateTime(:)];
    allTemp= [allTemp; double(temp(:))];
end

% Make sure data are chronological
[allTime, sortIndex] = sort(allTime);
allTemp = allTemp(sortIndex);

end