function [Y] = loadSensorNames(filePath)
% this function was written to load the sensor names for the MetaDatUse.txt
% file
fid = fopen(filePath);
tline = fgetl(fid);
%disp('Triying to read');
Z{1,:} = tline(1:5);

senCount = 2;
while ~feof(fid)
    tline = fgetl(fid);
    %disp(tline);
    Z{senCount,:} = tline(1:5);
    tline(1:5);
    senCount = senCount+1;
    %disp(senCount);
end

% conver it to a string array
Y = convertCharsToStrings(Z);
disp('All lines read')
end

