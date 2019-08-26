function [time_strct]=sec2struct(x)
% SEC2STRUCT splits seconds into hours, minutes, seconds and milliseconds 
% INPUT:
%	x -> number of seconds (scalar value, floating point)
% OUTPUT:
%	time_strct -> structure with the following fields
%		.hour
%		.minutes
%		.seconds
%		.milliseconds
%		.str: time as HH:MM:SS:MS string
%		.vec: time as HHMMSSMS string
%
% Author: Ulrich Schmidt, April 2011
%

time_strct.hour=int8(floor(x/3600.0));
time_strct.minutes=int8(floor(mod((x/60.0), 60.0)));
time_strct.seconds=int8(floor(mod(x,60.0)));
time_strct.milliseconds=int16(mod(x,1.0)*1000);

time_strct.str=sprintf('%02d:%02d:%02d:%03d',...
    time_strct.hour,...
    time_strct.minutes,...
    time_strct.seconds,...
    time_strct.milliseconds...
);

time_strct.vec=sprintf('%02d%02d%02d.%03d',...
    time_strct.hour,...
    time_strct.minutes,...
    time_strct.seconds,...
    time_strct.milliseconds...
);