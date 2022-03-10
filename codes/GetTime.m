function [hour, minutes, segondes] = GetTime(time)
        hour = floor(time/3600);
        hourRest = time/3600 - hour;
        minutes = floor(hourRest*60);
        minutesRest = hourRest*60-minutes;
        segondes = floor(minutesRest*60);
end