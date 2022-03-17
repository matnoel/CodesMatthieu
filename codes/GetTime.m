function temps = GetTime(time)
% convert time in segondes to xh: xm: xs
        hour = floor(time/3600);
        hourRest = time/3600 - hour;
        minutes = floor(hourRest*60);
        minutesRest = hourRest*60-minutes;
        segondes = floor(minutesRest*60);

        temps = ""+hour+"h: "+minutes+"m: "+segondes+"s";
%         temps = ('%2dh:%2.fm:%2.fs',hour, minutes, segondes);
end