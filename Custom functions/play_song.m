function [] = play_song(notes,duration, speed)

if isempty(notes)
notes       = [0 4 4 4 2 0 7 7 5 4 4 4 2 0 7];
end
if isempty(duration)
duration    = [8 4 4 4 8 8 2 8 8 4 4 4 8 8 2];
end
if isempty(speed)
    speed =1;
end

for j=1:numel(notes)

sound(sin(pi*(880*(2^(notes(j)/12)))*(1:(speed*200000/duration(j)))/200000)', 200000)
pause(speed/duration(j))
end

end