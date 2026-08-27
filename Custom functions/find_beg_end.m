function start_end = find_beg_end(boolean_val)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

start_end = [find(diff([0 boolean_val 0])==1)' find(diff([0 boolean_val 0])==-1)'-1];
end