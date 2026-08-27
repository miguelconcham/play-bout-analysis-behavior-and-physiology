function [run_events_table] = find_runs(all_rate,treshold,rate_wind,this_neuron_spikes, min_spikes, run_length_limits,max_coverage)
all_rate_time = (1:numel(all_rate))*rate_wind;
starting_event = 0;
run_events_table = [];
all_positions = 1:numel(all_rate);
while starting_event<numel(all_rate)

    current_event = min(find(all_rate'>=treshold & all_positions>starting_event));


    if ~isempty(current_event)

        current_run_range = [min(find(all_rate'>0  & all_positions<= current_event & all_positions>starting_event)) min(find(all_rate'==0  & all_positions>=current_event ))-1];
        if numel(current_run_range)==2
        [max_current_rate, loc] = max(all_rate(current_run_range(1):current_run_range(2)));
        max_loc = current_run_range(1)+loc-1;

        current_loc = max_loc;
        lower_bound_count = 0;
        searching = true;
        while searching && current_loc<numel(all_rate)
            if all_rate(current_loc)>=0.01*max_current_rate
                lower_bound_count=0;
            else
                lower_bound_count = lower_bound_count+1;
                if lower_bound_count>=.250/rate_wind
                    searching= false;
                end
            end

            current_loc = current_loc+1;
        end
        upper_end = current_loc-.250/rate_wind -1;
        % all_rate(upper_end)

        current_loc = max_loc;
        lower_bound_count = 0;
        searching = true;
        while searching && current_loc>0

            if all_rate(current_loc)>=0.01*max_current_rate
                lower_bound_count=0;
            else
                lower_bound_count = lower_bound_count+1;
            end

            if lower_bound_count>=.250/rate_wind
                searching= false;
            end
            current_loc = current_loc-1;
        end
        lower_end = current_loc+.250/rate_wind +1;
        % all_rate(lower_end)
        event_time_range = all_rate_time([lower_end upper_end]);
        if range(event_time_range)>run_length_limits(2)*2
            n_splits = floor(range(event_time_range));
            new_ranges =  linspace(event_time_range(1), event_time_range(2), n_splits+1);
            for sd = 1:numel(new_ranges)-1
                num_spikes = sum(this_neuron_spikes>=new_ranges(sd) & this_neuron_spikes<=new_ranges(sd+1));
                this_sub_run_indexes = find(all_rate_time>=new_ranges(sd) & all_rate_time<=new_ranges(sd+1));
                [max_current_rate, loc]=  max(all_rate(this_sub_run_indexes));
                loc = this_sub_run_indexes(loc);
                run_events_table = [run_events_table;[this_sub_run_indexes(1) this_sub_run_indexes(end) all_rate_time(max_loc) max_current_rate (new_ranges(sd+1)-new_ranges(sd)) num_spikes]];

            end
        else

            num_spikes = sum(this_neuron_spikes>=event_time_range(1) & this_neuron_spikes<=event_time_range(2));
            run_events_table = [run_events_table;[lower_end upper_end all_rate_time(max_loc) max_current_rate range(event_time_range) num_spikes]];
        end
        if upper_end<current_run_range(2)
            starting_event = current_run_range(2);
        else
            starting_event = upper_end;
        end
        else
            starting_event = Inf;
        end
    else
        starting_event = Inf;
    end


end
 
if isempty(run_events_table)
    run_events_table = nan(1,8);
    run_events_table = array2table(run_events_table);
    run_events_table.Properties.VariableNames = {'RunStart', 'RunEnd', 'MaxRateLoc', 'RunMaxRate', 'RunLength', 'RunNumSpikes', 'RunStartTime', 'RunEndTime'};

else
run_events_table = array2table(run_events_table);
run_events_table.Properties.VariableNames = {'RunStart', 'RunEnd', 'MaxRateLoc', 'RunMaxRate', 'RunLength', 'RunNumSpikes'};
run_events_table = run_events_table(run_events_table.RunLength>=run_length_limits(1) & run_events_table.RunNumSpikes>=min_spikes,:);
run_events_table.RunStartTime = all_rate_time(run_events_table.RunStart)';
run_events_table.RunEndTime = all_rate_time(run_events_table.RunEnd)';
end

if isempty(run_events_table)
    run_events_table = nan(1,8);
    run_events_table = array2table(run_events_table);
    run_events_table.Properties.VariableNames = {'RunStart', 'RunEnd', 'MaxRateLoc', 'RunMaxRate', 'RunLength', 'RunNumSpikes', 'RunStartTime', 'RunEndTime'};
end

if sum(run_events_table.RunLength)/range(all_rate_time)>max_coverage
    run_ranges = round(linspace(1,numel(all_rate), round(numel(all_rate)/(20*run_length_limits(1)/rate_wind)) +1));

    run_indexes =[run_ranges(1:end-1)', run_ranges(2:end)'];
    run_events_table = array2table(run_indexes);
    run_events_table.Properties.VariableNames ={'RunStart', 'RunEnd'};

    run_events_table.MaxRateLoc = nan(size(run_events_table,1),1);
    run_events_table.RunMaxRate = nan(size(run_events_table,1),1);
    run_events_table.RunLength = nan(size(run_events_table,1),1);
    run_events_table.RunNumSpikes = nan(size(run_events_table,1),1);
     run_events_table.RunStartTime = nan(size(run_events_table,1),1);
      run_events_table.RunEndTime = nan(size(run_events_table,1),1);

      for j=1:size(run_events_table,1)

          [run_events_table.RunMaxRate(j), run_events_table.MaxRateLoc(j)] = max(all_rate(run_events_table.RunStart(j):run_events_table.RunEnd(j)));
          run_events_table.MaxRateLoc(j) = run_events_table.MaxRateLoc(j)+ run_events_table.RunStart(j);
          run_events_table.RunEndTime(j) = all_rate_time(run_events_table.RunEnd(j));
          run_events_table.RunStartTime(j) = all_rate_time(run_events_table.RunStart(j));

          run_events_table.RunLength(j) = run_events_table.RunEndTime(j)-  run_events_table.RunStartTime(j);

           run_events_table.RunNumSpikes(j) = sum(this_neuron_spikes>=run_events_table.RunStartTime(j) & ...
               this_neuron_spikes<=run_events_table.RunEndTime(j));
      end
      run_events_table = run_events_table(run_events_table.RunLength>=run_length_limits(1) & run_events_table.RunNumSpikes>=min_spikes,:);

end







end