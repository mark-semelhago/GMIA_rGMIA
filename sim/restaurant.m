function neg_profit = restaurant(design_point, c, T, burnin, lambda, mu, rev, l, rep)
%Restaurant simulator simulates the negative profit (which we want to
%minimize) from running a restaurant, where you have a certain number of
%tables of different capacities. Groups of customers arrive according to a
%Poisson process with arrival rates lambda and service rates mu. Each group
%brings a certain amount of revenue and each table of different size has a
%per hour, per table cost to keep available and staffed.
%design_point is the number of tables for each capacity (dp_num x rep).
%c is the table capacities.
%T is the length of time the restaurant is open for.
%burnin is the time at which you begin counting the arrivals.
%lambda is the rates of groups (group/unit time).
%mu is the service rate associated with each group size.
%rev is the revenue associated with each group size (per group).
%l is the per hour, per table cost involved with having a table available.
%rep is the number of replications to produce.
    neg_profit = zeros(size(design_point,1),rep);
    for dp = 1:size(design_point,1)
        for r = 1:rep
            c_max = max(c);

            arrival_times = cell(c_max,1);
            service_times = cell(c_max,1);

            for i = 1:c_max
                t = exprnd(1/lambda(i));
                while (t < T)
                    arrival_times{i} = [arrival_times{i};t];
                    service_times{i} = [service_times{i};exprnd(1/mu(i))];
                    dt = exprnd(1/lambda(i));
                    t = t + dt;
                end
            end
            %timeline stores the events that occur
            %Timing stores the time at which the action occurs
            %Action is either 'A' for arrival or 'D' for departure
            %Group is the size of the group that has arrived if Action = 'A' or the
            %size of the table that has been freed if Action = 'D'
            arrivals_timeline = [];
            departure_timeline = [];
            for i = 1:c_max
                arrivals_timeline = [arrivals_timeline;arrival_times{i}, i*ones(length(arrival_times{i}),1), arrival_times{i} + service_times{i}];
            end


            tables_avail = [c, design_point(dp,:)'];
            profit = 0;

            %Only include arrivals that arrived after burnin time
            arrivals_timeline = arrivals_timeline(arrivals_timeline(:,1)>=burnin,:);

            while ~isempty(arrivals_timeline)
                [m_a, i_a] = min(arrivals_timeline(:,1));
                if ~isempty(departure_timeline)
                    [m_d, i_d] = min(departure_timeline(:,1));
                else
                    m_d = [];
                    i_d = [];
                end
                if m_d < m_a
                    tables_avail(find(tables_avail(:,1)==departure_timeline(i_d,2)),2) = tables_avail(find(tables_avail(:,1)==departure_timeline(i_d,2)),2)+1;
                    departure_timeline = [departure_timeline(1:i_d-1,:);departure_timeline(i_d+1:end,:)];
                else
                    remove_table = tables_avail(tables_avail(:,1)>=arrivals_timeline(i_a,2) & tables_avail(:,2) > 0,:);
                    if ~isempty(remove_table)
                        remove_table = remove_table(1,1);
                        tables_avail(find(tables_avail(:,1)==remove_table),2) = tables_avail(find(tables_avail(:,1)==remove_table),2) - 1;
                        departure_timeline = [departure_timeline;arrivals_timeline(i_a,3), remove_table];
                        profit = profit + rev(arrivals_timeline(i_a,2));
                    end
                    arrivals_timeline = [arrivals_timeline(1:i_a-1,:);arrivals_timeline(i_a+1:end,:)];
                end
            end
        profit = profit - T*l'*design_point(dp,:)';        
        neg_profit(dp,r) = -profit;
        end        





    end

end
