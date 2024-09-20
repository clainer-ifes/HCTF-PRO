function [Y_cumulative]=cumulative_distribution(Y)

    Y_cumulative = zeros(1, length(Y));
    for index=1:length(Y)
        if(index==1)
            TermToBeAdded = 0;
        else
            TermToBeAdded = Y_cumulative(index-1);
        end
        
        Y_cumulative(index) = TermToBeAdded + Y(index);
    end
    
end
