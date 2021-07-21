function display_progress(i, MaxIter)
    display(strcat("Percentage complete:", num2str(floor(i/MaxIter*100)), '%'))
