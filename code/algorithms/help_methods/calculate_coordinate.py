
def calculate_coordinate(fold, x, y):
     
        # compute x-value following the fold
        if fold == 1 or fold == -1:
            x_tmp = x + fold

        # keep current x-value if no move along x-axis 
        else:
            x_tmp = x

        # compute y-value following the fold    
        if fold == 2 or fold == -2:
            y_tmp = y + 0.5 * fold
        
        # keep current y-value if no move along y-axis
        else:
            y_tmp = y
        
        return int(x_tmp), int(y_tmp)