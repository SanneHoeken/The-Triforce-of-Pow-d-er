
def calculate_coordinate(fold, coordinate):

    co = list(coordinate)
    
    # adjust x-value following the fold
    if fold == 1 or fold == -1:
        co[0] += fold

    # adjust y-value following the fold    
    elif fold == 2:
        co[1] += 1

    elif fold == -2:
        co[1] -= 1

    # adjust z-value following the fold
    elif fold == 3:
        co[2] += 1
    
    elif fold == -3:
        co[2] -= 1
            
    return tuple(co)