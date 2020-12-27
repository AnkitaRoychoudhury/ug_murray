def get_colors(data_type,labels):
    '''data_type = string input (section the colors will be chosen from)
        options are [atp_controls]
    label = list of string values 
    output: colormap'''
    
    # Get data for atp controls
    if data_type == 'atp_controls':
        color_dict = {'Positive Control':'#1f77b3',
               'Negative Control':'#ff7e0e',
               'ATP Depletion':'#2ba02b',
               'No Buffer':'#d62628',
               'No Extract':'#9367bc',
               'Positive Control 2':'#1f77b3'}
    # Output proper cmap
    cmap = []
    for label in labels:
        cmap.append(color_dict[label])
        
    return cmap
    