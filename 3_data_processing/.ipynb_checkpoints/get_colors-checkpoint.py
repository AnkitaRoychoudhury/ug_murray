def get_colors(data_type,labels):
    '''data_type = string input (section the colors will be chosen from)
        options are [atp_controls, controls_descriptive]
    label = list of string values 
    output: colormap'''
    color_dict = {}
    # Get data for atp controls
    if data_type == 'atp_controls':
        color_dict = {'Positive Control':'#1f77b3',
               'Negative Control':'#ff7e0e',
               'ATP Depletion':'#2ba02b',
               'No Buffer':'#d62628',
               'No Extract':'#9367bc',
               'Positive Control 2':'#1f77b3'}
    elif data_type == 'controls_descriptive':
        color_dict = {'Positive Control (Buffer + Extract + DNA)':'#1f77b3',
               'Neg Control (Extract Only)':'#ff7e0e',
               'ATP Depletion (Buffer + Extract)':'#2ba02b',
               'No Energy (Extract + DNA)':'#d62628',
               'No TXTL (Buffer + DNA)':'#9367bc'}
        
    # Output proper cmap
    cmap = []
    for label in labels:
        print(label)
        cmap.append(color_dict[label])
        
    return cmap
    