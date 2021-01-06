def get_colors(data_type,labels):
    '''data_type = string input (section the colors will be chosen from)
        options are [atp_controls, controls_descriptive, atp added, beg/spike]
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
        
    elif data_type == 'atp added':
        color_dict = {'Positive Control':'#1f77b3',
                      'Positive Control (PC)':'#1f77b3',
                    'PC + 10 uL 10 mM ATP': '#980043',
                    'PC + 1 uL 10 mM ATP': '#df65b0',
                    'PC + 1 uL 1 mM ATP': '#d4b9da',
                    'Negative Control':'#d4b9da'}
        
    elif data_type == 'beg/spike':
        color_dict = {'Positive Control':'#e41a1c',
                      'Positive Control (PC)':'#1f77b3',
                     'Negative Control': '#377eb8',
                     'PC + 50 mM ATP':'#4daf4a' ,
                     'PC + 25 mM ATP': '#984ea3',
                      'PC + 1.4 M 3PGA': '#cab2d6',
                      'PC + Energy Buffer': '#a65628',
                      'PC + DNA': '#f781bf',
                      'PC + Extract': '#999999',
                      'PC + 175 mM NAD':'#ff7f00'}
        
#     elif data_type == 'atp assay':
#         color_dict = {'Positive Control':'#e41a1c' ,
#                      'Negative Control': '#377eb8',
#                      'PC + 50 mM ATP':'#4daf4a' ,
#                      'PC + 25 mM ATP': '#984ea3',
#                       'PC + 1.4 M 3PGA': '#ffff33',
#                       'PC + Energy Buffer': '#a65628',
#                       'PC + DNA': '#f781bf',
#                       'PC + Extract': '#999999',
#                       'PC + 175 mM NAD':'#ff7f00'}
                      
            
#         }
        
    # Output proper cmap
    cmap = []
    for label in labels:
        #print(label)
        cmap.append(color_dict[label])
        
    return cmap
    