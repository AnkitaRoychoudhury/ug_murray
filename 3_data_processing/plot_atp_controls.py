import holoviews as hv
def atp_controls(df, title = ''):
    '''df = dataframe with Time Point, Fluor Value, Descriptive Well Cols columns.
    title = title of outputted plot
    Outputs holoviews curve object'''

        
    try:
        curve = hv.Curve(data = df,
             kdims = ['Time Point', 'Fluor Value'],
             vdims = ['Descriptive Well Cols']
                        ).groupby('Descriptive Well Cols'
                        ).overlay("Descriptive Well Cols"
                        ).opts(title = title)
        return curve
        
    except:
        print("Incorrect input of dataframe. Make sure there are Time Point, Fluor Value, Descriptive Well Cols columns.")
    
def atp_controls_w_color(df, title = ''):
    ''' MAKES CUSTOM CMAP WITH get_colors FUNCTION. 
    DOESNT WORK BECAUSE HOLOVIEWS
    df = dataframe with Time Point, Fluor Value, Descriptive Well Cols columns.
    title = title of outputted plot,
    Outputs holoviews curve object'''
    try:
        cmap = get_colors('atp_controls',df['Descriptive Well Cols'])
        print(cmap)
    except:
        print('Could not retrieve color map')
        
    try:
        curve = hv.Curve(data = df,
             kdims = ['Time Point', 'Fluor Value'],
             vdims = ['Descriptive Well Cols']
                        ).groupby('Descriptive Well Cols'
                        ).overlay("Descriptive Well Cols"
                        ).opts(title = title,
                                  cmap = cmap)
        return curve
        
    except:
        print("Incorrect input of dataframe. Make sure there are Time Point, Fluor Value, Well Col columns.")
    