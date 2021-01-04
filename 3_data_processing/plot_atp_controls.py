import holoviews as hv
import pandas as pd
import bokeh.io
import bokeh.plotting
from get_colors import get_colors

def atp_controls(df, title = '', custom_df = False, x_axis_label = 'Time (hr)', y_axis_label = 'Raw Fluorescence Data', well_dict = {},
                legend_location = 'bottom_right', glyph = 'Curve'):
    '''df = dataframe with Time Point, Fluor Value, Descriptive Well Cols columns.
    title = title of outputted plot
    custom_df = False if df from murraylab_tools, True otherwise
    x_axis_label = for custom_df = False bokeh plot (default = 'Time (hr)')
    y_axis_label = for custom_df = False bokeh plot (default = 'Raw Fluorescence Data')
    well_dict = needed for custom_df = False (default = {}),
    replicates = 'No' or 'Yes' ('Yes' example in 20201107_exp10_analysis)
    Outputs holoviews curve object'''
    if custom_df == False:
        try:
            def y_vals(df, well):
                df2 = df.loc[df['Well'] == well, :]
                return df2['Measurement'].values
            curve = bokeh.plotting.figure(height = 350, width = 650,
                                         title = title,
                                         x_axis_label = x_axis_label, y_axis_label = y_axis_label)
            
            time_vals = df['Time (hr)'].unique()
            all_wells = list(well_dict.keys())
            # Get colors
            well_values = list(well_dict.values())
            cmap = get_colors('controls_descriptive', well_values)
            # Plot lines
            for i, curr_well in enumerate(all_wells):
                curve.line(time_vals, y_vals(df, curr_well), legend_label = well_dict[curr_well], color = cmap[i])
            curve.legend.click_policy = 'hide'
            curve.legend.location = legend_location
        except:
            print("Double check inputs")
            

        
    elif custom_df == True:  
        try:
            curve = hv.Curve(data = df,
                 kdims = ['Time Point', 'Fluor Value'],
                 vdims = ['Descriptive Well Cols']
                            ).groupby('Descriptive Well Cols'
                            ).overlay("Descriptive Well Cols"
                            ).opts(title = title)
      
        except:
            print("Incorrect input of dataframe. Make sure there are Time Point, Fluor Value, Descriptive Well Cols columns.")
            
    if glyph == "Points":
        try:
            curve = hv.Points(data = df,
                 kdims = ['Time Point', 'Fluor Value'],
                 vdims = ['Descriptive Well Cols']
                            ).groupby('Descriptive Well Cols'
                            ).overlay("Descriptive Well Cols"
                            ).opts(title = title,
                                  xlabel = 'Time (hrs)',
                                  ylabel = 'ATP Concentration [mM]')
        except:
            print("Something wrong with replicate part")
    return curve

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
    