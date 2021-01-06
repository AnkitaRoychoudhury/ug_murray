import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()
import holoviews as hv
hv.extension("bokeh")
import pandas as pd
from get_colors import get_colors

def plot_atp_assay(df, atp = 'ATP Concs', title = 'ATP Luminescence End Points'):
    '''Inputs:
    df = dataframe with ATP Values and Descriptive Well Cols as columns
    Outputs:
    holoviews scatter plot'''
    
    wells = df['Descriptive Well Cols'].unique()
    cmap = get_colors('beg/spike', wells)
    
    a = hv.Scatter(data = df,
         kdims = [atp],
          vdims = ['Descriptive Well Cols']).opts(
        color = 'Descriptive Well Cols',
        height = 450,
        width = 500,
        xlabel = atp,
        cmap = cmap,
        size=8,
        alpha = 1,
        show_legend = False,
        title = title,
)
    
    return a
    