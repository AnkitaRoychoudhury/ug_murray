import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()
import pandas as pd
import murraylab_tools.biotek as mt_biotek
import os
from get_colors import get_colors
import holoviews as hv
hv.extension("bokeh")

def y_vals(df, well):
    '''Used in plot_timeseries to extract desired data,
    helper function'''
    df2 = df.loc[df['Well'] == well, :]
    return df2['Measurement'].values

def get_tidy_df(data_folder, data_name, gain = 61):
    '''Use murray labtools to get tidy dataframe, 
    helper function (can also be used)'''
    inputfolder = data_folder
    outputfolder = data_folder
    data_name = data_name
    data_filename = os.path.join(inputfolder, data_name+".csv")
    tidy_filename = os.path.join(outputfolder,data_name+"_tidy.csv")
    mt_biotek.tidy_biotek_data(data_filename, volume = 10,
                               output_filename = tidy_filename)
    df = pd.read_csv(tidy_filename)
    
    df_gain = df.loc[df['Gain'] == gain, :]
    
    return df_gain

def plot_timeseries(data_folder, data_name,
                   well_dict, data_type = 'atp added', gain = 61, title = '', 
                   x_axis_label = 'Time (hrs)',
                   y_axis_label = 'Raw Fluorescence Data'):
    '''Inputs:
    data_folder = directory where data can be found and will be outputted
    data_name = name of the file
    well_dict = dictionary with well: descriptive well column
    data_type = str to be entered to get_colors
    gain = the gain to plot (default = 61)
    title = title of plot
    x_axis_label, y_axis_label = labels of plot
    Output:
    bokeh.plotting figure'''
    
    # Get in tidy format with murraylabtools
    df_gain = get_tidy_df(data_folder, data_name)
    
    # Get colors
    color_vals = get_colors(data_type = data_type, labels = list(well_dict.values()))
    
    # Plot
    p = bokeh.plotting.figure(height = 350, width = 650,
                             title = title,
                             x_axis_label = x_axis_label,
                             y_axis_label = y_axis_label)
    time_vals = df_gain['Time (hr)'].unique()
    
    all_wells = list(well_dict.keys())
    
    for i,curr_well in enumerate(all_wells):
        p.line(time_vals, y_vals(df_gain, curr_well), 
               legend_label = well_dict[curr_well],
              color = color_vals[i])
    p.legend.click_policy = 'hide'
    p.legend.location = 'bottom_right'
    
    return p


def combine_pre_post(data_folder, data_name_pre, data_name_post,
                    well_dict, data_type = 'atp added', gain = 61, title = '', 
                   x_axis_label = 'Time (hrs)',
                   y_axis_label = 'Raw Fluorescence Data'):
     
    '''Inputs:
    data_folder = directory where data can be found and will be outputted
    data_name_pre = name of the file for prespike data
    data_name_post = name of the file for postspike data
    well_dict = dictionary with well: descriptive well column
    data_type = str to be entered to get_colors
    gain = the gain to plot (default = 61)
    title = title of plot
    x_axis_label, y_axis_label = labels of plot
    Output:
    bokeh.plotting figure'''
    # PRESPIKE DATA
    # Get in tidy format with murraylabtools
    df_gain_pre = get_tidy_df(data_folder, data_name_pre)

    # Add time 2 column (to align pre and post data)
    df_gain_pre['Time2'] = df_gain_pre['Time (hr)']
    
    # POSTPSIKE DATA
    # Get in tidy format with murraylabtools
    df_gain_post = get_tidy_df(data_folder, data_name_post)
    
    # Get last time and make new column
    last_prespike_time = list(df_gain_pre['Time (hr)'])[-1]
    df_gain_post['Time2']= df_gain_post['Time (hr)'] + last_prespike_time
    
    # Concatenate pre and post lists
    df_both = pd.concat([df_gain_pre, df_gain_post])
    
    # Get colors and time
    color_vals = get_colors(data_type = data_type, labels = list(well_dict.values()))
    time_vals = df_both['Time2'].unique()

    # Plot
    p = bokeh.plotting.figure(height = 350, width = 650,
                             title = title,
                             x_axis_label = x_axis_label,
                             y_axis_label = y_axis_label)
    
    all_wells = list(well_dict.keys())
    
    for i,curr_well in enumerate(all_wells):
        p.line(time_vals, y_vals(df_both, curr_well), 
               legend_label = well_dict[curr_well],
              color = color_vals[i])
    p.legend.click_policy = 'hide'
    p.legend.location = 'bottom_right'
    
    return p

def plot_pre_post_endpt(df_pre, df_post, well_dict, data_type = 'beg/spike',
                       title = ''):
    '''Inputs:
    df_pre = dataframe with pre spike time
    df_post = dataframe with post spike time
    necessary columns - [Time (hr), Measurement]
    well_dict = dictionary with well: descriptive well column 
    data_type = for color identifier (default = 'beg/spike')
    Output:
    hv.Scatter plot'''

    # Get dataframes for last prespike and last postspike times
    last_pre_time = list(df_pre['Time (hr)'])[-1]
    last_post_time = list(df_post['Time (hr)'])[-1]

    df_last_pre = df_pre.loc[df_pre['Time (hr)'] == last_pre_time, :]
    df_last_post = df_post.loc[df_post['Time (hr)'] == last_post_time,:]

    # Add descriptive well col labels
    desc_list_pre = []
    desc_list_post = []
    for well in df_last_pre['Well']:
        desc_list_pre.append(well_dict[well])

    for well in df_last_post['Well']:
        desc_list_post.append(well_dict[well])

    df_last_pre['Descriptive Well Cols'] = desc_list_pre
    df_last_post['Descriptive Well Cols'] = desc_list_post

    # Get all average measurements and make summary df
    df_summary = pd.DataFrame()
    mean_list = []
    well_list = []
    for well in df_last_pre['Descriptive Well Cols']:
        curr_df = df_last_pre.loc[df_last_pre['Descriptive Well Cols'] == well, :]
        curr_mean = curr_df.mean()['Measurement']
        mean_list.append(curr_mean)
        well_list.append(well)

    for well in df_last_post['Descriptive Well Cols']:
        curr_df = df_last_post.loc[df_last_post['Descriptive Well Cols'] == well, :]
        curr_mean = curr_df.mean()['Measurement']
        mean_list.append(curr_mean)
        well_list.append(well)

    df_summary['Well'] = well_list
    df_summary['Avg Measurement'] = mean_list
    

    # Make plot
    color_vals = get_colors(data_type = data_type, labels = well_list)
    
    a = hv.Scatter(data = df_summary,
         kdims = ['Avg Measurement'],
          vdims = ['Well']).opts(
        color = 'Well',
        height = 450,
        width = 500,
        xlabel = 'Raw Fluor. Data',
        cmap = color_vals,
        size=8,
        alpha = 1,
        show_legend = False,
        title = title)
    
    return a
    
#def plot_init_slopes():
    
    
    
    
#def plot_endpt_vals():
    