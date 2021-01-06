import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()
import pandas as pd
import murraylab_tools.biotek as mt_biotek
import os
from get_colors import get_colors

def y_vals(df, well):
    '''Used in plot_timeseries to extract desired data'''
    df2 = df.loc[df['Well'] == well, :]
    return df2['Measurement'].values

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
    inputfolder = data_folder
    outputfolder = data_folder
    data_name = data_name
    data_filename = os.path.join(inputfolder, data_name+".csv")
    tidy_filename = os.path.join(outputfolder,data_name+"_tidy.csv")
    mt_biotek.tidy_biotek_data(data_filename, volume = 10,
                               output_filename = tidy_filename)
    df = pd.read_csv(tidy_filename)
    
    # Get subsection of gain
    df_gain = df.loc[df['Gain'] == gain]
    
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
    inputfolder_pre = data_folder
    outputfolder_pre = data_folder
    data_name_pre = data_name_pre
    data_filename_pre = os.path.join(inputfolder_pre, data_name_pre+".csv")
    tidy_filename_pre = os.path.join(outputfolder_pre,data_name_pre+"_tidy.csv")
    mt_biotek.tidy_biotek_data(data_filename_pre, volume = 10,
                               output_filename = tidy_filename_pre)
    df_pre = pd.read_csv(tidy_filename_pre)

    # Get subsection of gain and add time 2 column
    df_gain_pre = df_pre.loc[df_pre['Gain'] == gain]
    df_gain_pre['Time2'] = df_gain_pre['Time (hr)']
    
    # POSTPSIKE DATA
    # Get in tidy format with murraylabtools
    inputfolder_post = data_folder
    outputfolder_post = data_folder
    data_name_post = data_name_post
    data_filename_post = os.path.join(inputfolder_post, data_name_post+".csv")
    tidy_filename_post = os.path.join(outputfolder_post,data_name_post+"_tidy.csv")
    mt_biotek.tidy_biotek_data(data_filename_post, volume = 10,
                               output_filename = tidy_filename_post)
    df_post = pd.read_csv(tidy_filename_post)
    
    # Get subsection of gain
    df_gain_post = df_post.loc[df_post['Gain'] == gain]
    
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

    
#def plot_init_slopes():
    
    
    
    
#def plot_endpt_vals():
    