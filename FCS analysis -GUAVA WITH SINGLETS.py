# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:17:35 2024

@author: Eliza
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from flowio import FlowData
import pandas as pd
import glob
import os
import re  # Import regex to extract date and time from filenames

# Define global ellipse parameters
CELL_ELLIPSE_PARAMS = {
    'center': (3.7, 3.3),
    'width': 1.9,
    'height': 1.5,
    'angle': 50
}

SINGLETS_ELLIPSE_PARAMS = {
    'center': (3.55, 3.5),
    'width': 1.7,
    'height': 0.1,
    'angle': 45
}

def load_fcs_data(file_path):
    fcs_data = FlowData(file_path, ignore_offset_error="True")
    meta_data = fcs_data.text
    raw_data = np.array(fcs_data.events)
    n_channels = fcs_data.channel_count
    reshaped_data = np.reshape(raw_data, (-1, n_channels))

    # Extract the wellid from the meta_data
    sample_id = meta_data.get('wellid', 'Unknown')

    def get_channel_name(meta_data, i):
        name_key = f'p{i}n'
        if name_key in meta_data and meta_data[name_key] != 'NA':
            return meta_data[name_key]
        return f'Channel_{i}'

    channel_names = [get_channel_name(meta_data, i) for i in range(1, n_channels + 1)]
    
    # Return the DataFrame and the sample ID (wellid)
    return pd.DataFrame(reshaped_data, columns=channel_names), sample_id

def preprocess_data(df):
    if 'FSC-HLin' in df.columns and 'SSC-HLin' in df.columns and 'FSC-ALin' in df.columns:
        df['FSC-HLin (log)'] = np.log10(df['FSC-HLin'])
        df['SSC-HLin (log)'] = np.log10(df['SSC-HLin'])
        df['FSC-ALin (log)'] = np.log10(df['FSC-ALin'])
        
        if 'GRN-B-HLin' in df.columns and 'RED-G-HLin' in df.columns:
            df['GRN-B-HLin (log)'] = np.log10(df['GRN-B-HLin'])
            df['RED-G-HLin (log)'] = np.log10(df['RED-G-HLin'])
            df = df[np.isfinite(df['GRN-B-HLin (log)']) & np.isfinite(df['RED-G-HLin (log)'])] 
        else:
            print("'GRN-B-HLin' or 'RED-G-HLin' columns are missing from the DataFrame.")
    else:
        print("FSC-HLin, SSC-HLin, or FSC-ALin columns are missing from the DataFrame.")
    return df

def plot_fsc_ssc_with_gate(df, file_name):
    params = CELL_ELLIPSE_PARAMS
    plt.figure(figsize=(8, 6))
    plt.scatter(df['FSC-HLin (log)'], df['SSC-HLin (log)'], alpha=0.5, s=10)
    plt.xlabel('FSC-HLin (log scale)')
    plt.ylabel('SSC-HLin (log scale)')
    plt.title(f'FSC-HLin vs SSC-HLin All Events ({file_name})')

    ellipse = Ellipse(xy=params['center'], width=params['width'], height=params['height'], angle=params['angle'],
                      edgecolor='r', facecolor='none', linestyle='--', linewidth=2)
    plt.gca().add_patch(ellipse)
    plt.savefig(f'{file_name}_fs_sc.png')
    plt.show()
    plt.close()

def apply_cell_gate(df):
    params = CELL_ELLIPSE_PARAMS

    def in_ellipse(x, y, center, width, height, angle):
        x_translated = x - center[0]
        y_translated = y - center[1]
        angle_rad = np.radians(-angle)
        rotation_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])
        rotated_points = rotation_matrix @ np.vstack((x_translated, y_translated))
        x_rot, y_rot = rotated_points
        ellipse_eq = (x_rot**2 / (width / 2)**2) + (y_rot**2 / (height / 2)**2) <= 1
        return ellipse_eq

    cell_mask = in_ellipse(df['FSC-HLin (log)'].values, df['SSC-HLin (log)'].values, params['center'], params['width'], params['height'], params['angle'])
    return df[cell_mask]

def plot_fsca_vs_fsch_with_gate(cell_df, file_name):
    params = SINGLETS_ELLIPSE_PARAMS

    plt.figure(figsize=(8, 6))
    plt.scatter(cell_df['FSC-ALin (log)'], cell_df['FSC-HLin (log)'], alpha=0.5, color='r', s=10)
    plt.xlabel('FSC-ALin (log scale)')
    plt.ylabel('FSC-HLin (log scale)')
    plt.title(f'FSC-HLin vs FSC-ALin Gated on Cells ({file_name})')
    
    ellipse = Ellipse(xy=params['center'], width=params['width'], height=params['height'], angle=params['angle'],
                      edgecolor='b', facecolor='none', linestyle='--', linewidth=2)
    plt.gca().add_patch(ellipse)
    plt.savefig(f'{file_name}_fsca_h.png')
    plt.show()
    plt.close()

def apply_singlet_gate(df):
    params = SINGLETS_ELLIPSE_PARAMS

    def in_ellipse_fsca_fsch(x, y, center, width, height, angle):
        x_translated = x - center[0]
        y_translated = y - center[1]
        angle_rad = np.radians(-angle)
        rotation_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])
        rotated_points = rotation_matrix @ np.vstack((x_translated, y_translated))
        x_rot, y_rot = rotated_points
        ellipse_eq = (x_rot**2 / (width / 2)**2) + (y_rot**2 / (height / 2)**2) <= 1
        return ellipse_eq

    singlet_mask = in_ellipse_fsca_fsch(df['FSC-ALin (log)'].values, df['FSC-HLin (log)'].values, params['center'], params['width'], params['height'], params['angle'])
    return df[singlet_mask]

def plot_histograms(df, gated_grn_df, gated_red_df, grn_threshold, red_threshold, file_name):
    plt.figure(figsize=(14, 7))

    plt.subplot(1, 2, 1)
    plt.hist(df['GRN-B-HLin (log)'].dropna(), bins=50, alpha=0.7, color='g', edgecolor='black')
    plt.axvline(grn_threshold, color='r', linestyle='--', linewidth=1.5, label=f'Gating Threshold: {10**grn_threshold:.0f}')
    plt.xlabel('GRN-B-HLin (log scale)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of GRN-B-HLin (log) Gated on Singlets ({file_name})')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.hist(df['RED-G-HLin (log)'].dropna(), bins=50, alpha=0.7, color='y', edgecolor='black')
    plt.axvline(red_threshold, color='r', linestyle='--', linewidth=1.5, label=f'Gating Threshold: {10**red_threshold:.0f}')
    plt.xlabel('RED-G-HLin (log scale)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of RED-G-HLin (log) Gated on Singlets ({file_name})')
    plt.legend()

    plt.tight_layout()
    plt.savefig(f'{file_name}_hist.png')
    plt.show()
    plt.close()

def calculate_percentage(gated_df, total_df):
    total_count = len(total_df)
    gated_count = len(gated_df)
    percentage = (gated_count / total_count) * 100
    return percentage

# Define gates
grn_threshold = 2.4  # log10(10^4.5)
red_threshold = 1.6  # log10(10^4)

# Prepare to store results
results = []

# Get a list of FCS files to process
fcs_files = glob.glob('*.fcs')
last_file_date = None
last_file_time = None

for file_path in fcs_files:
    file_name = os.path.basename(file_path)
    
    # Extract date and time from the filename using regex
    match = re.search(r'(\d{4}-\d{2}-\d{2})_at_(\d{2}-\d{2}-\d{2}(?:am|pm))', file_name)
    if match:
        file_date = match.group(1)  # Extract the date part
        file_time = match.group(2)  # Extract the time part
    else:
        file_date = 'Unknown'
        file_time = 'Unknown'
    
    # Update the last processed file's date and time
    last_file_date = file_date
    last_file_time = file_time

    # Load and process FCS data
    df, sample_id = load_fcs_data(file_path)
    df = preprocess_data(df)
    
    # Plot FSC-HLin vs SSC-HLin
    #plot_fsc_ssc_with_gate(df, file_name)
    
    # Filter only cells
    cell_df = apply_cell_gate(df)
    
    # Check if there are fewer than 1000 events in the cell gate
    if len(cell_df) < 1000:
        print(f"Skipping {file_name}: fewer than 1000 cells after gating.")
        results.append({
            'File Name': file_name,
            'Sample ID': sample_id,
            'Date': file_date,
            'Time': file_time,
            'Percentage GRN-B-HLin': 'Skipped',
            'Percentage RED-G-HLin': 'Skipped',
            'Reason': 'Fewer than 1000 cells after gating'
        })
        continue  # Skip this file and move to the next one
        
    # Plot FSC-ALin vs FSC-HLin
    #plot_fsca_vs_fsch_with_gate(cell_df, file_name)
    
    # Filter out singlets
    singlet_df = apply_singlet_gate(cell_df)
    
    # Apply GRN-B-HLin and RED-G-HLin gating
    gated_grn_df = singlet_df[singlet_df['GRN-B-HLin (log)'] > grn_threshold]
    gated_red_df = singlet_df[singlet_df['RED-G-HLin (log)'] > red_threshold]
    
    # Plot histograms
    #plot_histograms(df, gated_grn_df, gated_red_df, grn_threshold, red_threshold, file_name)
    
    # Calculate percentages
    percentage_grn = calculate_percentage(gated_grn_df, singlet_df)
    percentage_red = calculate_percentage(gated_red_df, singlet_df)
    
    # Append results to the list, including date and time
    results.append({
        'File Name': file_name,
        'Sample ID': sample_id,
        'Date': file_date,
        'Time': file_time,
        'Percentage GRN-B-HLin': percentage_grn,
        'Percentage RED-G-HLin': percentage_red,
        'Reason': 'Processed Successfully'
    })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Save results to an Excel file with the last file's date and time in the filename
output_filename = f'gating_results_{last_file_date}_{last_file_time}.xlsx'
output_filename = output_filename.replace(":", "-")  # Replace colons in time with dashes for filename compatibility

with pd.ExcelWriter(output_filename) as writer:
    results_df.to_excel(writer, index=False)

print(f"Results have been saved to '{output_filename}'.")




