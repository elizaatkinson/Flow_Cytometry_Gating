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

# Define global ellipse parameters
CELL_ELLIPSE_PARAMS = {
    'center': (4.7, 4.7),
    'width': 1.8,
    'height': 0.6,
    'angle': 40
}

SINGLETS_ELLIPSE_PARAMS = {
    'center': (4.5, 4.5),
    'width': 1.5,
    'height': 0.1,
    'angle': 42
}

def load_fcs_data(file_path):
    fcs_data = FlowData(file_path)
    meta_data = fcs_data.text
    raw_data = np.array(fcs_data.events)
    n_channels = fcs_data.channel_count
    reshaped_data = np.reshape(raw_data, (-1, n_channels))

    def get_channel_name(meta_data, i):
        name_key = f'p{i}n'
        if name_key in meta_data and meta_data[name_key] != 'NA':
            return meta_data[name_key]
        return f'Channel_{i}'

    channel_names = [get_channel_name(meta_data, i) for i in range(1, n_channels + 1)]
    return pd.DataFrame(reshaped_data, columns=channel_names)

def preprocess_data(df):
    if 'FSC-H' in df.columns and 'SSC-H' in df.columns and 'FSC-A' in df.columns:
        df['FSC-H (log)'] = np.log10(df['FSC-H'])
        df['SSC-H (log)'] = np.log10(df['SSC-H'])
        df['FSC-A (log)'] = np.log10(df['FSC-A'])
        
        if 'BL1-H' in df.columns and 'YL2-H' in df.columns:
            df['BL1-H (log)'] = np.log10(df['BL1-H'])
            df['YL2-H (log)'] = np.log10(df['YL2-H'])
            df = df[np.isfinite(df['BL1-H (log)']) & np.isfinite(df['YL2-H (log)'])] 
        else:
            print("'BL1-H' or 'YL2-H' columns are missing from the DataFrame.")
    else:
        print("FSC-H, SSC-H, or FSC-A columns are missing from the DataFrame.")
    return df

def plot_fsc_ssc_with_gate(df, file_name):
    params = CELL_ELLIPSE_PARAMS
    plt.figure(figsize=(8, 6))
    plt.scatter(df['FSC-H (log)'], df['SSC-H (log)'], alpha=0.5, s=10)
    plt.xlabel('FSC-H (log scale)')
    plt.ylabel('SSC-H (log scale)')
    plt.title(f'FSC-H vs SSC-H All Events ({file_name})')

    ellipse = Ellipse(xy=params['center'], width=params['width'], height=params['height'], angle=params['angle'],
                      edgecolor='r', facecolor='none', linestyle='--', linewidth=2)
    plt.gca().add_patch(ellipse)
    plt.savefig(f'{file_name}_fsc_ssc.png')
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

    cell_mask = in_ellipse(df['FSC-H (log)'].values, df['SSC-H (log)'].values, params['center'], params['width'], params['height'], params['angle'])
    return df[cell_mask]

def plot_fsca_vs_fsch_with_gate(cell_df, file_name):
    params = SINGLETS_ELLIPSE_PARAMS

    plt.figure(figsize=(8, 6))
    plt.scatter(cell_df['FSC-A (log)'], cell_df['FSC-H (log)'], alpha=0.5, color='r', s=10)
    plt.xlabel('FSC-A (log scale)')
    plt.ylabel('FSC-H (log scale)')
    plt.title(f'FSC-H vs FSC-A Gated on Cells ({file_name})')
    
    ellipse = Ellipse(xy=params['center'], width=params['width'], height=params['height'], angle=params['angle'],
                      edgecolor='b', facecolor='none', linestyle='--', linewidth=2)
    plt.gca().add_patch(ellipse)
    plt.savefig(f'{file_name}_fsca_vs_fsch.png')
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

    singlet_mask = in_ellipse_fsca_fsch(df['FSC-A (log)'].values, df['FSC-H (log)'].values, params['center'], params['width'], params['height'], params['angle'])
    return df[singlet_mask]

def plot_histograms(df, gated_bl1_df, gated_yl2_df, bl1_threshold, yl2_threshold, file_name):
    plt.figure(figsize=(14, 7))

    plt.subplot(1, 2, 1)
    plt.hist(df['BL1-H (log)'].dropna(), bins=50, alpha=0.7, color='g', edgecolor='black')
    plt.axvline(bl1_threshold, color='r', linestyle='--', linewidth=1.5, label=f'Gating Threshold: {10**bl1_threshold:.0f}')
    plt.xlabel('BL1-H (log scale)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of BL1-H (log) Gated on Singlets ({file_name})')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.hist(df['YL2-H (log)'].dropna(), bins=50, alpha=0.7, color='y', edgecolor='black')
    plt.axvline(yl2_threshold, color='r', linestyle='--', linewidth=1.5, label=f'Gating Threshold: {10**yl2_threshold:.0f}')
    plt.xlabel('YL2-H (log scale)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of YL2-H (log) Gated on Singlets ({file_name})')
    plt.legend()

    plt.tight_layout()
    plt.savefig(f'{file_name}_histograms.png')
    plt.show()
    plt.close()

def calculate_percentage(gated_df, total_df):
    total_count = len(total_df)
    gated_count = len(gated_df)
    percentage = (gated_count / total_count) * 100
    return percentage

# Define gates
bl1_threshold = 4.0  # log10(10^4.5)
yl2_threshold = 4.0  # log10(10^4)

# Prepare a DataFrame to store results
results = []

# Get list of all FCS files in the current directory
file_list = glob.glob('*.fcs')

for file_path in file_list:
    file_name = os.path.splitext(os.path.basename(file_path))[0]  # Get file name without extension
    print(f"Processing file: {file_name}")
    
    df = load_fcs_data(file_path)
    df = preprocess_data(df)
    
    # Plot FSC vs SSC
    plot_fsc_ssc_with_gate(df, file_name)
    
    # Filter only cells
    cell_df = apply_cell_gate(df)
    
    # Plot FSC-A vs FSC-H
    plot_fsca_vs_fsch_with_gate(cell_df, file_name)
    
    # Filter out singlets
    singlet_df = apply_singlet_gate(cell_df)
    
    # Apply BL1-H and YL2-H gating
    gated_bl1_df = singlet_df[singlet_df['BL1-H (log)'] > bl1_threshold]
    gated_yl2_df = singlet_df[singlet_df['YL2-H (log)'] > yl2_threshold]
    
    # Plot histograms
    plot_histograms(df, gated_bl1_df, gated_yl2_df, bl1_threshold, yl2_threshold, file_name)
    
    
    # Calculate percentages
    percentage_bl1 = calculate_percentage(gated_bl1_df, singlet_df)
    percentage_yl2 = calculate_percentage(gated_yl2_df, singlet_df)
    
    # Append results to the list
    results.append({
        'File Name': file_name,
        'Percentage BL1-H': percentage_bl1,
        'Percentage YL2-H': percentage_yl2
    })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Save results to an Excel file
with pd.ExcelWriter('gating_results.xlsx') as writer:
    results_df.to_excel(writer, index=False)

print("Results have been saved to 'gating_results.xlsx'.")



