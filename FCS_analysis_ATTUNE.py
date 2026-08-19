# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:17:35 2024

@author: Eliza

Updated to use colored flow-cytometry-style filled density contour plots for gating figures.
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flowio import FlowData
from matplotlib.patches import Ellipse

SHOW_FIGS = True
PLOTTING_FIGS = True

# -----------------------------
# Global gating parameters
# -----------------------------
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

# Fluorescence thresholds on log10 scale
bl1_threshold = 4.0
yl2_threshold = 4.0

# Contour plotting controls
CONTOUR_BINS = 160
CONTOUR_LEVELS = 12
MIN_POINTS_FOR_CONTOUR = 200
SCATTER_FALLBACK_ALPHA = 0.25
SCATTER_FALLBACK_SIZE = 4
SAVE_DPI = 300


# -----------------------------
# Data loading and preprocessing
# -----------------------------
def load_fcs_data(file_path):
    fcs_data = FlowData(file_path, ignore_offset_error="True")
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


def preprocess_data(df, file_name=None):
    required_cols = ['FSC-H', 'SSC-H', 'FSC-A', 'BL1-H', 'YL2-H']
    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        print(f"Missing required columns in {file_name or 'file'}: {missing}")
        return pd.DataFrame(), {
            'total_events_loaded': len(df),
            'events_dropped_preprocess': len(df),
            'events_remaining_after_preprocess': 0,
            'drop_reason': f"Missing required columns: {missing}"
        }

    total_events_loaded = len(df)

    valid_mask = (
        (df['FSC-H'] > 0) &
        (df['SSC-H'] > 0) &
        (df['FSC-A'] > 0) &
        (df['BL1-H'] > 0) &
        (df['YL2-H'] > 0)
    )

    dropped_events = int((~valid_mask).sum())
    df = df.loc[valid_mask].copy()

    df['FSC-H (log)'] = np.log10(df['FSC-H'])
    df['SSC-H (log)'] = np.log10(df['SSC-H'])
    df['FSC-A (log)'] = np.log10(df['FSC-A'])
    df['BL1-H (log)'] = np.log10(df['BL1-H'])
    df['YL2-H (log)'] = np.log10(df['YL2-H'])

    finite_mask = np.isfinite(
        df[['FSC-H (log)', 'SSC-H (log)', 'FSC-A (log)', 'BL1-H (log)', 'YL2-H (log)']]
    ).all(axis=1)
    nonfinite_dropped = int((~finite_mask).sum())
    if nonfinite_dropped > 0:
        df = df.loc[finite_mask].copy()

    summary = {
        'total_events_loaded': total_events_loaded,
        'events_dropped_preprocess': dropped_events + nonfinite_dropped,
        'events_remaining_after_preprocess': len(df),
        'drop_reason': 'Non-positive or non-finite values in one or more required channels'
    }

    return df, summary


# -----------------------------
# Geometry helpers
# -----------------------------
def in_rotated_ellipse(x, y, center, width, height, angle):
    x_translated = x - center[0]
    y_translated = y - center[1]
    angle_rad = np.radians(-angle)
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad),  np.cos(angle_rad)]
    ])
    rotated_points = rotation_matrix @ np.vstack((x_translated, y_translated))
    x_rot, y_rot = rotated_points
    ellipse_eq = (x_rot ** 2 / (width / 2) ** 2) + (y_rot ** 2 / (height / 2) ** 2) <= 1
    return ellipse_eq


# -----------------------------
# Plotting helpers
# -----------------------------
def compute_density_grid(x, y, bins=CONTOUR_BINS):
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(xcenters, ycenters)
    Z = hist.T
    return X, Y, Z


def get_positive_contour_levels(Z, n_levels=CONTOUR_LEVELS):
    positive = Z[Z > 0]
    if positive.size == 0:
        return None

    zmin = positive.min()
    zmax = positive.max()

    if zmin == zmax:
        return np.array([zmin, zmax + 1e-9])

    return np.linspace(zmin, zmax, n_levels)


def plot_flow_contours(x, y, xlabel, ylabel, title, file_name, output_suffix,
                       ellipse_params=None, ellipse_edgecolor='r',
                       cmap='viridis', filled=True):
    plt.figure(figsize=(8, 6))

    if len(x) >= MIN_POINTS_FOR_CONTOUR:
        X, Y, Z = compute_density_grid(x, y, bins=CONTOUR_BINS)
        positive_mask = Z > 0

        if np.any(positive_mask):
            masked_Z = np.ma.masked_where(~positive_mask, Z)
            levels = get_positive_contour_levels(masked_Z.compressed(), n_levels=CONTOUR_LEVELS)

            if levels is not None and len(levels) > 0:
                if filled:
                    contourf = plt.contourf(X, Y, masked_Z, levels=levels, cmap=cmap, extend='max')
                    plt.contour(X, Y, masked_Z, levels=levels, colors='black', linewidths=0.35, alpha=0.35)
                    cbar = plt.colorbar(contourf)
                    cbar.set_label('Event density')
                else:
                    plt.contour(X, Y, masked_Z, levels=levels, cmap=cmap, linewidths=1.0)
            else:
                plt.scatter(x, y, alpha=SCATTER_FALLBACK_ALPHA, s=SCATTER_FALLBACK_SIZE, c='black', rasterized=True)
        else:
            plt.scatter(x, y, alpha=SCATTER_FALLBACK_ALPHA, s=SCATTER_FALLBACK_SIZE, c='black', rasterized=True)
    else:
        plt.scatter(x, y, alpha=SCATTER_FALLBACK_ALPHA, s=SCATTER_FALLBACK_SIZE, c='black', rasterized=True)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if ellipse_params is not None:
        ellipse = Ellipse(
            xy=ellipse_params['center'],
            width=ellipse_params['width'],
            height=ellipse_params['height'],
            angle=ellipse_params['angle'],
            edgecolor=ellipse_edgecolor,
            facecolor='none',
            linestyle='--',
            linewidth=2
        )
        plt.gca().add_patch(ellipse)

    plt.tight_layout()
    plt.savefig(f'{file_name}_{output_suffix}.png', dpi=SAVE_DPI)
    if SHOW_FIGS:
        plt.show()
    plt.close()


def plot_fsc_ssc_with_gate(df, file_name):
    plot_flow_contours(
        x=df['FSC-H (log)'].values,
        y=df['SSC-H (log)'].values,
        xlabel='FSC-H (log scale)',
        ylabel='SSC-H (log scale)',
        title=f'FSC-H vs SSC-H All Events ({file_name})',
        file_name=file_name,
        output_suffix='fsc_ssc',
        ellipse_params=CELL_ELLIPSE_PARAMS,
        ellipse_edgecolor='red',
        cmap='viridis',
        filled=True
    )


def apply_cell_gate(df):
    params = CELL_ELLIPSE_PARAMS
    cell_mask = in_rotated_ellipse(
        df['FSC-H (log)'].values,
        df['SSC-H (log)'].values,
        params['center'],
        params['width'],
        params['height'],
        params['angle']
    )
    return df[cell_mask].copy()


def plot_fsca_vs_fsch_with_gate(cell_df, file_name):
    plot_flow_contours(
        x=cell_df['FSC-A (log)'].values,
        y=cell_df['FSC-H (log)'].values,
        xlabel='FSC-A (log scale)',
        ylabel='FSC-H (log scale)',
        title=f'FSC-H vs FSC-A Gated on Cells ({file_name})',
        file_name=file_name,
        output_suffix='fsca_vs_fsch',
        ellipse_params=SINGLETS_ELLIPSE_PARAMS,
        ellipse_edgecolor='blue',
        cmap='plasma',
        filled=True
    )


def apply_singlet_gate(df):
    params = SINGLETS_ELLIPSE_PARAMS
    singlet_mask = in_rotated_ellipse(
        df['FSC-A (log)'].values,
        df['FSC-H (log)'].values,
        params['center'],
        params['width'],
        params['height'],
        params['angle']
    )
    return df[singlet_mask].copy()


def plot_histograms(df, gated_bl1_df, gated_yl2_df, bl1_threshold, yl2_threshold, file_name):
    plt.figure(figsize=(14, 7))

    plt.subplot(1, 2, 1)
    plt.hist(df['BL1-H (log)'].dropna(), bins=50, alpha=0.7, color='g', edgecolor='black')
    plt.axvline(
        bl1_threshold,
        color='r',
        linestyle='--',
        linewidth=1.5,
        label=f'Gating Threshold: {10 ** bl1_threshold:.0f}'
    )
    plt.xlabel('BL1-H (log scale)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of BL1-H (log) Gated on Singlets ({file_name})')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.hist(df['YL2-H (log)'].dropna(), bins=50, alpha=0.7, color='y', edgecolor='black')
    plt.axvline(
        yl2_threshold,
        color='r',
        linestyle='--',
        linewidth=1.5,
        label=f'Gating Threshold: {10 ** yl2_threshold:.0f}'
    )
    plt.xlabel('YL2-H (log scale)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of YL2-H (log) Gated on Singlets ({file_name})')
    plt.legend()

    plt.tight_layout()
    plt.savefig(f'{file_name}_histograms.png', dpi=SAVE_DPI)
    if SHOW_FIGS:
        plt.show()
    plt.close()


def calculate_percentage(gated_df, total_df):
    total_count = len(total_df)
    gated_count = len(gated_df)
    if total_count == 0:
        return np.nan
    percentage = (gated_count / total_count) * 100
    return percentage


# -----------------------------
# Batch processing
# -----------------------------
results = []
file_list = glob.glob('*.fcs')

for file_path in file_list:
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    print(f"Processing file: {file_name}")

    raw_df = load_fcs_data(file_path)
    df, preprocess_summary = preprocess_data(raw_df, file_name=file_name)

    if df.empty:
        print(f"Skipping {file_name}: no valid events after preprocessing.")
        results.append({
            'File Name': file_name,
            'Total Events Loaded': preprocess_summary['total_events_loaded'],
            'Dropped Events During Preprocessing': preprocess_summary['events_dropped_preprocess'],
            'Events Remaining After Preprocessing': preprocess_summary['events_remaining_after_preprocess'],
            'Cells After Gate': 0,
            'Singlets After Gate': 0,
            'BL1 Positive Events': 0,
            'YL2 Positive Events': 0,
            'Percentage BL1-H': 'Skipped',
            'Percentage YL2-H': 'Skipped',
            'Reason': preprocess_summary['drop_reason']
        })
        continue

    if PLOTTING_FIGS:
        plot_fsc_ssc_with_gate(df, file_name)

    cell_df = apply_cell_gate(df)

    if PLOTTING_FIGS:
        plot_fsca_vs_fsch_with_gate(cell_df, file_name)

    singlet_df = apply_singlet_gate(cell_df)

    gated_bl1_df = singlet_df[singlet_df['BL1-H (log)'] > bl1_threshold]
    gated_yl2_df = singlet_df[singlet_df['YL2-H (log)'] > yl2_threshold]

    if PLOTTING_FIGS:
        plot_histograms(singlet_df, gated_bl1_df, gated_yl2_df, bl1_threshold, yl2_threshold, file_name)

    percentage_bl1 = calculate_percentage(gated_bl1_df, singlet_df)
    percentage_yl2 = calculate_percentage(gated_yl2_df, singlet_df)

    results.append({
        'File Name': file_name,
        'Total Events Loaded': preprocess_summary['total_events_loaded'],
        'Dropped Events During Preprocessing': preprocess_summary['events_dropped_preprocess'],
        'Events Remaining After Preprocessing': preprocess_summary['events_remaining_after_preprocess'],
        'Cells After Gate': len(cell_df),
        'Singlets After Gate': len(singlet_df),
        'BL1 Positive Events': len(gated_bl1_df),
        'YL2 Positive Events': len(gated_yl2_df),
        'Percentage BL1-H': percentage_bl1,
        'Percentage YL2-H': percentage_yl2,
        'Reason': 'Processed Successfully'
    })

results_df = pd.DataFrame(results)

with pd.ExcelWriter('gating_results.xlsx') as writer:
    results_df.to_excel(writer, index=False, sheet_name='Summary')

print("Results have been saved to 'gating_results.xlsx'.")
