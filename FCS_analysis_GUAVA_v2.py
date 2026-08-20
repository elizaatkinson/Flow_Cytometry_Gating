# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 11:17:35 2024

@author: Eliza

"""

import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flowio import FlowData
from matplotlib.patches import Ellipse

PLOTTING_FIGS = True  # set to True to save gating figures
SHOW_FIGS = True

# -----------------------------
# Global gating parameters
# -----------------------------
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

# Define histogram gates
# 10**2.4 ≈ 251
# 10**1.6 ≈ 40
grn_threshold = 2.4
red_threshold = 1.6

# Standardized axis limits for scatter plots
FSC_SSC_XLIM = (2.0, 5.5)
FSC_SSC_YLIM = (2.0, 4.8)
FSCA_FSCH_XLIM = (2.5, 5.0)
FSCA_FSCH_YLIM = (2.5, 5.0)

# Standardized axis limits for histograms
GRN_HIST_XLIM = (0.0, 4.0)
RED_HIST_XLIM = (0.0, 3.0)

# Contour plotting controls
CONTOUR_BINS = 160
CONTOUR_LEVELS = 12
MIN_POINTS_FOR_CONTOUR = 200
SCATTER_FALLBACK_ALPHA = 0.25
SCATTER_FALLBACK_SIZE = 4
SAVE_DPI = 300
HIST_BINS = 50

# Label box styling
LABEL_BOX_KWARGS = {
    'boxstyle': 'round,pad=0.3',
    'facecolor': 'white',
    'edgecolor': 'black',
    'linewidth': 0.8,
    'alpha': 0.9
}


# -----------------------------
# Utility functions
# -----------------------------
def safe_filename(name):
    # Remove leading FCS3Exported_
    name = re.sub(r'^FCS3Exported_', '', name)

    # Convert patterns like:
    # 2023-11-09_at_11-36-46am-1
    # -> 231109_113646_1
    m = re.match(r'(\d{4})-(\d{2})-(\d{2})_at_(\d{2})-(\d{2})-(\d{2})(am|pm)-(\d+)', name)
    if m:
        year, month, day, hour, minute, second, ampm, suffix = m.groups()
        hour = int(hour)

        if ampm.lower() == 'pm' and hour != 12:
            hour += 12
        elif ampm.lower() == 'am' and hour == 12:
            hour = 0

        return f"{year[2:]}{month}{day}_{hour:02d}{minute}{second}_{suffix}"

    # Fallback: general cleanup
    name = re.sub(r'[<>:\"/\\|?*]', '_', name)
    name = re.sub(r'\s+', '_', name.strip())
    return name[:40].rstrip('._')


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
    required_cols = ['FSC-HLin', 'SSC-HLin', 'FSC-ALin', 'GRN-B-HLin', 'RED-G-HLin']
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
        (df['FSC-HLin'] > 0) &
        (df['SSC-HLin'] > 0) &
        (df['FSC-ALin'] > 0) &
        (df['GRN-B-HLin'] > 0) &
        (df['RED-G-HLin'] > 0)
    )

    dropped_events = int((~valid_mask).sum())
    df = df.loc[valid_mask].copy()

    df['FSC-HLin (log)'] = np.log10(df['FSC-HLin'])
    df['SSC-HLin (log)'] = np.log10(df['SSC-HLin'])
    df['FSC-ALin (log)'] = np.log10(df['FSC-ALin'])
    df['GRN-B-HLin (log)'] = np.log10(df['GRN-B-HLin'])
    df['RED-G-HLin (log)'] = np.log10(df['RED-G-HLin'])

    finite_mask = np.isfinite(
        df[['FSC-HLin (log)', 'SSC-HLin (log)', 'FSC-ALin (log)', 'GRN-B-HLin (log)', 'RED-G-HLin (log)']]
    ).all(axis=1)
    nonfinite_dropped = int((~finite_mask).sum())
    if nonfinite_dropped > 0:
        df = df.loc[finite_mask].copy()

    summary = {
        'total_events_loaded': total_events_loaded,
        'events_dropped_preprocess': dropped_events + nonfinite_dropped,
        'events_remaining_after_preprocess': len(df),
        'drop_reason': 'Non-positive or non-finite values in one or more log-transformed channels'
    }

    return df, summary


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


def add_boxed_label(ax, text, x=0.02, y=0.98, fontsize=10):
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        ha='left',
        va='top',
        fontsize=fontsize,
        bbox=LABEL_BOX_KWARGS,
        zorder=10
    )


def plot_flow_contours(x, y, xlabel, ylabel, title, file_name, output_suffix,
                       ellipse_params=None, ellipse_edgecolor='r',
                       cmap='viridis', filled=True, gate_fraction=None,
                       gate_label='In gate', xlim=None, ylim=None):
    fig, ax = plt.subplots(figsize=(8, 6))

    if len(x) >= MIN_POINTS_FOR_CONTOUR:
        X, Y, Z = compute_density_grid(x, y, bins=CONTOUR_BINS)
        positive_mask = Z > 0

        if np.any(positive_mask):
            masked_Z = np.ma.masked_where(~positive_mask, Z)
            levels = get_positive_contour_levels(masked_Z.compressed(), n_levels=CONTOUR_LEVELS)

            if levels is not None and len(levels) > 0:
                if filled:
                    contourf = ax.contourf(X, Y, masked_Z, levels=levels, cmap=cmap, extend='max')
                    ax.contour(X, Y, masked_Z, levels=levels, colors='black', linewidths=0.35, alpha=0.35)
                    cbar = fig.colorbar(contourf, ax=ax)
                    cbar.set_label('Event density')
                else:
                    ax.contour(X, Y, masked_Z, levels=levels, cmap=cmap, linewidths=1.0)
            else:
                ax.scatter(x, y, alpha=SCATTER_FALLBACK_ALPHA, s=SCATTER_FALLBACK_SIZE, c='black', rasterized=True)
        else:
            ax.scatter(x, y, alpha=SCATTER_FALLBACK_ALPHA, s=SCATTER_FALLBACK_SIZE, c='black', rasterized=True)
    else:
        ax.scatter(x, y, alpha=SCATTER_FALLBACK_ALPHA, s=SCATTER_FALLBACK_SIZE, c='black', rasterized=True)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

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
        ax.add_patch(ellipse)

    if gate_fraction is not None and not np.isnan(gate_fraction):
        add_boxed_label(ax, f'{gate_label}: {gate_fraction:.1f}%')

    fig.tight_layout()
    fig.savefig(f'{file_name}_{output_suffix}.png', dpi=SAVE_DPI)
    if SHOW_FIGS:
        plt.show()
    plt.close(fig)


def plot_fsc_ssc_with_gate(df, file_name, cell_fraction):
    plot_flow_contours(
        x=df['FSC-HLin (log)'].values,
        y=df['SSC-HLin (log)'].values,
        xlabel='FSC-HLin (log scale)',
        ylabel='SSC-HLin (log scale)',
        title=f'FSC-HLin vs SSC-HLin All Events ({file_name})',
        file_name=file_name,
        output_suffix='fsc_ssc',
        ellipse_params=CELL_ELLIPSE_PARAMS,
        ellipse_edgecolor='red',
        cmap='viridis',
        filled=True,
        gate_fraction=cell_fraction,
        gate_label='In cell gate',
        xlim=FSC_SSC_XLIM,
        ylim=FSC_SSC_YLIM
    )


def apply_cell_gate(df):
    params = CELL_ELLIPSE_PARAMS
    cell_mask = in_rotated_ellipse(
        df['FSC-HLin (log)'].values,
        df['SSC-HLin (log)'].values,
        params['center'],
        params['width'],
        params['height'],
        params['angle']
    )
    return df[cell_mask].copy()


def plot_fsca_vs_fsch_with_gate(cell_df, file_name, singlet_fraction):
    plot_flow_contours(
        x=cell_df['FSC-ALin (log)'].values,
        y=cell_df['FSC-HLin (log)'].values,
        xlabel='FSC-ALin (log scale)',
        ylabel='FSC-HLin (log scale)',
        title=f'FSC-HLin vs FSC-ALin Gated on Cells ({file_name})',
        file_name=file_name,
        output_suffix='fsca_vs_fsch',
        ellipse_params=SINGLETS_ELLIPSE_PARAMS,
        ellipse_edgecolor='blue',
        cmap='plasma',
        filled=True,
        gate_fraction=singlet_fraction,
        gate_label='In singlet gate',
        xlim=FSCA_FSCH_XLIM,
        ylim=FSCA_FSCH_YLIM
    )


def apply_singlet_gate(df):
    params = SINGLETS_ELLIPSE_PARAMS
    singlet_mask = in_rotated_ellipse(
        df['FSC-ALin (log)'].values,
        df['FSC-HLin (log)'].values,
        params['center'],
        params['width'],
        params['height'],
        params['angle']
    )
    return df[singlet_mask].copy()


def compute_histogram_ymax(values, xlim, bins=HIST_BINS):
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return 1

    counts, _ = np.histogram(finite_values, bins=bins, range=xlim)
    ymax = counts.max() if counts.size > 0 else 0
    return max(1, int(np.ceil(ymax * 1.1)))


def compute_global_histogram_limits(file_list):
    grn_ymax = 1
    red_ymax = 1

    for file_path in file_list:
        file_name = safe_filename(os.path.splitext(os.path.basename(file_path))[0])
        raw_df = load_fcs_data(file_path)
        df, _ = preprocess_data(raw_df, file_name=file_name)

        if df.empty:
            continue

        cell_df = apply_cell_gate(df)
        if cell_df.empty:
            continue

        singlet_df = apply_singlet_gate(cell_df)
        if singlet_df.empty:
            continue

        current_grn_ymax = compute_histogram_ymax(
            singlet_df['GRN-B-HLin (log)'].values,
            GRN_HIST_XLIM,
            bins=HIST_BINS
        )
        current_red_ymax = compute_histogram_ymax(
            singlet_df['RED-G-HLin (log)'].values,
            RED_HIST_XLIM,
            bins=HIST_BINS
        )

        grn_ymax = max(grn_ymax, current_grn_ymax)
        red_ymax = max(red_ymax, current_red_ymax)

    return grn_ymax, red_ymax


def plot_histograms(df, gated_grn_df, gated_red_df, grn_threshold, red_threshold, file_name,
                    percentage_grn, percentage_red, grn_hist_ylim, red_hist_ylim):
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    ax1 = axes[0]
    ax1.hist(df['GRN-B-HLin (log)'].dropna(), bins=HIST_BINS, range=GRN_HIST_XLIM,
             alpha=0.7, color='g', edgecolor='black')
    ax1.axvline(grn_threshold, color='r', linestyle='--', linewidth=1.5)
    ax1.set_xlabel('GRN-B-HLin (log scale)')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'Histogram of GRN-B-HLin (log) Gated on Singlets ({file_name})')
    ax1.set_xlim(GRN_HIST_XLIM)
    ax1.set_ylim(0, grn_hist_ylim)
    add_boxed_label(ax1, f'In gate: {percentage_grn:.1f}%', x=0.02, y=0.98)
    add_boxed_label(ax1, f'Threshold: {10 ** grn_threshold:.0f}', x=0.02, y=0.93)

    ax2 = axes[1]
    ax2.hist(df['RED-G-HLin (log)'].dropna(), bins=HIST_BINS, range=RED_HIST_XLIM,
             alpha=0.7, color='y', edgecolor='black')
    ax2.axvline(red_threshold, color='r', linestyle='--', linewidth=1.5)
    ax2.set_xlabel('RED-G-HLin (log scale)')
    ax2.set_ylabel('Frequency')
    ax2.set_title(f'Histogram of RED-G-HLin (log) Gated on Singlets ({file_name})')
    ax2.set_xlim(RED_HIST_XLIM)
    ax2.set_ylim(0, red_hist_ylim)
    add_boxed_label(ax2, f'In gate: {percentage_red:.1f}%', x=0.02, y=0.98)
    add_boxed_label(ax2, f'Threshold: {10 ** red_threshold:.0f}', x=0.02, y=0.93)

    fig.tight_layout()
    fig.savefig(f'{file_name}_histograms.png', dpi=SAVE_DPI)
    if SHOW_FIGS:
        plt.show()
    plt.close(fig)


def calculate_percentage(gated_df, total_df):
    total_count = len(total_df)
    gated_count = len(gated_df)
    if total_count == 0:
        return np.nan
    percentage = (gated_count / total_count) * 100
    return percentage


# -----------------------------
# Main processing
# -----------------------------
results = []
file_list = glob.glob('*.fcs')

# Compute global standardized histogram y-limits across all files
GRN_HIST_YMAX, RED_HIST_YMAX = compute_global_histogram_limits(file_list)

for file_path in file_list:
    file_name = safe_filename(os.path.splitext(os.path.basename(file_path))[0])
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
            'GRN Positive Events': 0,
            'RED Positive Events': 0,
            'Percentage GRN-B-HLin': 'Skipped',
            'Percentage RED-G-HLin': 'Skipped',
            'Reason': preprocess_summary['drop_reason']
        })
        continue

    cell_df = apply_cell_gate(df)
    cell_fraction = calculate_percentage(cell_df, df)

    if PLOTTING_FIGS:
        plot_fsc_ssc_with_gate(df, file_name, cell_fraction)

    if len(cell_df) < 1000:
        print(f"Skipping {file_name}: fewer than 1000 cells after gating.")
        results.append({
            'File Name': file_name,
            'Total Events Loaded': preprocess_summary['total_events_loaded'],
            'Dropped Events During Preprocessing': preprocess_summary['events_dropped_preprocess'],
            'Events Remaining After Preprocessing': preprocess_summary['events_remaining_after_preprocess'],
            'Cells After Gate': len(cell_df),
            'Singlets After Gate': 0,
            'GRN Positive Events': 0,
            'RED Positive Events': 0,
            'Percentage GRN-B-HLin': 'Skipped',
            'Percentage RED-G-HLin': 'Skipped',
            'Reason': 'Fewer than 1000 cells after gating'
        })
        continue

    singlet_df = apply_singlet_gate(cell_df)
    singlet_fraction = calculate_percentage(singlet_df, cell_df)

    if PLOTTING_FIGS:
        plot_fsca_vs_fsch_with_gate(cell_df, file_name, singlet_fraction)

    gated_grn_df = singlet_df[singlet_df['GRN-B-HLin (log)'] > grn_threshold]
    gated_red_df = singlet_df[singlet_df['RED-G-HLin (log)'] > red_threshold]

    percentage_grn = calculate_percentage(gated_grn_df, singlet_df)
    percentage_red = calculate_percentage(gated_red_df, singlet_df)

    if PLOTTING_FIGS:
        plot_histograms(
            singlet_df,
            gated_grn_df,
            gated_red_df,
            grn_threshold,
            red_threshold,
            file_name,
            percentage_grn,
            percentage_red,
            GRN_HIST_YMAX,
            RED_HIST_YMAX
        )

    results.append({
        'File Name': file_name,
        'Total Events Loaded': preprocess_summary['total_events_loaded'],
        'Dropped Events During Preprocessing': preprocess_summary['events_dropped_preprocess'],
        'Events Remaining After Preprocessing': preprocess_summary['events_remaining_after_preprocess'],
        'Cells After Gate': len(cell_df),
        'Singlets After Gate': len(singlet_df),
        'GRN Positive Events': len(gated_grn_df),
        'RED Positive Events': len(gated_red_df),
        'Percentage GRN-B-HLin': percentage_grn,
        'Percentage RED-G-HLin': percentage_red,
        'Reason': 'Processed Successfully'
    })

results_df = pd.DataFrame(results)

with pd.ExcelWriter('gating_results.xlsx') as writer:
    results_df.to_excel(writer, index=False, sheet_name='Summary')

print("Results have been saved to 'gating_results.xlsx'.")
