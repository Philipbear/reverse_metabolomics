import pandas as pd
import numpy as np
from .raw_data_utils import MSData
from .main import init_config
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_all_eic(df, file_path, out_dir, plots_per_page=5):
    """
    Plot EICs for adducts in the dataframe and additional related adducts on a PDF with an M x 1 layout per page and legend inside the plot

    Args:
        df (pd.DataFrame): DataFrame containing the data
        file_path (str): Path to the raw data file
        out_dir (str): Output directory for the PDF
        plots_per_page (int): Number of plots per page (default: 5)
    """
    config = init_config(mass_detect_int_tol=0)
    d = MSData()
    d.read_raw_data(file_path, config)

    file_name = file_path.split('/')[-1].split('.')[0]

    df = df[~pd.isnull(df['best_MS2_scan_idx'])].reset_index(drop=True)
    df = df[~df['discard_reason'].str.contains('MS2 explanation below cutoff')].reset_index(drop=True)

    _adduct_pos = {
        '[M+H]+': {'m': 1, 'charge': 1, 'mass': 1.00727645223},
        '[M+Na]+': {'m': 1, 'charge': 1, 'mass': 22.989220702},
        '[M+K]+': {'m': 1, 'charge': 1, 'mass': 38.9631579064},
        '[M+NH4]+': {'m': 1, 'charge': 1, 'mass': 18.03382555335},
        '[M-H2O+NH4]+': {'m': 1, 'charge': 1, 'mass': 0.023260869},
        '[M+H-H2O]+': {'m': 1, 'charge': 1, 'mass': -17.0032882318},
        '[M+H-2H2O]+': {'m': 1, 'charge': 1, 'mass': -35.01385291583},
        '[M+H-3H2O]+': {'m': 1, 'charge': 1, 'mass': -53.02441759986}
    }

    def get_adducts_to_plot(adducts):
        adducts_to_plot = set()
        for adduct in adducts:
            if '[M+H]+' in adduct:
                adducts_to_plot.update(['[M+H]+', '[M+H-H2O]+'])
            if '[M+H-H2O]+' in adduct:
                adducts_to_plot.update(['[M+H]+', '[M+H-H2O]+', '[M+H-2H2O]+'])
            if '[M+H-2H2O]+' in adduct:
                adducts_to_plot.update(['[M+H]+', '[M+H-H2O]+', '[M+H-2H2O]+'])
            if '[M+H-3H2O]+' in adduct:
                adducts_to_plot.update(['[M+H]+', '[M+H-H2O]+', '[M+H-2H2O]+', '[M+H-3H2O]+'])
            if '[M+NH4]+' in adduct:
                adducts_to_plot.update(['[M+H]+', '[M+NH4]+'])
            adducts_to_plot.add(adduct)  # Always include the original adduct

        # For all other cases, always include [M+H]+ and [M+H-H2O]+
        adducts_to_plot.update(['[M+H]+', '[M+H-H2O]+'])

        return list(adducts_to_plot)

    pdf_path = f"{out_dir}/{file_name}_eic_plots.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, axs = plt.subplots(plots_per_page, 1, figsize=(8.27, 11.69))  # A4 size
        if plots_per_page == 1:
            axs = [axs]  # Make axs always iterable
        fig.suptitle(f"EIC Plots - {file_name}", fontsize=14)
        plot_count = 0

        for compound_name, group in df.groupby('compound_name'):
            if plot_count == plots_per_page:
                # Save the current page and start a new one
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig(fig)
                plt.close(fig)

                # Start a new page
                fig, axs = plt.subplots(plots_per_page, 1, figsize=(8.27, 11.69))
                if plots_per_page == 1:
                    axs = [axs]
                fig.suptitle(f"EIC Plots - {file_name}", fontsize=14)
                plot_count = 0

            exact_mass = group['exact_mass'].iloc[0]

            ax = axs[plot_count]
            ax.set_title(f"{compound_name}", fontsize=10)
            ax.set_xlabel('RT (min)', fontsize=8)
            ax.set_ylabel('Intensity', fontsize=8)

            eic_data = {}
            adducts_to_plot = get_adducts_to_plot(group['t_adduct'].unique())

            for t_adduct in adducts_to_plot:
                if t_adduct in _adduct_pos:
                    t_mz = (exact_mass + _adduct_pos[t_adduct]['mass']) / _adduct_pos[t_adduct]['charge']
                    eic_rt, eic_int, _, _ = d.get_eic_data(t_mz, mz_tol=max(0.005, t_mz * 1e-5))
                    line, = ax.plot(eic_rt, eic_int, label=f"{t_adduct} (m/z {t_mz:.4f})", linewidth=0.8, alpha=0.7)
                    eic_data[t_adduct] = {'color': line.get_color(), 'rt': eic_rt, 'int': eic_int}

            for _, row in group.iterrows():
                ms2_rt = d.scans[int(row['best_MS2_scan_idx']) - 1].rt
                if row['t_adduct'] in eic_data:
                    color = eic_data[row['t_adduct']]['color']
                    eic_rt = eic_data[row['t_adduct']]['rt']
                    eic_int = eic_data[row['t_adduct']]['int']

                    prev_rt_index = np.searchsorted(eic_rt, ms2_rt) - 1
                    if prev_rt_index >= 0:
                        prev_rt = eic_rt[prev_rt_index]
                        intensity = eic_int[prev_rt_index]

                        # Plot dots for MS2 scan
                        ax.plot(ms2_rt, intensity, 'o', color=color, markersize=3)

                        # Add label for RT
                        ax.text(ms2_rt, intensity, f'{prev_rt:.2f}',
                                color=color, fontsize=5, ha='right', va='bottom',
                                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0))

                        # Add cross if not selected
                        if not row.get('selected', True):  # Default to True if 'selected' column doesn't exist
                            ax.plot(ms2_rt, intensity, 'x', color='black', markersize=4.5, linewidth=0.4)

            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.legend(fontsize=6, loc='upper right', bbox_to_anchor=(1, 1), ncol=1)
            plot_count += 1

        # Save the last page
        if plot_count > 0:
            for i in range(plot_count, plots_per_page):
                axs[i].axis('off')
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"EIC plots saved to {pdf_path}")


def plot_all_ms2(df, file_path, out_dir):
    """
    Plot all MS2 spectra on a PDF with a 5x2 layout per page, including unselected spectra
    """
    config = init_config(mass_detect_int_tol=0)
    d = MSData()
    d.read_raw_data(file_path, config)

    file_name = file_path.split('/')[-1].split('.')[0]

    pdf_path = f"{out_dir}/{file_name}_ms2_spectra.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, axs = plt.subplots(5, 2, figsize=(11, 17))  # A4 size in inches
        fig.suptitle("MS2 Spectra", fontsize=16)
        plot_count = 0

        for i, row in df.iterrows():
            if not pd.isnull(row['best_MS2_scan_idx']):
                this_scan = d.scans[int(row['best_MS2_scan_idx']) - 1]
                if this_scan is None or this_scan.peaks is None:
                    continue
                mz_arr, int_arr = this_scan.peaks[:, 0], this_scan.peaks[:, 1]

                ax = axs[plot_count // 2, plot_count % 2]

                # Create selection status and discard reason
                status = "Selected" if row['selected'] else "Unselected"
                discard_reason = row['discard_reason'] if not pd.isnull(row['discard_reason']) else None

                plot_single_ms2(ax, mz_arr, int_arr, row['compound_name'], row['t_adduct'],
                                this_scan.precursor_mz, row['RT'], row['best_MS2_scan_idx'],
                                row['ms2_explained_intensity'], status, discard_reason)

                plot_count += 1

                if plot_count == 10:
                    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                    pdf.savefig(fig)
                    plt.close(fig)
                    fig, axs = plt.subplots(5, 2, figsize=(11, 17))
                    fig.suptitle("MS2 Spectra", fontsize=16)
                    plot_count = 0

        # Save the last page if there are any remaining plots
        if plot_count > 0:
            for i in range(plot_count, 10):
                axs[i // 2, i % 2].axis('off')
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"MS2 spectra saved to {pdf_path}")


def plot_single_ms2(ax, mz_arr, int_arr, name, adduct, precursor_mz, rt, scan_idx, explained_intensity, status,
                    discard_reason=None):
    """
    Plot single MS2 spectrum on a given axes
    """
    ax.vlines(mz_arr, 0, int_arr, color='0.4', linewidth=1)
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')

    title = f'{name} - {adduct}\n'
    title += f'm/z {precursor_mz:.4f}, RT: {rt:.2f} min, Scan {scan_idx}\n'
    title += f'Explained: {explained_intensity:.2f}, Status: {status}'

    if discard_reason:
        title += f'\nDiscard reason: {discard_reason}'

    ax.set_title(title, fontsize=8)
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.set_xlim(right=precursor_mz * 1.025)


def plot_mz_rt(df, mzml_name, out_dir):
    """
    Plot m/z vs RT for all compounds in the dataframe, handling multiple labels per spectrum
    """
    df = df[~pd.isnull(df['best_MS2_scan_idx'])].reset_index(drop=True)

    # Function to determine status, prioritizing 'Selected'
    def get_status(reasons):
        if pd.isnull(reasons) or reasons == '':
            return 'Selected'
        reasons_list = reasons.split(';')
        return 'Selected' if 'Selected' in reasons_list else reasons_list[0]

    # Apply the function to create the 'status' column
    df['status'] = df['discard_reason'].apply(get_status)

    # Sort by status, first 'Selected', then others
    df['status_sort'] = df['status'].map({'Selected': 0}).fillna(1)
    df = df.sort_values(['status_sort', 'status']).drop('status_sort', axis=1)

    # Dereplicate the best_MS2_scan_idx column
    df = df.drop_duplicates(subset='best_MS2_scan_idx', keep='first')

    # Define a set of distinct colors
    distinct_colors = {
        'Selected': '#1f77b4',
        'No core adduct': '#d62728',
        'MS2 explanation below cutoff': '#2ca02c',
        'Matched to doubly charged ion': '#ff7f0e',
        'Other': '#9467bd'  # For any other status that might appear
    }

    # Create the plot
    plt.figure(figsize=(12, 8))

    for status in df['status'].unique():
        subset = df[df['status'] == status]
        color = distinct_colors.get(status, distinct_colors['Other'])
        plt.scatter(subset['RT'], subset['m/z'], c=[color], label=status, alpha=0.7)

    plt.xlabel('RT (min)')
    plt.ylabel('m/z')
    plt.title('m/z vs Retention Time')

    # Adjust the plot area to make room for the legend
    plt.subplots_adjust(right=0.75)

    # Place the legend outside the plot area
    plt.legend(title='Status', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Save the plot
    output_path = f"{out_dir}/{mzml_name}_mz_rt_plot.svg"
    plt.savefig(output_path, transparent=True, bbox_inches='tight')
    plt.close()

    print(f"m/z vs RT plot saved to {output_path}")

    # Print some statistics for verification
    print(f"Total number of data points: {len(df)}")
    print(df['status'].value_counts())
