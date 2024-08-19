import pandas as pd

from .raw_data_utils import MSData
from .main import init_config
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_all_eic(df, file_path, out_dir, plots_per_page=5):
    """
    Plot all EIC on a PDF with an M x 1 layout per page

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

    df = df[(df['selected']) & (~pd.isnull(df['best_MS2_scan_idx']))].reset_index(drop=True)

    pdf_path = f"{out_dir}/{file_name}_eic_plots.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, axs = plt.subplots(plots_per_page, 1, figsize=(8.27, 11.69))  # A4 size in inches
        fig.suptitle(f"EIC Plots - {file_name}", fontsize=14)
        plot_count = 0

        for compound_name, group in df.groupby('compound_name'):
            unique_adducts = group[['t_adduct', 't_mz']].drop_duplicates()

            ax = axs[plot_count] if plots_per_page > 1 else axs
            ax.set_title(f"{compound_name}", fontsize=10)
            ax.set_xlabel('RT (min)', fontsize=8)
            ax.set_ylabel('Intensity', fontsize=8)

            # Dictionary to store EIC colors
            eic_colors = {}

            for _, row in unique_adducts.iterrows():
                t_mz = row['t_mz']
                t_adduct = row['t_adduct']

                eic_rt, eic_int, _, _ = d.get_eic_data(t_mz, mz_tol=0.005)
                line, = ax.plot(eic_rt, eic_int, label=f"{t_adduct} (m/z {t_mz:.4f})", linewidth=0.8, alpha=0.7)
                eic_colors[t_adduct] = line.get_color()  # Store the color for this adduct

            # Get the y-axis limits for positioning the MS2 labels
            y_min, y_max = ax.get_ylim()
            y_range = y_max - y_min

            for _, row in group.iterrows():
                if not pd.isnull(row['best_MS2_scan_idx']):
                    ms2_rt = d.scans[int(row['best_MS2_scan_idx']) - 1].rt
                    # Use the same color as the EIC for this adduct
                    ax.axvline(x=ms2_rt, color=eic_colors[row['t_adduct']], linestyle='--', linewidth=0.8)

                    # # Position the label within the plot, slightly above the x-axis
                    # label_y = y_max - 0.05 * y_range
                    # ax.text(ms2_rt, label_y, row['t_adduct'],
                    #         rotation=0, va='bottom', ha='center', fontsize=8,
                    #         bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0),
                    #         color=eic_colors[row['t_adduct']])  # Use the same color for the text

            ax.legend(fontsize=8, loc='upper right')
            ax.tick_params(axis='both', which='major', labelsize=8)

            plot_count += 1

            if plot_count == plots_per_page:
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig(fig)
                plt.close(fig)
                fig, axs = plt.subplots(plots_per_page, 1, figsize=(8.27, 11.69))
                fig.suptitle(f"EIC Plots - {file_name}", fontsize=14)
                plot_count = 0

        # Save the last page if there are any remaining plots
        if plot_count > 0:
            for i in range(plot_count, plots_per_page):
                if plots_per_page > 1:
                    axs[i].axis('off')
                else:
                    axs.axis('off')
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"EIC plots saved to {pdf_path}")


def plot_all_ms2(df, file_path, out_dir):
    """
    Plot all MS2 spectra on a PDF with a 5x2 layout per page
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
            if row['selected'] and not pd.isnull(row['best_MS2_scan_idx']):
                this_scan = d.scans[int(row['best_MS2_scan_idx']) - 1]
                mz_arr, int_arr = this_scan.peaks[:, 0], this_scan.peaks[:, 1]

                ax = axs[plot_count // 2, plot_count % 2]
                plot_single_ms2(ax, mz_arr, int_arr, row['compound_name'], row['t_adduct'],
                                this_scan.precursor_mz, row['RT'], row['best_MS2_scan_idx'],
                                row['ms2_explained_intensity'])

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


def plot_single_ms2(ax, mz_arr, int_arr, name, adduct, precursor_mz, rt, scan_idx, explained_intensity):
    """
    Plot single MS2 spectrum on a given axes
    """
    ax.vlines(mz_arr, 0, int_arr, color='0.8', linewidth=1)
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_title(f'{name} {adduct}\nm/z {precursor_mz:.4f}, RT: {rt:.2f} min\nScan {scan_idx}, Explained {explained_intensity:.2f}',
                 fontsize=8)
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.set_xlim(right=precursor_mz * 1.025)
