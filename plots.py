import os
import re
import shutil
import filecmp
import json
import random
import itertools
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from datetime import datetime
import MDAnalysis as mda
from scipy.stats import gaussian_kde
import glob
from joblib import Parallel, delayed

# ============================================================================
# CENTRALIZED COLOR MAPPINGS
# ============================================================================
# These color mappings are used consistently across all plotting functions
# to ensure visual consistency across the entire codebase.

# Use matplotlib's built-in colorblind-friendly style across all figures.
try:
    plt.style.use('tableau-colorblind10')
except OSError:
    # Fallback to default style if the tableau style is unavailable.
    pass

_mpl_cycle = plt.rcParams.get('axes.prop_cycle')
_mpl_colors = _mpl_cycle.by_key().get('color', []) if _mpl_cycle is not None else []
if len(_mpl_colors) < 7:
    _mpl_colors = [
        '#006BA4', '#FF800E', '#ABABAB', '#595959', '#5F9ED1',
        '#C85200', '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF'
    ]

# Colorblind-friendly palettes for consistent categorical plotting.
COLORBLIND_PALETTE = _mpl_colors[:7]
EXTENDED_COLORBLIND_PALETTE = list(dict.fromkeys(
    _mpl_colors + sns.color_palette('colorblind', n_colors=10).as_hex() + ['#7F7F7F', '#BCBD22', '#17BECF']
))
sns.set_palette(COLORBLIND_PALETTE)

# Standard axial ligand combination color mapping
AXIAL_LIGAND_COLOR_MAP = {
    'CYS-HOH': COLORBLIND_PALETTE[0],
    'HIS-HIS': COLORBLIND_PALETTE[1],
    'HIS-HOH': COLORBLIND_PALETTE[2],
    'HIS-MET': COLORBLIND_PALETTE[3],
    'HIS-OXY': COLORBLIND_PALETTE[4]
}

# Standard charge-multiplicity combination color mapping
CHARGE_MULTIPLICITY_COLOR_MAP = {
    '0,1': COLORBLIND_PALETTE[0],
    '0_1': COLORBLIND_PALETTE[0],
    '0-1': COLORBLIND_PALETTE[0],
    '0,2': COLORBLIND_PALETTE[1],
    '0_2': COLORBLIND_PALETTE[1],
    '0-2': COLORBLIND_PALETTE[1],
    '0,3': COLORBLIND_PALETTE[2],
    '0_3': COLORBLIND_PALETTE[2],
    '0-3': COLORBLIND_PALETTE[2],
    '0,5': COLORBLIND_PALETTE[1],
    '0_5': COLORBLIND_PALETTE[1],
    '0-5': COLORBLIND_PALETTE[1],
    '1,2': COLORBLIND_PALETTE[2],
    '1_2': COLORBLIND_PALETTE[2],
    '1-2': COLORBLIND_PALETTE[2],
    '1,6': COLORBLIND_PALETTE[3],
    '1_6': COLORBLIND_PALETTE[3],
    '1-6': COLORBLIND_PALETTE[3],
    '1,5': COLORBLIND_PALETTE[4],
    '1_5': COLORBLIND_PALETTE[4],
    '1-5': COLORBLIND_PALETTE[4],
    '-1,2': COLORBLIND_PALETTE[5],
    '-1_2': COLORBLIND_PALETTE[5],
    '-1-2': COLORBLIND_PALETTE[5],
    '-1,4': COLORBLIND_PALETTE[6],
    '-1_4': COLORBLIND_PALETTE[6],
    '-1-4': COLORBLIND_PALETTE[6],
    '-1,6': COLORBLIND_PALETTE[0],
    '-1_6': COLORBLIND_PALETTE[0],
    '-1-6': COLORBLIND_PALETTE[0]
}

def _format_charge_multiplicity_number(value):
    """Normalize numeric charge/multiplicity tokens to compact strings."""
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return str(value).strip()
    if numeric_value.is_integer():
        return str(int(numeric_value))
    return str(numeric_value)

def canonicalize_charge_multiplicity_combo(combo):
    """Return a canonical 'charge,multiplicity' label for common combo formats."""
    try:
        if pd.isna(combo):
            return None
    except TypeError:
        pass

    if isinstance(combo, (tuple, list)) and len(combo) == 2:
        charge, multiplicity = combo
        return (
            f"{_format_charge_multiplicity_number(charge)},"
            f"{_format_charge_multiplicity_number(multiplicity)}"
        )

    combo_str = str(combo).strip()
    match = re.match(
        r'^\s*(?:q\s*=\s*)?(-?\d+(?:\.\d+)?)\s*(?:,|_|-|/)\s*(?:m\s*=\s*)?(-?\d+(?:\.\d+)?)\s*$',
        combo_str
    )
    if not match:
        return combo_str

    return (
        f"{_format_charge_multiplicity_number(match.group(1))},"
        f"{_format_charge_multiplicity_number(match.group(2))}"
    )

def charge_multiplicity_sort_key(value):
    """Sort charge-multiplicity labels consistently across supported formats."""
    canonical = canonicalize_charge_multiplicity_combo(value)
    if canonical is None:
        return (float('inf'), float('inf'), str(value))

    try:
        charge_str, multiplicity_str = canonical.split(',', 1)
        return (float(charge_str), float(multiplicity_str), str(value))
    except (TypeError, ValueError):
        return (float('inf'), float('inf'), str(value))

# Helper function to get colors for axial ligand combinations
def get_axial_ligand_colors(combinations, extended=False):
    """
    Get consistent colors for axial ligand combinations.
    
    Parameters:
    -----------
    combinations : list
        List of axial ligand combination strings
    extended : bool, default=False
        Whether to use extended palette for many combinations
    
    Returns:
    --------
    dict : Dictionary mapping combinations to colors
    """
    color_map = {}
    palette = EXTENDED_COLORBLIND_PALETTE if extended else COLORBLIND_PALETTE
    
    for i, combo in enumerate(combinations):
        if combo in AXIAL_LIGAND_COLOR_MAP:
            color_map[combo] = AXIAL_LIGAND_COLOR_MAP[combo]
        else:
            color_map[combo] = palette[i % len(palette)]
    
    return color_map

# Helper function to get colors for charge-multiplicity combinations
def get_charge_multiplicity_colors(combinations, extended=False):
    """
    Get consistent colors for charge-multiplicity combinations.
    
    Parameters:
    -----------
    combinations : list
        List of charge-multiplicity combination strings
    extended : bool, default=False
        Whether to use extended palette for many combinations
    
    Returns:
    --------
    dict : Dictionary mapping combinations to colors
    """
    color_map = {}
    palette = EXTENDED_COLORBLIND_PALETTE if extended else COLORBLIND_PALETTE
    
    for i, combo in enumerate(combinations):
        canonical_combo = canonicalize_charge_multiplicity_combo(combo)
        if canonical_combo in CHARGE_MULTIPLICITY_COLOR_MAP:
            color_map[combo] = CHARGE_MULTIPLICITY_COLOR_MAP[canonical_combo]
        elif combo in CHARGE_MULTIPLICITY_COLOR_MAP:
            color_map[combo] = CHARGE_MULTIPLICITY_COLOR_MAP[combo]
        else:
            color_map[combo] = palette[i % len(palette)]
    
    return color_map

def get_color_sequence(n, offset=0):
    """Return n deterministic colors from the shared colorblind palette."""
    if n <= 0:
        return []
    return [EXTENDED_COLORBLIND_PALETTE[(offset + i) % len(EXTENDED_COLORBLIND_PALETTE)] for i in range(n)]

class plots_class:
    """
    A class for plotting HOMO/LUMO values (the last 10 alpha_occ and first 10 alpha_virt)
    from multiple Gaussian-log-derived JSON files, with support for grouping by presence in axial positions.
    """

    def __init__(self):
        self.json_dir = "/home/pbuser/Desktop/PhD_WORK/heme/jsons/"
        self.split_legend_flag = False
        self.legend_split_variables = ["axials"]

    def add_concatenated_column(self, df, col1, col2, new_col_name, separator=" "):
        df[new_col_name] = df[col1].astype(str) + separator + df[col2].astype(str)
        return df

    def create_charge_multiplicity_column(self, df):
        """
        Create a combined charge_multiplicity column from separate charge and multiplicity columns.
        Format: f'{charge}_{multiplicity}'
        """
        if 'charge' in df.columns and 'multiplicity' in df.columns:
            df['charge_multiplicity'] = df['charge'].astype(str) + '_' + df['multiplicity'].astype(str)
        return df

    def preprocess_dataframe(self, df, decode_dicts):
        """
        Decode only the relevant columns in decode_dicts that actually exist in df,
        then concatenate axial1 and axial2 into 'axials' if present.
        """
        # If already decoded once, skip
        if df.attrs.get("decoded", False):
            return df
        # Decode mappings only for columns present in df
        if decode_dicts and isinstance(decode_dicts, dict):
            for col, mapping in decode_dicts.items():
                if col not in df.columns:
                    continue
                # Map values directly (replacements are already in {code: string} format)
                # No need to invert - mapping is already code -> string
                df[col] = df[col].map(lambda x: mapping.get(x, x))
        # Only create 'axials' if both source columns exist
        if 'axial1' in df.columns and 'axial2' in df.columns:
            self.add_concatenated_column(
                df,
                col1="axial1",
                col2="axial2",
                new_col_name="axials",
                separator="-"
            )
        # Mark as decoded
        df.attrs["decoded"] = True
        return df

    def filter_negative_one(self, df, col_x, col_y):
        return df[(df[col_x] != -1) & (df[col_y] != -1)]

    def plot_categorical_categorical_heatmap(
        self,
        df,
        col_x,
        col_y,
        decode_dicts=None,
        ax=None,
        save=False,
        save_dir='plots',
        filename=None,
        cmap='cividis',
        annot=True,
        fmt='d',
        row_order=None,
        col_order=None,
        include_all_levels=False,
        cell_size=0.8,
        min_figsize=(4.0, 3.0),
        filter_col=None,
        filter_value=True,
        drop_negative_one=True
    ):
        """
        Plot a categorical-vs-categorical count heatmap (crosstab), adapted from
        scripts/old_plots.py.

        Parameters
        ----------
        df : pandas.DataFrame
            Input dataframe.
        col_x, col_y : str
            Categorical columns to cross-tabulate.
        decode_dicts : dict, optional
            Optional decode mapping used by preprocess_dataframe.
        ax : matplotlib.axes.Axes, optional
            Existing axis to draw on. If None, creates a new figure.
        save : bool, default=False
            If True and this method created its own figure, save plot to disk.
        save_dir : str, default='plots'
            Directory used when save=True.
        filename : str, optional
            Output filename. Defaults to heatmap_{col_x}_{col_y}.png.
        cmap : str, default='cividis'
            Matplotlib colormap for heatmap fill.
        annot : bool, default=True
            Show numeric counts on cells.
        fmt : str, default='d'
            Number format for annotations.
        row_order, col_order : list, optional
            Explicit order (and full set) of heatmap rows/columns.
        include_all_levels : bool, default=False
            If True, include all category levels (showing zero-count cells too).
            When decode_dicts has mappings for col_x/col_y, mapped levels are used.
        cell_size : float, default=0.8
            Approximate size (inches) of each heatmap cell when creating own figure.
            Lower values make each field smaller.
        min_figsize : tuple, default=(4.0, 3.0)
            Minimum figure size used with cell_size auto-sizing.
        filter_col : str, optional
            If set, filter dataframe to rows where filter_col == filter_value.
            Useful for converged-only heatmaps.
        filter_value : any, default=True
            Value used with filter_col.
        drop_negative_one : bool, default=True
            If True, drop rows where col_x or col_y equals -1.

        Returns
        -------
        tuple
            (crosstab_dataframe, axis)
        """
        if col_x not in df.columns or col_y not in df.columns:
            raise ValueError(f"Columns '{col_x}' and/or '{col_y}' not found in dataframe.")

        if filter_col is not None and filter_col not in df.columns:
            raise ValueError(f"Filter column '{filter_col}' not found in dataframe.")

        os.makedirs(save_dir, exist_ok=True)

        df_plot = df.copy()
        if decode_dicts is not None and not df_plot.attrs.get("decoded", False):
            df_plot = self.preprocess_dataframe(df_plot, decode_dicts)

        # Keep a copy before filtering so include_all_levels can retain full grid.
        df_levels = df_plot.copy()

        if filter_col is not None:
            df_plot = df_plot[df_plot[filter_col] == filter_value]

        if drop_negative_one:
            df_plot = self.filter_negative_one(df_plot, col_x, col_y)

        df_plot = df_plot.dropna(subset=[col_x, col_y])

        # Build full level lists when requested, so zero-count cells still appear.
        def _decoded_levels(col_name, explicit_order):
            if explicit_order is not None:
                return list(explicit_order)
            if not include_all_levels:
                return None
            if decode_dicts is not None and col_name in decode_dicts:
                mapping = decode_dicts[col_name]
                if isinstance(mapping, dict):
                    return [mapping[k] for k in sorted(mapping.keys())]
            values = df_levels[col_name].dropna()
            if drop_negative_one:
                values = values[values != -1]
            return sorted(values.astype(str).unique().tolist())

        levels_x = _decoded_levels(col_x, row_order)
        levels_y = _decoded_levels(col_y, col_order)

        x_data = df_plot[col_x]
        y_data = df_plot[col_y]
        if levels_x is not None:
            x_data = pd.Categorical(x_data, categories=levels_x, ordered=True)
        if levels_y is not None:
            y_data = pd.Categorical(y_data, categories=levels_y, ordered=True)

        ct = pd.crosstab(x_data, y_data, dropna=False)

        if ct.empty:
            print(f"Warning: Crosstab for {col_x} vs {col_y} is empty. Skipping plot.")
            return ct, ax

        own_fig = False
        if ax is None:
            own_fig = True
            fig_w = max(min_figsize[0], ct.shape[1] * cell_size)
            fig_h = max(min_figsize[1], ct.shape[0] * cell_size)
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        sns.heatmap(ct, annot=annot, fmt=fmt, cmap=cmap, ax=ax, cbar=True)
        ax.set_xlabel(col_y, fontsize=12)
        ax.set_ylabel(col_x, fontsize=12)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        if own_fig:
            plt.tight_layout()
            if save:
                out_name = filename or f"heatmap_{col_x}_{col_y}.png"
                fig.savefig(os.path.join(save_dir, out_name), dpi=600, bbox_inches='tight')
                print(f"Heatmap saved as {out_name}")
                plt.close(fig)
            else:
                plt.show()

        return ct, ax

    def plot_convergence_heatmap(
        self,
        df,
        col_x='axial1',
        col_y='axial2',
        converged_col='converged',
        converged_value=True,
        decode_dicts=None,
        ax=None,
        save=False,
        save_dir='plots',
        filename=None,
        cmap='cividis',
        annot=True,
        fmt='d',
        row_order=None,
        col_order=None,
        include_all_levels=True,
        cell_size=0.8,
        min_figsize=(4.0, 3.0),
        drop_negative_one=True
    ):
        """
        Convenience wrapper for a converged-structures heatmap.
        """
        out_name = filename or f"heatmap_converged_{col_x}_{col_y}.png"
        return self.plot_categorical_categorical_heatmap(
            df=df,
            col_x=col_x,
            col_y=col_y,
            decode_dicts=decode_dicts,
            ax=ax,
            save=save,
            save_dir=save_dir,
            filename=out_name,
            cmap=cmap,
            annot=annot,
            fmt=fmt,
            row_order=row_order,
            col_order=col_order,
            include_all_levels=include_all_levels,
            cell_size=cell_size,
            min_figsize=min_figsize,
            filter_col=converged_col,
            filter_value=converged_value,
            drop_negative_one=drop_negative_one
        )

    def plot_x_random(self, df, split=False, charge_color=False, mult_color=False, axial_color=False, decode_dicts=None):
        """
        Line-plots the HOMO/LUMO columns in a single figure.

        Supports both the original 'homo[i]'/'lumo[i]' naming **and**
        the renamed 'HOMO'/'HOMO-#' & 'LUMO+%' conventions from DataPreprocessor.
        
        Parameters:
        -----------
        df : pandas.DataFrame
            DataFrame containing HOMO/LUMO energy data
        split : bool, default=False
            If True, sample 50% of data for plotting
        charge_color : bool, default=False
            If True, color lines by charge values
        mult_color : bool, default=False
            If True, color lines by multiplicity values
        axial_color : bool, default=False
            If True, color lines by axial ligand combinations (axial1-axial2)
        decode_dicts : dict, default=None
            Dictionary for decoding categorical variables from numeric codes to strings.
            Expected format is the encoding_dict returned by preprocessor.process().
            When axial_color=True, this enables meaningful axial ligand labels.
        """
        # 1) subsample if needed
        if split:
            df_plot = df.sample(frac=0.5, random_state=42)
        else:
            df_plot = df.copy()

        # 1.5) decode categorical variables if needed
        if decode_dicts is not None and not df_plot.attrs.get("decoded", False):
            df_plot = self.preprocess_dataframe(df_plot, decode_dicts)

        # 2) discover HOMO/LUMO column names & order them
        if any(col.startswith('homo[') for col in df_plot.columns):
            homo_cols = [f'homo[{i}]' for i in range(10) if f'homo[{i}]' in df_plot.columns]
            lumo_cols = [f'lumo[{i}]' for i in range(10) if f'lumo[{i}]' in df_plot.columns]
        else:
            homo_pat = re.compile(r'^HOMO(?:-(\d+))?$')
            lumo_pat = re.compile(r'^LUMO(?:\+(\d+))?$')
            pairs = []
            for c in df_plot.columns:
                m = homo_pat.match(c)
                if m:
                    off = int(m.group(1)) if m.group(1) else 0
                    idx = 9 - off
                    pairs.append((idx, c))
            homo_cols = [c for _, c in sorted(pairs)]
            pairs = []
            for c in df_plot.columns:
                m = lumo_pat.match(c)
                if m:
                    idx = int(m.group(1)) if m.group(1) else 0
                    pairs.append((idx, c))
            lumo_cols = [c for _, c in sorted(pairs)]
        if not homo_cols or not lumo_cols:
            raise ValueError(f"Could not find any HOMO/LUMO columns; looked for {homo_cols=} and {lumo_cols=}")
        plot_cols = homo_cols + lumo_cols

        # 3) decide on hue grouping and label formatting
        hue = None
        legend_title = None
        
        # Determine grouping priority: axial > charge-multiplicity > none
        if axial_color:
            if 'axial1' not in df_plot.columns or 'axial2' not in df_plot.columns:
                raise ValueError("df must contain 'axial1' and 'axial2' columns for axial coloring.")
            
            # Create axials column if it doesn't exist
            if 'axials' not in df_plot.columns:
                df_plot = self.add_concatenated_column(
                    df_plot, 
                    col1="axial1", 
                    col2="axial2", 
                    new_col_name="axials", 
                    separator="-"
                )
            
            # Handle combined axial + charge-multiplicity grouping
            if charge_color or mult_color:
                if 'charge' not in df_plot.columns or 'multiplicity' not in df_plot.columns:
                    raise ValueError("df must contain 'charge' and 'multiplicity' columns for combined coloring.")
                df_plot = self.create_charge_multiplicity_column(df_plot)
                df_plot["_group"] = df_plot["axials"] + " (" + df_plot["charge_multiplicity"] + ")"
                legend_title = "Axial + Charge–Multiplicity"
                
                # Define ordered axial ligands and their color families
                axial_order = ['HIS-HOH', 'HIS-HIS', 'HIS-OXY', 'HIS-MET', 'CYS-HOH']
                axial_color_families = {
                    'HIS-HOH': 'Blues',      # Blue family
                    'HIS-HIS': 'Greens',     # Green family
                    'HIS-OXY': 'Grays',    # Yellow family
                    'HIS-MET': 'Purples',    # Purple family
                    'CYS-HOH': 'Reds'        # Red family
                }
                
                # Sort groups by axial order, then by charge-multiplicity
                unique_groups = df_plot["_group"].unique()
                sorted_groups = []
                for axial in axial_order:
                    axial_groups = [g for g in unique_groups if g.startswith(axial)]
                    sorted_groups.extend(sorted(axial_groups))
                
                # Create custom color palette using color families
                custom_colors = []
                for group in sorted_groups:
                    axial = group.split(' (')[0]  # Extract axial part
                    if axial in axial_color_families:
                        color_family = axial_color_families[axial]
                        # Get colors from the color family
                        axial_groups_count = len([g for g in sorted_groups if g.startswith(axial)])
                        # Generate shades from the color family
                        cmap = plt.cm.get_cmap(color_family)
                        # Start from 0.4 to avoid too light colors, go to 0.9 for good contrast
                        shade_idx = sorted([g for g in sorted_groups if g.startswith(axial)]).index(group)
                        color_intensity = 0.4 + (0.5 * shade_idx / max(1, axial_groups_count - 1))
                        custom_colors.append(cmap(color_intensity))
                    else:
                        custom_colors.append('gray')  # Fallback color

                # Set the color palette
                sns.set_palette(custom_colors)
                
                # Reorder dataframe groups to match our desired order
                df_plot["_group"] = pd.Categorical(df_plot["_group"], categories=sorted_groups, ordered=True)
                
            else:
                df_plot["_group"] = df_plot["axials"]
                legend_title = "Axial Ligands"
                
                # Define ordered axial ligands for simple axial coloring
                axial_order = ['HIS-HOH', 'HIS-HIS', 'HIS-OXY', 'HIS-MET', 'CYS-HOH']
                # Filter to only existing axials in data
                existing_axials = [ax for ax in axial_order if ax in df_plot["axials"].values]
                df_plot["_group"] = pd.Categorical(df_plot["_group"], categories=existing_axials, ordered=True)
                
            hue = "_group"
            
        elif charge_color or mult_color:
            if 'charge' not in df_plot.columns or 'multiplicity' not in df_plot.columns:
                raise ValueError("df must contain 'charge' and 'multiplicity' columns for coloring.")
            df_plot = self.create_charge_multiplicity_column(df_plot)
            df_plot["_group"] = df_plot["charge_multiplicity"]
            legend_title = "Charge–Multiplicity"
            hue = "_group"

        # 4) melt to long form
        df_melt = (
            df_plot
            .reset_index(drop=False)
            .rename(columns={"index": "_row_id"})
            .melt(
                id_vars=["_row_id"] + ([hue] if hue else []),
                value_vars=plot_cols,
                var_name="orbital",
                value_name="energy",
            )
        )

        # 5) draw
        plt.figure(figsize=(12, 6))
        if hue:
            hue_order = None
            palette = None
            group_series = df_plot[hue].dropna()
            if hasattr(group_series.dtype, "categories"):
                hue_order = [group for group in group_series.dtype.categories if group in set(group_series.astype(str))]
            else:
                hue_order = group_series.astype(str).unique().tolist()

            if axial_color and not (charge_color or mult_color):
                palette = get_axial_ligand_colors(hue_order)
            elif (charge_color or mult_color) and not axial_color:
                hue_order = sorted(hue_order, key=charge_multiplicity_sort_key)
                palette = get_charge_multiplicity_colors(hue_order)

            ax = sns.lineplot(
                data=df_melt,
                x="orbital", y="energy",
                hue=hue,
                hue_order=hue_order,
                palette=palette,
                units="_row_id",
                estimator=None,
                alpha=0.8, linewidth=1.5
            )
            # manually rebuild legend with appropriate title and place it outside the plot area
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles, labels,
                title=legend_title if legend_title else "Charge–Multiplicity",
                fontsize=14,
                title_fontsize=16,
                bbox_to_anchor=(1.05, 1),
                loc='upper left'
            )
        else:
            sns.lineplot(
                data=df_melt,
                x="orbital", y="energy",
                units="_row_id", estimator=None,
                color="gray", alpha=0.8, linewidth=1,
                legend=False
            )

        # 6) highlight the frontier orbitals
        if 'HOMO' in plot_cols:
            h0 = plot_cols.index('HOMO')
        else:
            h0 = plot_cols.index('homo[9]')
        if 'LUMO' in plot_cols:
            l0 = plot_cols.index('LUMO')
        else:
            l0 = plot_cols.index('lumo[0]')
        plt.axvspan(h0 - 0.5, h0 + 0.5, color="lightblue", alpha=0.3)
        plt.axvspan(l0 - 0.5, l0 + 0.5, color="lightblue", alpha=0.3)

        # 7) tidy up
        plt.xticks(ticks=np.arange(len(plot_cols)), labels=plot_cols, rotation=45, ha="right", fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("", fontsize=14)
        plt.ylabel("Energy [a.u.]", fontsize=16)
        #title = "HOMO/LUMO energy profiles"
        #if split:
        #    title += " (50% sample)"
        #plt.title(title, fontsize=18)
        plt.tight_layout()
        plt.show()

    def plot_grouped_histogram_with_contours(self, df, column, group_by='axial', bins=10, 
                                           save=False, save_dir='.', opacity=0.6, 
                                           color_palette=None, figsize=(12, 8),
                                           contour_linewidth=2, kde_bandwidth=None):
        """
        Plot histogram of a single column grouped by either axial ligand combinations 
        or charge-multiplicity combinations, with colored contour lines outlining each group.
        
        Parameters
        ----------
        df : pandas.DataFrame
            The dataframe containing the data
        column : str
            The column name to plot histogram for
        group_by : str, default 'axial'
            Grouping method: 'axial' for axial1-axial2 combinations, 
            'charge' for charge-multiplicity combinations
        bins : int, default 10
            Number of histogram bins
        save : bool, default False
            Whether to save the plot
        save_dir : str, default '.'
            Directory to save the plot
        opacity : float, default 0.6
            Transparency of histogram bars
        color_palette : list, optional
            Custom color palette for groups
        figsize : tuple, default (12, 8)
            Figure size (width, height) in inches
        contour_linewidth : float, default 2
            Width of the contour lines
        kde_bandwidth : float, optional
            Bandwidth for KDE estimation. If None, uses automatic selection
        """
        # Decoding dictionary for encoded categorical variables
        REPLACEMENTS = {
            'axial1': {0: 'CYS', 1: 'HIS'},
            'axial2': {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'},
            'function': {
                0: 'electron transporter',
                1: 'hydrolase',
                2: 'oxidoreductase',
                3: 'oxygen carrier',
                4: 'transferase',
                5: 'unclassified'
            }
        }
        
        # Color palette for histograms (colorblind-friendly)
        COLOR_PALETTE = COLORBLIND_PALETTE
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        # Check if column exists
        if column not in df.columns:
            print(f"Column '{column}' not found in dataframe")
            return
        
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2', 'function'])
        
        # Filter out missing values (-1 and NaN)
        df_plot = df_dec[(df_dec[column] != -1) & df_dec[column].notna()].copy()
        
        if df_plot.empty:
            print(f"No valid data found for column '{column}'")
            return
        
        # Set up grouping variables
        if group_by == 'axial':
            df_plot['group_combo'] = df_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
            group_title = 'Axial Ligand Combinations'
        elif group_by == 'charge':
            df_plot['group_combo'] = df_plot.apply(lambda r: f"q={r.charge},m={r.multiplicity}", axis=1)
            group_title = 'Charge-Multiplicity Combinations'
        else:
            raise ValueError("group_by must be either 'axial' or 'charge'")
        
        # Get unique combinations
        if group_by == 'axial':
            unique_groups = sorted(df_plot['group_combo'].unique())
            color_map = (
                {group: color_palette[i % len(color_palette)] for i, group in enumerate(unique_groups)}
                if color_palette else get_axial_ligand_colors(unique_groups)
            )
        else:
            unique_groups = sorted(df_plot['group_combo'].unique(), key=charge_multiplicity_sort_key)
            color_map = (
                {group: color_palette[i % len(color_palette)] for i, group in enumerate(unique_groups)}
                if color_palette else get_charge_multiplicity_colors(unique_groups)
            )
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Store histogram data for contour calculation
        histogram_data = {}
        
        # Create histogram for each group
        for group in unique_groups:
            group_data = df_plot[df_plot['group_combo'] == group][column]
            
            if not group_data.empty:
                color = color_map[group]
                
                # Create histogram
                n, bins_edges, patches = ax.hist(group_data, bins=bins, alpha=opacity, 
                                                color=color, label=group, 
                                                edgecolor='black', linewidth=0.5)
                
                # Store data for contour calculation
                histogram_data[group] = {
                    'counts': n,
                    'bin_edges': bins_edges,
                    'color': color,
                    'data': group_data
                }
        
        # Add contour lines for each group
        for group, hist_info in histogram_data.items():
            color = hist_info['color']
            counts = hist_info['counts']
            bin_edges = hist_info['bin_edges']
            
            # Calculate bin centers
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            # Create a smooth line through the histogram peaks
            # Use interpolation to create smooth contour
            from scipy.interpolate import interp1d
            from scipy.ndimage import gaussian_filter1d
            
            # Smooth the histogram counts
            smoothed_counts = gaussian_filter1d(counts, sigma=0.8)
            
            # Create more points for smoother line
            if len(bin_centers) > 3:  # Need at least 4 points for interpolation
                try:
                    # Create interpolation function
                    f = interp1d(bin_centers, smoothed_counts, kind='cubic', 
                               fill_value=0, bounds_error=False)
                    
                    # Create finer x points for smooth line
                    x_smooth = np.linspace(bin_centers.min(), bin_centers.max(), 200)
                    y_smooth = f(x_smooth)
                    
                    # Plot the contour line
                    ax.plot(x_smooth, y_smooth, color=color, linewidth=contour_linewidth,
                           alpha=0.9, linestyle='-')
                           
                except Exception as e:
                    # Fallback to simple line plot if interpolation fails
                    ax.plot(bin_centers, smoothed_counts, color=color, 
                           linewidth=contour_linewidth, alpha=0.9, linestyle='-')
            else:
                # For few points, just connect them directly
                ax.plot(bin_centers, smoothed_counts, color=color, 
                       linewidth=contour_linewidth, alpha=0.9, linestyle='-')
        
        # Set labels and title
        if column == 'homo_lumo_gap':
            ax.set_xlabel('LUMO-HOMO [eV]', fontsize=20)
            
            # Create secondary x-axis with nm conversion and custom tick spacing
            ax2 = ax.secondary_xaxis('top')
            
            # Get the eV range and create exactly 10 round nm values
            ev_min, ev_max = ax.get_xlim()
            nm_max = 1239.84 / ev_min if ev_min > 0 else 2000  # Avoid division by zero
            nm_min = 1239.84 / ev_max if ev_max > 0 else 400
            
            # Create exactly 10 nice round numbers across the nm range
            # Find appropriate step size for 10 ticks
            nm_range = nm_max - nm_min
            base_step = nm_range / 9  # 9 intervals for 10 points
            
            # Round step to nice numbers (powers of 10, 2, 5)
            magnitude = 10 ** np.floor(np.log10(base_step))
            normalized_step = base_step / magnitude
            if normalized_step <= 1:
                nice_step = 1 * magnitude
            elif normalized_step <= 2:
                nice_step = 2 * magnitude
            elif normalized_step <= 5:
                nice_step = 5 * magnitude
            else:
                nice_step = 10 * magnitude
            
            # Generate exactly 10 ticks
            start_tick = np.ceil(nm_min / nice_step) * nice_step
            nm_ticks = [start_tick + i * nice_step for i in range(10)]
            
            # Filter to keep only ticks within reasonable range
            nm_ticks = [nm for nm in nm_ticks if nm_min <= nm <= nm_max * 1.1][:10]
            
            # Convert nm ticks back to eV positions
            ev_positions = [1239.84 / nm for nm in nm_ticks if nm > 0]
            
            # Set custom ticks and labels with scientific notation for large numbers
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(ev_positions)
            
            # Format labels - use scientific notation if numbers are large
            if max(nm_ticks) > 10000:
                labels = [f'{nm/1000:.1f}k' for nm in nm_ticks]
            else:
                labels = [f'{nm:.0f}' for nm in nm_ticks]
            
            ax2.set_xticklabels(labels)
            ax2.set_xlabel('Wavelength [nm]', fontsize=18)
            
            # Rotate labels if they might overlap
            ax2.tick_params(axis='x', rotation=45, labelsize=14)
            
            # Add vertical grid lines for nm axis intervals
            for ev_pos in ev_positions:
                ax.axvline(x=ev_pos, color='lightblue', linestyle='--', alpha=0.7, linewidth=0.8)
            
        else:
            ax.set_xlabel(column, fontsize=20)
        
        ax.set_ylabel('Frequency', fontsize=20)
        
        # Add legend (without title)
        ax.legend(loc='upper right', bbox_to_anchor=(0.75, 1), 
                 frameon=False, fancybox=False, shadow=False, fontsize=16)
        
        # Set tick font sizes
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save if requested
        if save:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            filename = f"{column}_grouped_contours_{group_by}.png"
            fig.savefig(os.path.join(save_dir, filename), bbox_inches='tight')
            print(f"Plot saved as {filename}")
        
        plt.show()

    def plot_combined_homo_lumo_gap_histogram(self, df, column='homo_lumo_gap', bins=10, 
                                            save=False, save_dir='plots', opacity=0.6, 
                                            color_palette=None, figsize=(16, 6),
                                            contour_linewidth=2, kde_bandwidth=None):
        """
        Create combined plot with both axial and charge-multiplicity grouped HOMO-LUMO gap histograms.
        Saves as homo-lumo-gap_histogram.png
        """
        # Decoding dictionary for encoded categorical variables
        REPLACEMENTS = {
            'axial1': {0: 'CYS', 1: 'HIS'},
            'axial2': {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'},
        }
        
        # Color palette for histograms (colorblind-friendly)
        COLOR_PALETTE = COLORBLIND_PALETTE
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        # Check if column exists
        if column not in df.columns:
            print(f"Column '{column}' not found in dataframe")
            return
        
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2'])
        
        # Filter out missing values (-1 and NaN)
        df_plot = df_dec[(df_dec[column] != -1) & df_dec[column].notna()].copy()
        
        if df_plot.empty:
            print(f"No valid data found for column '{column}'")
            return
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        def create_subplot(ax, group_by_type, group_title):
            # Set up grouping variables
            if group_by_type == 'axial':
                df_plot['group_combo'] = df_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
            elif group_by_type == 'charge':
                df_plot['group_combo'] = df_plot.apply(lambda r: f"q={r.charge},m={r.multiplicity}", axis=1)
            
            # Get unique combinations
            if group_by_type == 'axial':
                unique_groups = sorted(df_plot['group_combo'].unique())
                color_map = (
                    {group: color_palette[i % len(color_palette)] for i, group in enumerate(unique_groups)}
                    if color_palette else get_axial_ligand_colors(unique_groups)
                )
            else:
                unique_groups = sorted(df_plot['group_combo'].unique(), key=charge_multiplicity_sort_key)
                color_map = (
                    {group: color_palette[i % len(color_palette)] for i, group in enumerate(unique_groups)}
                    if color_palette else get_charge_multiplicity_colors(unique_groups)
                )
            
            # Store histogram data for contour calculation
            histogram_data = {}
            
            # Create histogram for each group
            for group in unique_groups:
                group_data = df_plot[df_plot['group_combo'] == group][column]
                
                if not group_data.empty:
                    color = color_map[group]
                    
                    # Create histogram
                    n, bins_edges, patches = ax.hist(group_data, bins=bins, alpha=opacity, 
                                                    color=color, label=group, 
                                                    edgecolor='black', linewidth=0.5)
                    
                    # Store data for contour calculation
                    histogram_data[group] = {
                        'counts': n,
                        'bin_edges': bins_edges,
                        'color': color,
                        'data': group_data
                    }
            
            # Add contour lines for each group
            for group, hist_info in histogram_data.items():
                color = hist_info['color']
                counts = hist_info['counts']
                bin_edges = hist_info['bin_edges']
                
                # Calculate bin centers
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                
                # Create a smooth line through the histogram peaks
                from scipy.interpolate import interp1d
                from scipy.ndimage import gaussian_filter1d
                
                # Smooth the histogram counts
                smoothed_counts = gaussian_filter1d(counts, sigma=0.8)
                
                # Create more points for smoother line
                if len(bin_centers) > 3:  # Need at least 4 points for interpolation
                    try:
                        # Create interpolation function
                        f = interp1d(bin_centers, smoothed_counts, kind='cubic', 
                                   fill_value=0, bounds_error=False)
                        
                        # Create finer x points for smooth line
                        x_smooth = np.linspace(bin_centers.min(), bin_centers.max(), 200)
                        y_smooth = f(x_smooth)
                        
                        # Plot the contour line
                        ax.plot(x_smooth, y_smooth, color=color, linewidth=contour_linewidth,
                               alpha=0.9, linestyle='-')
                               
                    except Exception as e:
                        # Fallback to simple line plot if interpolation fails
                        ax.plot(bin_centers, smoothed_counts, color=color, 
                               linewidth=contour_linewidth, alpha=0.9, linestyle='-')
                else:
                    # For few points, just connect them directly
                    ax.plot(bin_centers, smoothed_counts, color=color, 
                           linewidth=contour_linewidth, alpha=0.9, linestyle='-')
            
            # Set labels
            ax.set_xlabel('LUMO-HOMO [eV]', fontsize=20)
            ax.set_ylabel('Frequency', fontsize=20)
            
            # Create secondary x-axis with nm conversion
            ax2_top = ax.secondary_xaxis('top')
            
            # Get the eV range and create exactly 10 round nm values
            ev_min, ev_max = ax.get_xlim()
            nm_max = 1239.84 / ev_min if ev_min > 0 else 2000
            nm_min = 1239.84 / ev_max if ev_max > 0 else 400
            
            # Create exactly 10 nice round numbers across the nm range
            nm_range = nm_max - nm_min
            base_step = nm_range / 9
            
            # Round step to nice numbers
            magnitude = 10 ** np.floor(np.log10(base_step))
            normalized_step = base_step / magnitude
            if normalized_step <= 1:
                nice_step = 1 * magnitude
            elif normalized_step <= 2:
                nice_step = 2 * magnitude
            elif normalized_step <= 5:
                nice_step = 5 * magnitude
            else:
                nice_step = 10 * magnitude
            
            # Generate exactly 10 ticks
            start_tick = np.ceil(nm_min / nice_step) * nice_step
            nm_ticks = [start_tick + i * nice_step for i in range(10)]
            nm_ticks = [nm for nm in nm_ticks if nm_min <= nm <= nm_max * 1.1][:10]
            
            # Convert nm ticks back to eV positions
            ev_positions = [1239.84 / nm for nm in nm_ticks if nm > 0]
            
            # Set custom ticks and labels
            ax2_top.set_xlim(ax.get_xlim())
            ax2_top.set_xticks(ev_positions)
            
            # Format labels
            if max(nm_ticks) > 10000:
                labels = [f'{nm/1000:.1f}k' for nm in nm_ticks]
            else:
                labels = [f'{nm:.0f}' for nm in nm_ticks]
            
            ax2_top.set_xticklabels(labels)
            ax2_top.set_xlabel('Wavelength [nm]', fontsize=18)
            ax2_top.tick_params(axis='x', rotation=45, labelsize=14)
            
            # Add vertical grid lines
            for ev_pos in ev_positions:
                ax.axvline(x=ev_pos, color='lightblue', linestyle='--', alpha=0.7, linewidth=0.8)
            
            # Add legend
            ax.legend(loc='upper right', bbox_to_anchor=(0.75, 1), 
                     frameon=False, fancybox=False, shadow=False, fontsize=16)
            
            # Set tick font sizes
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
        
        # Create both subplots
        create_subplot(ax1, 'axial', 'Axial Ligand Combinations')
        
        create_subplot(ax2, 'charge', 'Charge-Multiplicity Combinations')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the combined plot
        if save:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            filename = "homo-lumo-gap_histogram.png"
            fig.savefig(os.path.join(save_dir, filename), bbox_inches='tight')
            print(f"Combined HOMO-LUMO plot saved as {filename}")
        
        plt.show()

    def plot_lft_orbitals_histograms(self, df, keyword='lft_orbitals', bins=10,
                                      save=False, save_dir='.', group_by_charge=False):
        """
        Plot histograms for columns containing a specific keyword (e.g., 'lft_orbitals')
        or 'lft_occupancy_delta', grouped by axial1/axial2 or charge/multiplicity.
        Uses linear y-axis, aligns x-axes to same data range across subplots.
        """
        # Decoding dictionary for encoded categorical variables
        REPLACEMENTS = {
            'axial1': {0: 'CYS', 1: 'HIS'},
            'axial2': {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'},
            'function': {
                0: 'electron transporter',
                1: 'hydrolase',
                2: 'oxidoreductase',
                3: 'oxygen carrier',
                4: 'transferase',
                5: 'unclassified'
            }
        }
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2'])
        group_vars = ('charge', 'multiplicity') if group_by_charge else ('axial1', 'axial2')

        if save and not os.path.exists(save_dir):
            os.makedirs(save_dir)

        numeric_cols = df_dec.select_dtypes(include=[np.number]).columns
        lft_cols = [col for col in numeric_cols if keyword in col or col == 'lft_occupancy_delta']

        combos = df_dec[list(group_vars)].drop_duplicates().to_records(index=False)
        n = len(combos)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))

        # Legend proxies for stats lines
        from matplotlib.lines import Line2D
        legend_proxies = [
            Line2D([0], [0], color='red', linestyle='--', linewidth=1, label='Mean'),
            Line2D([0], [0], color='blue', linestyle=':', linewidth=1, label='25th/75th'),
            Line2D([0], [0], color='green', linestyle='-.', linewidth=1, label='Median')
        ]

        for col in lft_cols:
            # compute global x-range
            all_data = df_dec[col][df_dec[col] != -1].dropna()
            if all_data.empty:
                continue
            x_min, x_max = all_data.min(), all_data.max()

            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(4*ncols, 3*nrows),
                                     sharex=True)
            axes_flat = axes.flatten() if isinstance(axes, np.ndarray) else [axes]

            for ax, combo in zip(axes_flat, combos):
                val1, val2 = combo
                mask = (df_dec[group_vars[0]] == val1) & (df_dec[group_vars[1]] == val2) & (df_dec[col] != -1)
                data = df_dec.loc[mask, col].dropna()

                ax.set_xlim(x_min, x_max)

                if data.empty:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                else:
                    ax.hist(data, bins=bins, edgecolor='black')

                    mean_val = data.mean()
                    q25, q50, q75 = data.quantile([0.25, 0.5, 0.75])
                    ax.axvline(mean_val, color='red', linestyle='--', linewidth=1)
                    ax.axvline(q25, color='blue', linestyle=':', linewidth=1)
                    ax.axvline(q50, color='green', linestyle='-.', linewidth=1)
                    ax.axvline(q75, color='blue', linestyle=':', linewidth=1)


            for ax in axes_flat[n:]:
                ax.set_visible(False)

            grouping_name = 'charge_multiplicity' if group_by_charge else 'axial'
            fig.legend(handles=legend_proxies, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
            plt.tight_layout(rect=(0, 0.03, 1, 0.95))
            if save:
                fname = f"{col}_lft_orbital_{grouping_name}_linear.png"
                fig.savefig(os.path.join(save_dir, fname), bbox_inches='tight')
            plt.close(fig)

    def plot_dual_lft_scatter_with_violin(self, df1, df2,
                                        x_col='lft_energy_delta',
                                        y_col='iron_natural_charge',
                                        color_by_left='axial',  # Color scheme for left plot
                                        color_by_right='charge_multiplicity',  # Color scheme for right plot
                                        color_by_violin='axial',  # Color scheme for violin plot
                                        save=False,
                                        save_dir='.',
                                        palette=None,
                                        figsize=(20, 8),
                                        marker_size=30,
                                        draw_dividers=True,
                                        light_lines=False,
                                        number_quadrants=False,
                                        violin_alpha=0.7,
                                        violin_width=0.8,
                                        color_grading=True,
                                        gray_background_condition=None,
                                        titles=None,
                                        convert_to_ev=False):
        """
        Create dual scatter plots with a violin plot in the middle showing density distributions.
        Optionally highlights background regions in light gray based on specified conditions.

        Parameters
        ----------
        df1, df2 : pandas.DataFrame
            DataFrames for left and right scatter plots
        x_col, y_col : str
            Column names for x and y axes
        color_by_left : str
            Grouping strategy for left plot ('axial', 'charge_multiplicity', or 'both')
        color_by_right : str
            Grouping strategy for right plot ('axial', 'charge_multiplicity', or 'both')
        color_by_violin : str
            Grouping strategy for violin plot ('axial', 'charge_multiplicity', or 'both')
        save : bool
            Whether to save the figure
        save_dir : str
            Directory to save the figure
        palette : list
            Color palette to use
        figsize : tuple
            Figure size (width, height) in inches
        marker_size : int
            Size of scatter points
        draw_dividers : bool
            Whether to draw reference lines
        light_lines : bool
            Whether to use light colored reference lines
        number_quadrants : bool
            Whether to number quadrants/regions
        violin_alpha : float
            Transparency of violin plots (0-1)
        violin_width : float
            Width scaling factor for violin plots
        color_grading : bool
            Whether to apply color intensity grading
        gray_background_condition : callable or None
            Function that takes (x, y) and returns True for gray background regions.
            Default: lambda x, y: (x < 0) or (y < 0)
        titles : list or None
            List of titles for [left_plot, violin_plot, right_plot]. If None, uses defaults.
        convert_to_ev : bool
            If True, convert x_col values from Hartree to eV (multiply by 27.211386245988).
            Default: False (keep values in Hartree)
        """
        # Decoding dictionary for encoded categorical variables
        REPLACEMENTS = {
            'axial1': {0: 'CYS', 1: 'HIS'},
            'axial2': {0: 'HIS', 1: 'HOH', 2: 'MET', 3: 'OXY'},
            'function': {
                0: 'electron transporter',
                1: 'hydrolase',
                2: 'oxidoreductase',
                3: 'oxygen carrier',
                4: 'transferase',
                5: 'unclassified'
            }
        }
        
        # Color palette for categorical plotting (colorblind-friendly)
        COLOR_PALETTE = EXTENDED_COLORBLIND_PALETTE
        MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X']
        
        # Default gray background condition: x < 0 or y < 0
        if gray_background_condition is None:
            gray_background_condition = lambda x, y: (x < 0) or (y < 0)
        
        # Default titles
        if titles is None:
            titles = ['Left Dataset', 'Density\nDistribution', 'Right Dataset']
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        df1_dec = _decode_df(df1, vars_to_decode=['axial1', 'axial2'])
        df2_dec = _decode_df(df2, vars_to_decode=['axial1', 'axial2'])
        pal = palette or COLOR_PALETTE

        # Prepare data for both datasets
        cols = [x_col, y_col, 'axial1', 'axial2', 'charge', 'multiplicity', 'total_energies[0]']
        
        df1_plot = df1_dec.loc[
            (df1_dec[x_col] != -1) & (df1_dec[y_col] != -1), cols
        ].dropna()
        df2_plot = df2_dec.loc[
            (df2_dec[x_col] != -1) & (df2_dec[y_col] != -1), cols
        ].dropna()

        df1_plot['axial_combo'] = df1_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
        df1_plot['charge_combo'] = df1_plot.apply(lambda r: f"{r.charge}-{r.multiplicity}", axis=1)
        df2_plot['axial_combo'] = df2_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
        df2_plot['charge_combo'] = df2_plot.apply(lambda r: f"{r.charge}-{r.multiplicity}", axis=1)

        # Convert x_col from Hartree to eV if requested
        if convert_to_ev:
            HARTREE_TO_EV = 27.211386245988  # CODATA 2018 value
            df1_plot[x_col] = df1_plot[x_col] * HARTREE_TO_EV
            df2_plot[x_col] = df2_plot[x_col] * HARTREE_TO_EV

        # Combine datasets to get consistent color mapping
        combined_df = pd.concat([df1_plot, df2_plot], ignore_index=True)

        # Create color mappings for all schemes using centralized color maps
        groups_axial = sorted(combined_df['axial_combo'].unique())
        color_map_axial = get_axial_ligand_colors(groups_axial)

        groups_charge = sorted(combined_df['charge_combo'].unique())
        color_map_charge = get_charge_multiplicity_colors(groups_charge)
        marker_map = {g: MARKERS[i % len(MARKERS)] for i, g in enumerate(groups_charge)}

        # Create figure with subplots: left scatter, violin, right scatter
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 3, width_ratios=[3, 1, 3], wspace=0.15)
        
        # Left scatter plot
        ax_left = fig.add_subplot(gs[0, 0])
        
        # Violin plot (middle)
        ax_violin = fig.add_subplot(gs[0, 1])
        
        # Right scatter plot
        ax_right = fig.add_subplot(gs[0, 2])

        # Create marker mapping for axial combinations
        axial_marker_map = {
            'CYS-HOH': 'o',    # circle
            'HIS-HIS': 's',    # square
            'HIS-HOH': '^',    # triangle up
            'HIS-MET': 'D',    # diamond
            'HIS-OXY': 'v'     # triangle down
        }
        
        # Function to add gray background
        def add_gray_background(ax, df_plot):
            if len(df_plot) == 0:
                return
            
            # Get axis limits
            x_min, x_max = df_plot[x_col].min(), df_plot[x_col].max()
            y_min, y_max = df_plot[y_col].min(), df_plot[y_col].max()
            
            # Expand limits slightly for background
            x_range = x_max - x_min
            y_range = y_max - y_min
            x_min -= 0.1 * x_range
            x_max += 0.1 * x_range
            y_min -= 0.1 * y_range
            y_max += 0.1 * y_range
            
            # Create a grid to determine gray regions
            x_grid = np.linspace(x_min, x_max, 200)
            y_grid = np.linspace(y_min, y_max, 200)
            X, Y = np.meshgrid(x_grid, y_grid)
            
            # Apply condition to determine gray regions
            gray_mask = np.zeros_like(X, dtype=bool)
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    gray_mask[i, j] = gray_background_condition(X[i, j], Y[i, j])
            
            # Add gray background where condition is True
            ax.contourf(X, Y, gray_mask.astype(int), levels=[0.5, 1.5], colors=['lightgray'], alpha=0.3)
        
        # Function to calculate energy stats and plot scatter points
        def plot_scatter_data(ax, df_plot, title, color_by_scheme):
            import matplotlib.colors as mcolors
            
            # Add gray background first
            add_gray_background(ax, df_plot)
            
            # Calculate energy stats for color grading
            axial_combinations = df_plot['axial_combo'].unique()
            energy_stats = {}
            
            for combo in axial_combinations:
                combo_data = df_plot[df_plot['axial_combo'] == combo]
                if not combo_data.empty:
                    energy_min = combo_data['total_energies[0]'].min()
                    energy_max = combo_data['total_energies[0]'].max()
                    energy_range = energy_max - energy_min
                    energy_stats[combo] = {
                        'min': energy_min,
                        'max': energy_max,
                        'range': energy_range
                    }
            
            # Plot scatter points
            for _, row in df_plot.iterrows():
                if color_by_scheme == 'axial':
                    base_color = color_map_axial[row['axial_combo']]
                    m = axial_marker_map.get(row['axial_combo'], 'o')
                elif color_by_scheme == 'charge_multiplicity':
                    base_color = color_map_charge[row['charge_combo']]
                    m = marker_map[row['charge_combo']]
                else:  # 'both'
                    base_color = color_map_axial[row['axial_combo']]
                    m = axial_marker_map.get(row['axial_combo'], 'o')
                
                # Calculate intensity based on total_energy
                if color_grading:
                    axial_combo = row['axial_combo']
                    if axial_combo in energy_stats and energy_stats[axial_combo]['range'] > 0:
                        stats = energy_stats[axial_combo]
                        # Amplify small energy differences
                        amplified_diff = 10 * (row['total_energies[0]'] - stats['min'])
                        amplified_range = 10 * stats['range']
                        normalized_energy = amplified_diff / amplified_range
                        normalized_energy = min(1.0, max(0.0, normalized_energy))
                        intensity = 0.1 + 0.9 * (normalized_energy ** 2)
                    else:
                        intensity = 1.0
                    
                    # Convert base color to RGB and apply intensity
                    base_rgb = mcolors.to_rgb(base_color)
                    final_color = tuple(intensity * c + (1 - intensity) * 1.0 for c in base_rgb)
                else:
                    final_color = base_color
                
                ax.scatter(row[x_col], row[y_col], color=final_color, marker=m,
                          alpha=0.8, edgecolor='black', s=marker_size)
            
            # Optional reference lines
            if draw_dividers:
                line_color = 'lightgrey' if light_lines else 'black'
                ax.axvline(0, color=line_color, linestyle='--', linewidth=1)
                for y in (0.0, 0.3, 0.6, 0.9):
                    ax.axhline(y, color=line_color, linestyle='--', linewidth=1)

                if number_quadrants:
                    xlim = ax.get_xlim()
                    ylim = ax.get_ylim()
                    x_cuts = [xlim[0], 0, xlim[1]]
                    y_cuts = [ylim[0], 0.2, 0.6, 0.85, ylim[1]]
                    label = 1
                    for j in range(len(y_cuts)-1):
                        for i in range(len(x_cuts)-1):
                            x_center = (x_cuts[i] + x_cuts[i+1]) / 2
                            y_center = (y_cuts[j] + y_cuts[j+1]) / 2
                            ax.text(x_center, y_center, str(label),
                                   ha='center', va='center', fontsize=14, alpha=0.7)
                            label += 1
            
            # Set xlabel based on convert_to_ev flag
            x_unit = "eV" if convert_to_ev else "Hartree"
            ax.set_xlabel(f"Δ [{x_unit}]", fontsize=16)
            if ax == ax_left:
                ax.set_ylabel("Iron Charge (NBO) [a.u.]", fontsize=16)
        
        # Plot both scatter plots with different color schemes
        plot_scatter_data(ax_left, df1_plot, "", color_by_left)
        plot_scatter_data(ax_right, df2_plot, "", color_by_right)
        
        # Ensure consistent y-axis limits
        y_min = min(df1_plot[y_col].min(), df2_plot[y_col].min()) if len(df1_plot) > 0 and len(df2_plot) > 0 else 0
        y_max = max(df1_plot[y_col].max(), df2_plot[y_col].max()) if len(df1_plot) > 0 and len(df2_plot) > 0 else 1
        y_range = y_max - y_min
        y_min -= 0.05 * y_range
        y_max += 0.05 * y_range
        
        ax_left.set_ylim(y_min, y_max)
        ax_right.set_ylim(y_min, y_max)
        ax_violin.set_ylim(y_min, y_max)
        
        # Create violin plots from combined data
        if color_by_violin == 'axial':
            grouping_col = 'axial_combo'
            groups = groups_axial
            color_map = color_map_axial
        elif color_by_violin == 'charge_multiplicity':
            grouping_col = 'charge_combo'
            groups = groups_charge
            color_map = color_map_charge
        else:  # both
            grouping_col = 'axial_combo'
            groups = groups_axial
            color_map = color_map_axial

        # Prepare data for violin plots
        violin_data = []
        violin_labels = []
        violin_colors = []
        
        # Add group-specific violin plots
        for group in sorted(groups):
            group_data = combined_df[combined_df[grouping_col] == group][y_col].dropna()
            if len(group_data) >= 2:
                violin_data.append(group_data)
                violin_labels.append(group)
                violin_colors.append(color_map[group])
        
        # Add total density violin plot
        total_data = combined_df[y_col].dropna()
        if len(total_data) >= 2:
            violin_data.append(total_data)
            violin_labels.append('Total')
            violin_colors.append('gray')

        # Create violin plots
        if len(violin_data) > 0:
            positions = np.arange(len(violin_data))
            
            violin_parts = ax_violin.violinplot(violin_data, positions=positions, 
                                              widths=violin_width, vert=True)
            
            # Color the violin plots
            for i, (pc, color) in enumerate(zip(violin_parts['bodies'], violin_colors)):
                pc.set_facecolor(color)
                pc.set_alpha(violin_alpha)
                pc.set_edgecolor('black')
                pc.set_linewidth(1)
            
            # Style the violin plot components
            for partname in ('cbars', 'cmins', 'cmaxes'):
                if partname in violin_parts:
                    violin_parts[partname].set_color('black')
                    violin_parts[partname].set_linewidth(1)
            
            ax_violin.set_xticks(positions)
            ax_violin.set_xticklabels(violin_labels, rotation=45, ha='right')
            if color_by_violin == 'axial':
                ax_violin.set_xlabel('Axial Ligands', fontsize=14)
            elif color_by_violin == 'charge_multiplicity':
                ax_violin.set_xlabel('Electronic State', fontsize=14)
            else:
                ax_violin.set_xlabel('Groups', fontsize=14)

        # Add reference lines to violin plot
        if draw_dividers:
            line_color = 'lightgrey' if light_lines else 'black'
            for y in (0.0, 0.3, 0.6, 0.9):
                ax_violin.axhline(y, color=line_color, linestyle='--', linewidth=1, alpha=0.5)
        
        ax_violin.tick_params(axis='y', which='both', labelleft=False, labelright=False)
        
        # Remove y-axis labels from right plot
        ax_right.tick_params(axis='y', which='both', labelleft=False, labelright=False)
        
        # Create legends for both plots
        from matplotlib.lines import Line2D
        
        # Legend for left plot (axial)
        legend_handles_left = []
        if color_by_left in ('axial', 'both'):
            legend_handles_left += [mpatches.Patch(color=color_map_axial[g], label=g)
                                   for g in color_map_axial]
        if color_by_left in ('charge_multiplicity', 'both'):
            legend_handles_left += [Line2D([0], [0], marker=marker_map[g], color='w',
                                          markerfacecolor=color_map_charge[g], markersize=8,
                                          markeredgecolor='black', label=g)
                                   for g in marker_map]
        
        legend_title_left = None
        if color_by_left == 'charge_multiplicity':
            legend_title_left = 'Charge-Multiplicity'
        elif color_by_left == 'axial':
            legend_title_left = 'Axial Ligands'
        elif color_by_left == 'both':
            legend_title_left = 'Groups'
        
        ax_left.legend(handles=legend_handles_left, loc='upper right',
                      bbox_to_anchor=(0.98, 0.98), frameon=False, fancybox=False, shadow=False,
                      title=legend_title_left, fontsize=10, title_fontsize=12)
        
        # Legend for right plot (charge_multiplicity)
        legend_handles_right = []
        if color_by_right in ('axial', 'both'):
            legend_handles_right += [mpatches.Patch(color=color_map_axial[g], label=g)
                                    for g in color_map_axial]
        if color_by_right in ('charge_multiplicity', 'both'):
            legend_handles_right += [Line2D([0], [0], marker=marker_map[g], color='w',
                                           markerfacecolor=color_map_charge[g], markersize=8,
                                           markeredgecolor='black', label=g)
                                    for g in marker_map]
        
        legend_title_right = None
        if color_by_right == 'charge_multiplicity':
            legend_title_right = 'Charge-Multiplicity'
        elif color_by_right == 'axial':
            legend_title_right = 'Axial Ligands'
        elif color_by_right == 'both':
            legend_title_right = 'Groups'
        
        ax_right.legend(handles=legend_handles_right, loc='upper right',
                       bbox_to_anchor=(0.98, 0.98), frameon=False, fancybox=False, shadow=False,
                       title=legend_title_right, fontsize=10, title_fontsize=12)

        # Set tick font sizes for all subplots
        ax_left.tick_params(axis='x', labelsize=12)
        ax_left.tick_params(axis='y', labelsize=12)
        ax_violin.tick_params(axis='x', labelsize=12)
        ax_right.tick_params(axis='x', labelsize=12)

        plt.tight_layout()
        if save:
            os.makedirs(save_dir, exist_ok=True)
            unit_suffix = "_eV" if convert_to_ev else "_Hartree"
            fname = f"plots/{x_col}_{y_col}{unit_suffix}.png"
            fig.savefig(os.path.join(save_dir, fname), bbox_inches='tight', dpi=300)
        plt.show()
        plt.close(fig)

    def plot_iron_d_orbital_energy_differences(self, df, group_by='charge_multiplicity', 
                                               save_dir='plots', decode_dicts=None,
                                               figsize=(16, 12), show_groupings='separate', verbose=False):
        """
        Plot iron d-orbital distortion-type counts using two stacked barplots.
        
        The figure contains only two panels:
        1. Distortion-type counts stacked by charge-multiplicity
        2. Distortion-type counts stacked by axial ligands

        Colors are assigned from the centralized color-blind mappings
        (`CHARGE_MULTIPLICITY_COLOR_MAP` and `AXIAL_LIGAND_COLOR_MAP`) so
        this plot stays consistent with the rest of the codebase.
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing orbital energy difference columns
        group_by : str
            Retained for backward compatibility (not used for panel selection)
        save_dir : str
            Directory to save plots
        decode_dicts : dict
            Dictionary for decoding categorical variables
        figsize : tuple
            Figure size (width, height)
        show_groupings : str
            Retained for backward compatibility; the function always renders
            only the two top grouped barplots
        verbose : bool
            Whether to print detailed statistics and debugging information
        """
        os.makedirs(save_dir, exist_ok=True)
        
        # Enable proper math text rendering for subscripts and superscripts
        import matplotlib
        matplotlib.rcParams['text.usetex'] = False  # Use matplotlib's built-in math renderer
        matplotlib.rcParams['mathtext.default'] = 'regular'
        matplotlib.rcParams['font.family'] = 'DejaVu Sans'  # Use font that supports math text
        
        if decode_dicts is not None and not df.attrs.get("decoded", False):
            df = self.preprocess_dataframe(df, decode_dicts)
        
        # Create charge_multiplicity column if not present
        if 'charge_multiplicity' not in df.columns:
            df = self.create_charge_multiplicity_column(df)
        
        # Extract PDB-ID from filename if not present
        if 'PDB-ID' not in df.columns and 'file_name' in df.columns:
            df['PDB-ID'] = df['file_name'].str[:4]
        
        # Create axials column if not present
        if 'axials' not in df.columns:
            if 'axial1' in df.columns and 'axial2' in df.columns:
                df['axials'] = df['axial1'].astype(str) + '-' + df['axial2'].astype(str)
            else:
                df['axials'] = 'Unknown'
        
        # Define the orbital energy difference columns
        energy_diff_cols = [
            'dx2y2_dz2_diff', 'dz2_dxy_diff', 'dxy_dxz_dyz_diff', 'eg_t2g_diff'
        ]
        
        # Check if required columns exist
        missing_cols = [col for col in energy_diff_cols if col not in df.columns]
        if missing_cols:
            if verbose:
                print(f"Missing required orbital energy difference columns: {missing_cols}")
                print("Available columns with 'diff' in name:")
                diff_cols = [col for col in df.columns if 'diff' in col.lower()]
                for col in diff_cols[:10]:  # Show first 10 to avoid overwhelming output
                    print(f"  - {col}")
            return
        
        # Filter out rows where orbital differences are missing
        df_valid = df.dropna(subset=energy_diff_cols)
        
        if len(df_valid) == 0:
            if verbose:
                print("No valid orbital energy difference data found after filtering missing values")
            return
        
        # Keep this plot focused on the two top grouping panels only:
        # charge-multiplicity and axial ligand stacked barplots.
        fig, (ax_charge, ax_axial) = plt.subplots(1, 2, figsize=figsize)
        fig.subplots_adjust(wspace=0.4)
        if verbose and show_groupings != 'separate':
            print(f"DEBUG: show_groupings='{show_groupings}' ignored; rendering two grouped panels only.")
        
        # Add comprehensive column validation and debugging
        if verbose:
            print(f"DEBUG: Total rows in df_valid: {len(df_valid)}")
            print(f"DEBUG: Available columns: {list(df_valid.columns)}")
            
            # Check key columns
            key_columns = ['orbital_distortion_type', 'charge_multiplicity', 'axials', 'axial1', 'axial2']
            for col in key_columns:
                if col in df_valid.columns:
                    print(f"DEBUG: Column '{col}' found, sample values: {df_valid[col].dropna().head(3).tolist()}")
                    print(f"DEBUG: Column '{col}' unique count: {df_valid[col].nunique()}")
                else:
                    print(f"DEBUG: Column '{col}' NOT found")
        
        # Check for orbital distortion column with flexible naming
        distortion_col = None
        possible_distortion_names = ['orbital_distortion_type', 'distortion_type', 'orbital_distortion', 'distortion']
        for col_name in possible_distortion_names:
            if col_name in df_valid.columns:
                distortion_col = col_name
                break
        
        # Decode distortion types if available
        # Support both numeric and string distortion type values
        numeric_distortion_map = {0: 'Octahedral', 1: 'Tetr. Elong.', 2: 'Tetr. Compr.', 3: 'Square Planar'}
        string_distortion_map = {
            'octahedral': 'Octahedral', 
            'tetrag_elong': 'Tetr. Elong.', 
            'tetrag_compr': 'Tetr. Compr.', 
            'square_planar': 'Square Planar'
        }
        has_distortion_data = False
        
        if distortion_col:
            # Check if values are numeric or string and apply appropriate mapping
            sample_value = df_valid[distortion_col].dropna().iloc[0] if len(df_valid[distortion_col].dropna()) > 0 else None
            if sample_value is not None:
                if pd.api.types.is_numeric_dtype(df_valid[distortion_col]):
                    df_valid['distortion_decoded'] = df_valid[distortion_col].map(numeric_distortion_map)
                else:
                    # Convert to lowercase for case-insensitive matching
                    df_valid['distortion_decoded'] = df_valid[distortion_col].str.lower().map(string_distortion_map)
            else:
                df_valid['distortion_decoded'] = 'Unknown'
            has_distortion_data = True
            if verbose:
                print(f"DEBUG: Using distortion column: {distortion_col}")
                print(f"DEBUG: Distortion value counts: {df_valid['distortion_decoded'].value_counts().to_dict()}")
        else:
            if verbose:
                print("DEBUG: No distortion column found, will create fallback categories")
            # Create a fallback single category for plotting
            df_valid['distortion_decoded'] = 'All Types'
            has_distortion_data = False
        
        # Helper functions for consistent ordering and centralized color mapping.
        def charge_mult_sort_key(value):
            match = re.match(r'^\s*(-?\d+(?:\.\d+)?)\D+(-?\d+(?:\.\d+)?)\s*$', str(value))
            if not match:
                return (float('inf'), float('inf'), str(value))
            return (float(match.group(1)), float(match.group(2)), str(value))

        def get_ordered_groups(values, group_col):
            unique_groups = sorted(pd.Series(values).dropna().astype(str).unique())
            if group_col == 'axials':
                preferred = [g for g in AXIAL_LIGAND_COLOR_MAP.keys() if g in unique_groups]
                remaining = [g for g in unique_groups if g not in preferred]
                return preferred + remaining
            if group_col == 'charge_multiplicity':
                return sorted(unique_groups, key=charge_mult_sort_key)
            return unique_groups

        def get_group_color_map(groups, group_col):
            if group_col == 'axials':
                return get_axial_ligand_colors(groups, extended=True)
            if group_col == 'charge_multiplicity':
                return get_charge_multiplicity_colors(groups, extended=True)
            return {group: color for group, color in zip(groups, get_color_sequence(len(groups)))}

        # Helper function to create stacked bar plot for groupings.
        def create_stacked_barplot(ax, distortion_counts, group_col):
            """Create stacked barplot showing breakdown by group_col for each distortion type"""
            unique_groups = get_ordered_groups(df_valid[group_col], group_col)
            
            verbose=False
            if verbose:
                print(f"DEBUG: Creating stacked barplot for {group_col}")
                print(f"DEBUG: Unique groups found: {unique_groups}")
                print(f"DEBUG: Distortion counts: {distortion_counts}")
                print(f"DEBUG: Distortion labels: {list(distortion_counts.index)}")
            
            if len(unique_groups) == 0:
                print(f"DEBUG: No unique groups found for {group_col}")
                return
            
            group_color_map = get_group_color_map(unique_groups, group_col)
            
            # Create data structure for stacked bars
            stacked_data = {}
            for distortion in distortion_counts.index:
                stacked_data[distortion] = {}
                subset = df_valid[df_valid['distortion_decoded'] == distortion]
                group_counts = subset[group_col].value_counts()
                if verbose:
                    print(f"DEBUG: For {distortion}, subset size: {len(subset)}")
                    print(f"DEBUG: Group counts: {group_counts.to_dict()}")
                for group in unique_groups:
                    stacked_data[distortion][group] = group_counts.get(group, 0)
            
            # FIXED: Use numeric positions for bars
            x_pos = np.arange(len(distortion_counts))
            distortion_labels = list(distortion_counts.index)
            
            if verbose:
                print(f"DEBUG: x_pos: {x_pos}")
                print(f"DEBUG: distortion_labels: {distortion_labels}")
            
            # Create stacked bars with proper positioning
            bottom = np.zeros(len(distortion_counts))
            for i, group in enumerate(unique_groups):
                values = [stacked_data[dist][group] for dist in distortion_counts.index]
                if verbose:
                    print(f"DEBUG: Group {group}, values: {values}, bottom: {bottom}")
                if sum(values) > 0:  # Only plot if there's data
                    bars = ax.bar(x_pos, values, bottom=bottom, width=0.6,
                                 color=group_color_map[group], label=str(group), alpha=0.8)
                    if verbose:
                        print(f"DEBUG: Created bars for {group}: {[bar.get_height() for bar in bars]}")
                    bottom += np.array(values)
            
            # Set proper x-axis labels
            ax.set_xticks(x_pos)
            ax.set_xticklabels(distortion_labels, rotation=45, fontsize=12)
            ax.set_ylabel('Count', fontsize=14)
            ax.tick_params(axis='y', labelsize=12)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, title_fontsize=14)
            
            # Add total count labels on bars
            total_all_distortions = distortion_counts.sum()
            for i, (distortion, total) in enumerate(distortion_counts.items()):
                percentage = (100.0 * total / total_all_distortions) if total_all_distortions > 0 else 0.0
                ax.annotate(f'{int(total)} ({percentage:.1f}%)',
                           xy=(x_pos[i], total),
                           xytext=(0, 3), textcoords="offset points",
                           ha='center', va='bottom', fontsize=14)
            
            # Set ylim to leave space at the top (reduce height but leave vertical space)
            max_val = distortion_counts.max() if not distortion_counts.empty else 1
            ax.set_ylim(0, max_val * 1.4)  # extra space for count + percentage labels

            if verbose:               
                print(f"DEBUG: Final bottom values: {bottom}")
        
        # Validate required columns for grouping plots
        required_grouping_cols = ['charge_multiplicity', 'axials']
        missing_grouping_cols = [col for col in required_grouping_cols if col not in df_valid.columns]
        
        if missing_grouping_cols:
            if verbose: 
                print(f"DEBUG: Missing required grouping columns: {missing_grouping_cols}")
                print("DEBUG: Cannot create top grouping plots.")
            for ax in (ax_charge, ax_axial):
                ax.text(0.5, 0.5, 'Grouping data not available', ha='center', va='center',
                        transform=ax.transAxes)
                ax.set_axis_off()
        else:
            # We have the required grouping columns, proceed with barplots
            distortion_counts = df_valid['distortion_decoded'].value_counts()
            
            if verbose:
                print(f"DEBUG: Distortion counts for plotting: {distortion_counts.to_dict()}")
                print(f"DEBUG: Charge-multiplicity counts: {df_valid['charge_multiplicity'].value_counts().to_dict()}")
                print(f"DEBUG: Axials counts: {df_valid['axials'].value_counts().to_dict()}")
            create_stacked_barplot(ax_charge, distortion_counts, 'charge_multiplicity')
            ax_charge.set_title('Charge-Multiplicity Breakdown', fontsize=16)
            create_stacked_barplot(ax_axial, distortion_counts, 'axials')
            ax_axial.set_title('Axial Ligand Breakdown', fontsize=16)
        
        
        # Save plot
        filename = f'{save_dir}/iron_d-orbs_symmetry.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()
        
        if verbose: print(f"Iron d-orbital energy analysis plot saved to: {filename}")
        
        if verbose:
            # Print summary statistics
            print("\nSummary Statistics:")
            print(f"Total structures analyzed: {len(df_valid)}")
            if 'orbital_distortion_type' in df_valid.columns:
                print("\nDistortion type distribution:")
                for distortion, count in df_valid['distortion_decoded'].value_counts().items():
                    print(f"  {distortion}: {count}")
                
                print("\nCharge-multiplicity breakdown by distortion type:")
                for distortion in df_valid['distortion_decoded'].unique():
                    subset = df_valid[df_valid['distortion_decoded'] == distortion]
                    charge_counts = subset['charge_multiplicity'].value_counts()
                    print(f"  {distortion}:")
                    for charge, count in charge_counts.items():
                        print(f"    {charge}: {count}")
                
                print("\nAxial ligand breakdown by distortion type:")
                for distortion in df_valid['distortion_decoded'].unique():
                    subset = df_valid[df_valid['distortion_decoded'] == distortion]
                    axial_counts = subset['axials'].value_counts()
                    print(f"  {distortion}:")
                    for axial, count in axial_counts.head(5).items():  # Top 5 most common
                        print(f"    {axial}: {count}")
            
            print(f"\nEnergy difference ranges (a.u.):")
            for col in energy_diff_cols:
                if col in df_valid.columns:
                    col_data = df_valid[col].dropna()
                    print(f"  {col}: {col_data.min():.4f} to {col_data.max():.4f} (mean: {col_data.mean():.4f})")

    def plot_reduction_potentials_comprehensive(self, reduction_potentials_file='tables/reduction_potentials.csv',
                                            processed_data_file='tables/processed_output.csv',
                                            primary_grouping='charge_multiplicity', save_dir='plots',
                                            decode_dicts=None, figsize=(20, 16),
                                            potential_range=None, focus_pdbs=None,
                                            extra_plots=False):
        """
        Creates visualizations of reduction potential data.
        Optionally includes extra comparison plots if `extra_plots=True`.
        """
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        import os
        verbose=False
        os.makedirs(save_dir, exist_ok=True)

        try:
            # Load data
            redox_df = pd.read_csv(reduction_potentials_file)
            processed_df = pd.read_csv(processed_data_file)
            
            # Extract structure ID
            processed_df['structure_id'] = processed_df['file_name'].str[:4]
            
            # Merge metadata
            merged_df = redox_df.merge(
                processed_df[['structure_id', 'charge', 'multiplicity', 'axial1', 'axial2', 'function', 'PDB-ID']].drop_duplicates('structure_id'),
                on='structure_id', how='left'
            )
            
            # Create grouping columns
            merged_df['charge_multiplicity'] = merged_df['charge'].astype(str) + '_' + merged_df['multiplicity'].astype(str)
            if 'axial1' in merged_df.columns and 'axial2' in merged_df.columns:
                merged_df['axials'] = merged_df['axial1'].astype(str) + '-' + merged_df['axial2'].astype(str)
            
            # Decode categories if provided
            if decode_dicts is not None:
                for col, mapping in decode_dicts.items():
                    if col in merged_df.columns:
                        merged_df[col + '_decoded'] = merged_df[col].map(mapping)
                        if col in ['axial1', 'axial2'] and 'axial1_decoded' in merged_df.columns and 'axial2_decoded' in merged_df.columns:
                            merged_df['axials'] = merged_df['axial1_decoded'].astype(str) + '-' + merged_df['axial2_decoded'].astype(str)
            
            # Apply potential filter
            if potential_range is not None:
                merged_df = merged_df[
                    (merged_df['E°approx (V)'] >= potential_range[0]) &
                    (merged_df['E°approx (V)'] <= potential_range[1])
                ]
            
            # Filter extremes
            merged_df = merged_df[
                (merged_df['E°approx (V)'] < 5.0) & 
                (merged_df['E°approx (V)'] >= -5.0)
            ]
            
            # Create figure layout
            if extra_plots:
                fig = plt.figure(figsize=figsize)
                gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.3)
                ax1 = fig.add_subplot(gs[0, :])
            else:
                fig, ax1 = plt.subplots(figsize=(figsize[0], figsize[1] / 2))
            
            # Color mapping
            charge_groups = sorted(merged_df['charge_multiplicity'].dropna().unique(), key=charge_multiplicity_sort_key)
            colors_cm = get_charge_multiplicity_colors(charge_groups, extended=True)
            colors_axial = get_color_sequence(len(merged_df['axials'].unique()) if 'axials' in merged_df.columns else 4, offset=2)
            colors_function = get_color_sequence(len(merged_df['function'].unique()) if 'function' in merged_df.columns else 6, offset=4)
            
            # Determine grouping variable
            if primary_grouping == 'charge_multiplicity':
                scatter_colors = colors_cm
                group_var = 'charge_multiplicity'
            elif primary_grouping == 'axials' and 'axials' in merged_df.columns:
                scatter_colors = colors_axial
                group_var = 'axials'
            else:
                scatter_colors = colors_function
                group_var = 'function'
            
            # Sort samples by axial ligand combination for better ordering
            if 'axials' in merged_df.columns:
                # Define preferred order for axial ligand combinations
                axial_order = ['HIS-HOH', 'HIS-HIS', 'HIS-OXY', 'HIS-MET', 'CYS-HOH']
                # Get unique axials in data and sort them according to preferred order
                unique_axials = merged_df['axials'].unique()
                sorted_axials = [ax for ax in axial_order if ax in unique_axials]
                # Add any axials not in the preferred order at the end
                sorted_axials.extend([ax for ax in unique_axials if ax not in axial_order])
                
                # Create a categorical column for sorting
                merged_df['axials_cat'] = pd.Categorical(merged_df['axials'], categories=sorted_axials, ordered=True)
                # Sort the dataframe by axial ligand combinations, then by charge_multiplicity
                merged_df_sorted = merged_df.sort_values(['axials_cat', 'charge_multiplicity']).reset_index(drop=True)
            else:
                # Fallback: sort by the primary grouping variable
                merged_df_sorted = merged_df.sort_values(group_var).reset_index(drop=True)
            
            # Scatter plot with error bars using sorted data - color by axial ligands
            if 'axials' in merged_df_sorted.columns:
                # Color by axial ligand combinations using consistent color mapping
                unique_axials = merged_df_sorted['axials'].unique()
                axial_color_map = get_axial_ligand_colors(sorted(unique_axials), extended=True)
                
                # Plot each point with axial ligand color
                for idx, row in merged_df_sorted.iterrows():
                    axial_group = row['axials']
                    color = axial_color_map[axial_group]
                    ax1.errorbar(idx, row['E°approx (V)'], 
                                yerr=row['std_dev (V)'], fmt='o',
                                color=color, alpha=0.7, capsize=3, markersize=6)
                
                # Create legend for axial ligands using Line2D import
                from matplotlib.lines import Line2D
                legend_handles = []
                for axial in sorted(unique_axials):
                    legend_handles.append(Line2D([0], [0], marker='o', color='w', 
                                                markerfacecolor=axial_color_map[axial], 
                                                markersize=8, label=axial, alpha=0.7))
                ax1.legend(handles=legend_handles, title='Axial Ligands', 
                          fontsize=12, title_fontsize=14, loc='upper right')
            else:
                # Fallback to original grouping if axials not available
                if group_var == 'charge_multiplicity':
                    ordered_groups = sorted(merged_df_sorted[group_var].dropna().unique(), key=charge_multiplicity_sort_key)
                    scatter_color_map = colors_cm
                else:
                    ordered_groups = merged_df_sorted[group_var].dropna().unique()
                    scatter_color_map = {
                        group: scatter_colors[i % len(scatter_colors)]
                        for i, group in enumerate(ordered_groups)
                    }

                for group in ordered_groups:
                    group_data = merged_df_sorted[merged_df_sorted[group_var] == group]
                    ax1.errorbar(group_data.index, group_data['E°approx (V)'], 
                                yerr=group_data['std_dev (V)'], fmt='o',
                                color=scatter_color_map[group], 
                                label=group, alpha=0.7, capsize=3, markersize=6)
                ax1.legend(title=group_var.replace('_', ' ').title(), 
                          fontsize=12, title_fontsize=14, loc='upper right')
            
            # Highlight focus PDBs
            if focus_pdbs is not None:
                focus_data = merged_df[merged_df['structure_id'].isin(focus_pdbs)]
                ax1.scatter(focus_data.index, focus_data['E°approx (V)'],
                        s=100, facecolors='none', edgecolors='red', linewidths=2, label='Focus PDBs')
            
            ax1.set_xlabel('Structure Index', fontsize=16)
            ax1.set_ylabel('E°approx (V)', fontsize=16)
            ax1.grid(True, alpha=0.3)

            # --- Optional bottom plots ---
            if extra_plots:
                # 2. Box plots by axial ligands  
                ax3 = fig.add_subplot(gs[1, 0])
                if 'axials' in merged_df.columns:
                    axial_counts = merged_df['axials'].value_counts()
                    top_axials = axial_counts.head(6).index
                    axial_data, axial_labels = [], []
                    for axial_combo in top_axials:
                        subset = merged_df[merged_df['axials'] == axial_combo]['E°approx (V)'].dropna()
                        if len(subset) > 0:
                            axial_data.append(subset)
                            axial_labels.append(f'{axial_combo}\n(n={len(subset)})')
                    if axial_data:
                        bp2 = ax3.boxplot(axial_data, labels=axial_labels, patch_artist=True)
                        for patch, color in zip(bp2['boxes'], colors_axial[:len(axial_data)]):
                            patch.set_facecolor(color)
                            patch.set_alpha(0.7)
                        ax3.set_ylabel('E°approx (V)', fontsize=14)
                        plt.setp(ax3.get_xticklabels(), rotation=45, fontsize=12)

                # 3. High standard deviation bar plot
                ax8 = fig.add_subplot(gs[1, 1])
                high_std_pdbs = merged_df.nlargest(8, 'std_dev (V)')
                if len(high_std_pdbs) > 0:
                    pdb_ids = high_std_pdbs['structure_id'].tolist()
                    std_devs = high_std_pdbs['std_dev (V)'].tolist()
                    mean_potentials = high_std_pdbs['E°approx (V)'].tolist()
                    x_pos = range(len(pdb_ids))
                    bars = ax8.bar(x_pos, std_devs, alpha=0.7, color=get_color_sequence(len(pdb_ids), offset=1))
                    for bar, mean_pot in zip(bars, mean_potentials):
                        height = bar.get_height()
                        ax8.text(bar.get_x() + bar.get_width()/2., height + 0.01, f'{mean_pot:.2f}V', 
                                ha='center', va='bottom', fontsize=14)
                    ax8.set_xlabel('PDB ID', fontsize=14)
                    ax8.set_ylabel('Standard Deviation (V)', fontsize=14)
                    ax8.set_xticks(x_pos)
                    ax8.set_xticklabels(pdb_ids, rotation=45, fontsize=12)
                    ax8.grid(True, alpha=0.3)
            # -----------------------------

            # Set tick font sizes for all axes
            ax1.tick_params(axis='x', labelsize=12)
            ax1.tick_params(axis='y', labelsize=12)
            if extra_plots:
                if 'ax3' in locals():
                    ax3.tick_params(axis='x', labelsize=12)
                    ax3.tick_params(axis='y', labelsize=12)
                if 'ax8' in locals():
                    ax8.tick_params(axis='y', labelsize=12)

            # Save
            filename = f'{save_dir}/OERP.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            plt.show()

            if verbose: print(f"Plot saved to: {filename}")
            return merged_df, fig

        except Exception as e:
            print(f"Error in plotting reduction potentials: {e}")
            return None, None

    def create_distortion_modes_histogram(self, df, output_file="plots/porphyrin-distortions.png", decode_dicts=None):
        """
        Create histograms of porphyrin distortion modes (ruffling, saddling, doming) 
        with axial ligand color coding. Data should come from preprocessed dataframe.
        
        Parameters:
        -----------
        df : pd.DataFrame
            Preprocessed dataframe containing distortion data
        output_file : str
            Path to save the plot
        decode_dicts : dict, optional
            Dictionary for decoding encoded columns
        """
        # Preprocess the dataframe if decoding is needed
        df = self.preprocess_dataframe(df.copy(), decode_dicts)
        
        # Check if required distortion columns exist
        required_cols = ['ruffling', 'saddling', 'doming']
        if df.empty or not all(col in df.columns for col in required_cols):
            print("No distortion mode data to plot!")
            return

        # Get axial ligand data for coloring
        plot_df = df.copy()
        use_axial_colors = False
        axial_color_map = None
        axial_combinations = None
        
        # Check if we have axial combination data
        if 'axials' in plot_df.columns and plot_df['axials'].notna().sum() > 0:
            use_axial_colors = True
            # Only use the 5 most common combinations
            allowed_combos = ['CYS-HOH', 'HIS-HIS', 'HIS-HOH', 'HIS-MET', 'HIS-OXY']
            axial_combinations = [combo for combo in allowed_combos 
                                if combo in plot_df['axials'].values]
            
            # Use centralized color mapping for axial combinations
            axial_color_map = AXIAL_LIGAND_COLOR_MAP
            
            print(f"Axial ligand coloring enabled with {len(axial_combinations)} combinations:")
            for combo in axial_combinations:
                count = (plot_df['axials'] == combo).sum()
                print(f"  {combo}: {count} structures")

        # Create figure with subplots for each distortion mode
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 6))
        
        # Define distortion parameters: (column_name, x_range, title)
        distortion_params = [
            ('ruffling', (-1.0, 1.0), 'Ruffling'),
            ('saddling', (-1.0, 0.5), 'Saddling'), 
            ('doming', (-0.5, 0.5), 'Doming')
        ]
        
        axes = [ax1, ax2, ax3]
        
        # Function to create histogram for given distortion mode
        def create_single_histogram(ax, distortion_col, x_range, title):
            distortions = plot_df[distortion_col].dropna()
            
            if distortions.empty:
                ax.set_xlabel(f'{title} (Å)', fontsize=14)
                ax.set_ylabel('Frequency', fontsize=14)
                return
            
            x_min, x_max = x_range
            
            if use_axial_colors and axial_color_map is not None and axial_combinations is not None:
                # Create stacked histogram by axial combination
                combo_data = {}
                for combo in axial_combinations:
                    combo_mask = plot_df['axials'] == combo
                    combo_values = plot_df.loc[combo_mask, distortion_col].dropna()
                    combo_data[combo] = combo_values.values
                
                # Create stacked histogram - filter out empty combinations
                valid_combos = [combo for combo in axial_combinations if len(combo_data[combo]) > 0]
                data_arrays = [combo_data[combo] for combo in valid_combos]
                colors = [axial_color_map[combo] for combo in valid_combos]
                labels = [f"{combo} (n={len(combo_data[combo])})" for combo in valid_combos]
                
                # Create stacked histogram
                ax.hist(data_arrays, bins=np.linspace(x_min, x_max, 51), color=colors, alpha=0.8, 
                        edgecolor='black', label=labels, stacked=True)
                
                # Add legend only to the first subplot
                if ax == ax1:
                    ax.legend(title="Axial Ligand Combinations", fontsize=10, title_fontsize=12,
                             frameon=False, fancybox=False, shadow=False, loc='upper right')
                
            else:
                # Fall back to simple histogram
                ax.hist(distortions, bins=np.linspace(x_min, x_max, 51), alpha=0.7, color='lightcoral', edgecolor='black')
            
            # Formatting
            ax.set_xlabel(f'{title} (Å)', fontsize=14)
            ax.set_ylabel('Frequency', fontsize=14)
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlim(x_min, x_max)
            ax.grid(True, alpha=0.3)
            
            # Print statistics
            print(f"\n{title} Statistics:")
            print(f"n={len(distortions)}, mean={np.mean(distortions):.3f}, std={np.std(distortions):.3f}")
        
        # Create histograms for each distortion mode
        for i, (ax, (col, x_range, title)) in enumerate(zip(axes, distortion_params)):
            create_single_histogram(ax, col, x_range, title)

        # Calculate overall statistics for each distortion mode
        for col, _, title in distortion_params:
            if col in plot_df.columns:
                all_values = plot_df[col].dropna()
                if not all_values.empty:
                    print(f"\nOverall {title} Statistics:")
                    print(f"Min: {all_values.min():.3f} Å, Max: {all_values.max():.3f} Å")
                    print(f"Mean: {all_values.mean():.3f} ± {all_values.std():.3f} Å")
                    print(f"Total measurements: {len(all_values)}")

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()

        print(f"Distortion modes histogram saved to {output_file}")
