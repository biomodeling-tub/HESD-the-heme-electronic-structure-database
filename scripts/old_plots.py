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

    def plot_continuous_continuous(self, df, col_x, col_y,
                                   annotate_counts=False, color_intensity=False,
                                   group_counts=False, decode_dicts=None, ax=None):
        """
        Plots continuous vs. continuous variables using a scatter plot.
        If group_counts=True, groups by the 'charge-multiplicity' column.
        """
        os.makedirs('plots', exist_ok=True)

        own_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
            own_fig = True

        if decode_dicts is not None and not df.attrs.get("decoded", False):
            df = self.preprocess_dataframe(df, decode_dicts)
        df = self.filter_negative_one(df, col_x, col_y)

        if group_counts:
            if 'charge' not in df.columns or 'multiplicity' not in df.columns:
                raise ValueError("DataFrame must contain 'charge' and 'multiplicity' columns for group_counts option.")
            df = self.create_charge_multiplicity_column(df)
            groups = df.groupby('charge_multiplicity')
            for combo, group in groups:
                count = len(group)
                label = f"{combo} (n={count})"
                ax.scatter(group[col_x], group[col_y],
                           s=60, edgecolor='w', alpha=0.8, label=label)
            ax.set_xlabel(col_x, fontsize=14)
            ax.set_ylabel(col_y, fontsize=14)
            # layout legend
            n_items = len(groups)
            if self.split_legend_flag or (col_x in self.legend_split_variables or col_y in self.legend_split_variables):
                ncol_val = int(np.ceil(n_items / 2.0))
                legend_obj = ax.legend(title="Charge–Multiplicity", ncol=ncol_val)
            else:
                legend_obj = ax.legend(title="Charge–Multiplicity")
            for text in legend_obj.get_texts():
                if text.get_text().strip().lower() == "nan":
                    text.set_text("")
        elif color_intensity:
            aggregated = df.groupby([col_x, col_y]).size().reset_index(name='count')
            sc = ax.scatter(aggregated[col_x], aggregated[col_y],
                            c=aggregated['count'], cmap='viridis', s=60,
                            edgecolor='w', alpha=0.8)
            plt.colorbar(sc, ax=ax, label="Count")
            ax.set_xlabel(col_x, fontsize=14)
            ax.set_ylabel(col_y, fontsize=14)
        elif annotate_counts:
            aggregated = df.groupby([col_x, col_y]).size().reset_index(name='count')
            ax.scatter(aggregated[col_x], aggregated[col_y],
                       s=60, color='blue', edgecolor='w', alpha=0.8)
            for _, row in aggregated.iterrows():
                ax.text(row[col_x], row[col_y], str(row['count']),
                        fontsize=12, color='black',
                        verticalalignment='bottom', horizontalalignment='right')
            ax.set_xlabel(col_x, fontsize=14)
            ax.set_ylabel(col_y, fontsize=14)
        else:
            sns.scatterplot(data=df, x=col_x, y=col_y, s=60,
                            edgecolor='w', alpha=0.8, ax=ax)
            ax.set_xlabel(col_x, fontsize=14)
            ax.set_ylabel(col_y, fontsize=14)

        if own_fig:
            filename = f"plots/scatter_{col_x}_{col_y}.png"
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)

    def plot_continuous_categorical(self, df, col_x, col_y,
                                    decode_dicts=None, absolute_counts=False,
                                    ax=None, axials=False, log_y=False):
        """
        Plots one continuous and one categorical variable using a KDE plot.
        Saves to plots/kde_{col_x}_{col_y}.png
        """
        os.makedirs('plots', exist_ok=True)

        own_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
            own_fig = True

        if decode_dicts is not None and not df.attrs.get("decoded", False):
            df = self.preprocess_dataframe(df, decode_dicts)

        if not axials:
            def is_categorical(series):
                if pd.api.types.is_numeric_dtype(series):
                    if pd.api.types.is_integer_dtype(series) or (
                       pd.api.types.is_float_dtype(series) and (series.dropna() % 1 == 0).all()):
                        return series.nunique() <= 10
                    return False
                return True
            if is_categorical(df[col_x]):
                continuous_var = col_y
                categorical_var = col_x
            else:
                continuous_var = col_x
                categorical_var = col_y
        else:
            continuous_var = col_y

        if axials:
            df_clean = df[df[continuous_var] != -1]
        else:
            df_clean = self.filter_negative_one(df, continuous_var, categorical_var)

        if axials:
            all_axials = pd.unique(df_clean[['axial1', 'axial2']].values.ravel())
            groups = sorted([g for g in all_axials if pd.notna(g) and g != ""])
            plotted_any = False
            for grp in groups:
                subset = df_clean[(df_clean['axial1'] == grp) | (df_clean['axial2'] == grp)]
                if subset[continuous_var].nunique() <= 1:
                    print(f"Skipping KDE for axial '{grp}' due to zero variance.")
                    continue
                if absolute_counts:
                    data = subset[continuous_var].dropna()
                    kde = gaussian_kde(data)
                    x_min, x_max = data.min(), data.max()
                    x_vals = np.linspace(x_min, x_max, 200)
                    density = kde(x_vals) * len(data)
                    ax.plot(x_vals, density, label=str(grp))
                else:
                    sns.kdeplot(subset[continuous_var], label=str(grp), warn_singular=False, ax=ax)
                plotted_any = True
            ax.set_xlabel(continuous_var, fontsize=14)
            ax.set_ylabel("Density", fontsize=14)
            if plotted_any:
                ncol_val = int(np.ceil(len(groups) / 2.0)) if self.split_legend_flag or axials else len(groups)
                legend_obj = ax.legend(title="Axial", ncol=ncol_val)
                for text in legend_obj.get_texts():
                    if text.get_text().strip().lower() == "nan":
                        text.set_text("")
        else:
            plotted_any = False
            for cat in sorted(df_clean[categorical_var].dropna().unique()):
                subset = df_clean[df_clean[categorical_var] == cat]
                if subset[continuous_var].nunique() <= 1:
                    print(f"Skipping KDE for category '{cat}' in '{continuous_var}' due to zero variance.")
                    continue
                if absolute_counts:
                    data = subset[continuous_var].dropna()
                    kde = gaussian_kde(data)
                    x_min, x_max = data.min(), data.max()
                    x_vals = np.linspace(x_min, x_max, 200)
                    density = kde(x_vals) * len(data)
                    ax.plot(x_vals, density, label=str(cat))
                else:
                    sns.kdeplot(subset[continuous_var], label=str(cat), warn_singular=False, ax=ax)
                plotted_any = True
            ax.set_xlabel(continuous_var, fontsize=14)
            ax.set_ylabel("Density", fontsize=14)
            if plotted_any:
                unique_vals = sorted(df_clean[categorical_var].dropna().unique())
                ncol_val = int(np.ceil(len(unique_vals) / 2.0)) if (self.split_legend_flag or (categorical_var in self.legend_split_variables)) else len(unique_vals)
                legend_obj = ax.legend(title=categorical_var, ncol=ncol_val)
                for text in legend_obj.get_texts():
                    if text.get_text().strip().lower() == "nan":
                        text.set_text("")

        if log_y:
            ax.set_yscale('log')

        if own_fig:
            filename = f"plots/kde_{col_x}_{col_y}.png"
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)

    def plot_continuous_categorical_histogram(self, df, col_x, col_y,
                                              decode_dicts=None, ax=None, bin_width=None):
        """
        Plots one continuous and one categorical variable using per-category histograms.
        Saves to plots/hist_{col_x}_{col_y}.png

        If bin_width is None, it will auto-compute using Freedman-Diaconis rule per category.
        """
        os.makedirs('plots', exist_ok=True)
        own_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
            own_fig = True
        if decode_dicts is not None and not df.attrs.get("decoded", False):
            df = self.preprocess_dataframe(df, decode_dicts)
        # Determine which is continuous vs categorical
        def is_categorical(series):
            if pd.api.types.is_numeric_dtype(series):
                if pd.api.types.is_integer_dtype(series) or (
                   pd.api.types.is_float_dtype(series) and (series.dropna() % 1 == 0).all()):
                    return series.nunique() <= 10
                return False
            return True
        if is_categorical(df[col_x]):
            continuous_var = col_y
            categorical_var = col_x
        else:
            continuous_var = col_x
            categorical_var = col_y
        df_clean = self.filter_negative_one(df, continuous_var, categorical_var)
        unique_categories = sorted(df_clean[categorical_var].dropna().unique())
        plotted_any = False
        # Compute global bin edges if bin_width is given; otherwise, each category will auto-bin
        if bin_width is not None:
            global_min = df_clean[continuous_var].min()
            global_max = df_clean[continuous_var].max()
            bins = np.arange(global_min, global_max + bin_width, bin_width)
        else:
            bins = None  # let matplotlib auto-bin per call
        for cat in unique_categories:
            subset = df_clean[df_clean[categorical_var] == cat]
            if subset[continuous_var].nunique() == 0:
                continue

            data = subset[continuous_var].dropna()
            if data.empty:
                continue

            # If bin_width is None, call hist without explicit bins
            if bins is not None:
                ax.hist(data, bins=bins, alpha=0.5, edgecolor='black', label=str(cat))
            else:
                ax.hist(data, alpha=0.5, edgecolor='black', label=str(cat))
            plotted_any = True

        ax.set_xlabel(continuous_var, fontsize=14)
        ax.set_ylabel("Count", fontsize=14)
        if plotted_any:
            ncol_val = int(np.ceil(len(unique_categories) / 2.0)) if self.split_legend_flag or (categorical_var in self.legend_split_variables) else len(unique_categories)
            legend_obj = ax.legend(title=categorical_var, ncol=ncol_val)
            for text in legend_obj.get_texts():
                if text.get_text().strip().lower() == "nan":
                    text.set_text("")

        if own_fig:
            filename = f"plots/hist_{col_x}_{col_y}.png"
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)

    def plot_categorical_categorical(self, df, col_x, col_y, decode_dicts=None, ax=None):
        """
        Plots two categorical variables using a heatmap of counts.
        Saves to plots/heatmap_{col_x}_{col_y}.png
        """
        os.makedirs('plots', exist_ok=True)

        own_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
            own_fig = True

        if decode_dicts is not None and not df.attrs.get("decoded", False):
            df = self.preprocess_dataframe(df, decode_dicts)
        df = self.filter_negative_one(df, col_x, col_y)

        ct = pd.crosstab(df[col_x], df[col_y])
        if ct.empty:
            print(f"Warning: Crosstab for {col_x} vs {col_y} is empty. Skipping plot.")
        else:
            sns.heatmap(ct, annot=True, fmt='d', cmap="YlGnBu", ax=ax)
            ax.set_xlabel(col_y, fontsize=14)
            ax.set_ylabel(col_x, fontsize=14)

        if own_fig and not ct.empty:
            filename = f"plots/heatmap_{col_x}_{col_y}.png"
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)

    def infer_col_type(self, series):
        """
        Helper method to infer whether a column is 'continuous' or 'categorical'.
        """
        if pd.api.types.is_numeric_dtype(series):
            if (pd.api.types.is_integer_dtype(series) or
                (pd.api.types.is_float_dtype(series) and (series.dropna() % 1 == 0).all())) and series.nunique() <= 30:
                return 'categorical'
            else:
                return 'continuous'
        else:
            return 'categorical'

    def plot_all_scatter_combinations(
        self,
        df,
        annotate_counts=False,
        color_intensity=False,
        group_counts=False,
        numeric_only=False,
        decode_dicts=None,
        selected_columns=None,
        cc_plot="kde"   # NEW: choose "kde" or "hist" for continuous‐categorical combos
        ):
        """
        Iterates over all unique pairs of columns in the DataFrame, sorts them into three categories:
        - continuous-continuous,
        - continuous-categorical, and
        - categorical-categorical,
        and calls the appropriate plotting function for each combination.

        New parameter:
        cc_plot : str, either "kde" (default) or "hist"
            Determines whether continuous‐categorical pairs use KDE (plot_continuous_categorical)
            or per‐category histogram (plot_continuous_categorical_histogram).
        """
        df = self.preprocess_dataframe(df, decode_dicts)
        decode_dicts = None

        if selected_columns is not None:
            cols = [col for col in selected_columns if col in df.columns]
        else:
            if numeric_only:
                cols = df.select_dtypes(include='number').columns.tolist()
            else:
                cols = list(df.columns)

        col_types = {}
        for col in cols:
            series = df[col]
            if pd.api.types.is_numeric_dtype(series):
                if (pd.api.types.is_integer_dtype(series) or
                    (pd.api.types.is_float_dtype(series) and (series.dropna() % 1 == 0).all())) and series.nunique() <= 30:
                    col_types[col] = 'categorical'
                else:
                    col_types[col] = 'continuous'
            else:
                col_types[col] = 'categorical'

        for col_x, col_y in itertools.combinations(cols, 2):
            print(f"Plotting {col_x} vs {col_y} with types {col_types[col_x]} and {col_types[col_y]}")
            if col_types[col_x] == 'continuous' and col_types[col_y] == 'continuous':
                self.plot_continuous_continuous(
                    df, col_x, col_y,
                    annotate_counts, color_intensity,
                    group_counts, decode_dicts
                )

            elif {col_types[col_x], col_types[col_y]} == {'continuous', 'categorical'}:
                # If user asked for histograms on continuous‐categorical:
                if cc_plot.lower() == "hist":
                    self.plot_continuous_categorical_histogram(
                        df, col_x, col_y, decode_dicts
                    )
                else:
                    # default to KDE
                    self.plot_continuous_categorical(
                        df, col_x, col_y, decode_dicts
                    )

            else:
                self.plot_categorical_categorical(df, col_x, col_y, decode_dicts)

    def get_plot_func(self, plot_type, df, col_x, col_y, **kwargs):
        """
        Returns a function that generates the specified plot on an Axes,
        and optionally overrides axis labels via kwargs:
         - xlabel: custom x-axis label
         - ylabel: custom y-axis label
        Remaining kwargs passed to underlying plotting routines.
        """
        # Determine plot type automatically if needed
        if plot_type is None or plot_type == "auto":
            type_x = self.infer_col_type(df[col_x])
            type_y = self.infer_col_type(df[col_y])
            if type_x == "continuous" and type_y == "continuous":
                plot_type = "continuous_continuous"
            elif (type_x == "continuous" and type_y == "categorical") or \
                 (type_x == "categorical" and type_y == "continuous"):
                plot_type = "continuous_categorical"
            else:
                plot_type = "categorical_categorical"
        # Common label overrides
        xlabel = kwargs.pop('xlabel', None)
        ylabel = kwargs.pop('ylabel', None)

        if plot_type == "continuous_continuous":
            def plot_func(ax):
                self.plot_continuous_continuous(
                    df, col_x, col_y,
                    annotate_counts=kwargs.get("annotate_counts", False),
                    color_intensity=kwargs.get("color_intensity", False),
                    group_counts=kwargs.get("group_counts", False),
                    decode_dicts=kwargs.get("decode_dicts", None),
                    ax=ax
                )
                if xlabel:
                    ax.set_xlabel(xlabel)
                if ylabel:
                    ax.set_ylabel(ylabel)
            return plot_func

        elif plot_type == "continuous_categorical":
            def plot_func(ax):
                # Check whether caller explicitly passed 'cc_plot' in kwargs:
                if kwargs.get("cc_plot", "kde").lower() == "hist":
                    self.plot_continuous_categorical_histogram(
                        df, col_x, col_y,
                        decode_dicts=kwargs.get("decode_dicts", None),
                        ax=ax,
                        bin_width=kwargs.get("bin_width", None)
                    )
                else:
                    self.plot_continuous_categorical(
                        df, col_x, col_y,
                        decode_dicts=kwargs.get("decode_dicts", None),
                        absolute_counts=kwargs.get("absolute_counts", False),
                        ax=ax,
                        axials=kwargs.get("axials", False),
                        log_y=kwargs.get("log_y", False)
                    )
                if xlabel:
                    ax.set_xlabel(xlabel)
                if ylabel:
                    ax.set_ylabel(ylabel)
            return plot_func

        elif plot_type == "categorical_categorical":
            def plot_func(ax):
                self.plot_categorical_categorical(
                    df, col_x, col_y,
                    decode_dicts=kwargs.get("decode_dicts", None),
                    ax=ax
                )
                if xlabel:
                    ax.set_xlabel(xlabel)
                if ylabel:
                    ax.set_ylabel(ylabel)
            return plot_func

        else:
            raise ValueError(f"Unknown plot_type: {plot_type}")

    def create_multiplot(self, plot_funcs, grid_rows=3, grid_cols=3, filename='multiplot.png'):
        """
        Arrange multiple plots (each provided as a function that accepts an Axes object)
        in a grid layout and save to a PNG file.
        Supports an instance-level flag `force_single_row_legend`:
          - If True, legend always appears in one row.
        """
        # set up figure and axes
        fig, axes = plt.subplots(grid_rows, grid_cols,
                                 figsize=(grid_cols * 5, grid_rows * 4))
        # normalize to 1D array even if 1x1
        axes = np.atleast_1d(axes).flatten()
        # draw each subplot
        for i, plot_func in enumerate(plot_funcs):
            if i < len(axes):
                plot_func(axes[i])
            else:
                print("Warning: More plot functions provided than available subplots.")
                break
        # turn off unused axes
        for j in range(len(plot_funcs), len(axes)):
            axes[j].axis('off')
        # collect a single global legend
        global_handles, global_labels = None, None
        for ax in axes:
            legend = ax.get_legend()
            if legend is not None:
                handles, labels = ax.get_legend_handles_labels()
                if handles and global_handles is None:
                    global_handles, global_labels = handles, labels
                ax.get_legend().remove()
        if global_handles is not None:
            # clean up 'nan' labels
            global_labels = ["" if str(l).strip().lower() == "nan" else l for l in global_labels]
            n_items = len(global_labels)
            # determine columns for legend
            if self.force_single_row_legend:
                ncol_val = n_items
            else:
                split_flag = self.split_legend_flag or ('axials' in self.legend_split_variables)
                ncol_val = int(np.ceil(n_items / 2.0)) if split_flag else n_items
            # place the legend
            global_legend = fig.legend(
                global_handles, global_labels,
                loc='upper center',
                bbox_to_anchor=(0.5, 1.05),
                ncol=ncol_val,
                title="Groups",
                borderaxespad=0.1
            )
            # remove any empty labels
            for text in global_legend.get_texts():
                if text.get_text().strip() == "":
                    text.set_text("")
            if global_legend.get_title():
                global_legend.get_title().set_fontsize(9)
        # style axes
        for ax in axes:
            ax.title.set_fontsize(16)
            ax.xaxis.label.set_fontsize(14)
            ax.yaxis.label.set_fontsize(14)
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(12)
        plt.tight_layout(pad=3)
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)

    def create_multiplot_from_list(self, plot_specs, grid_rows=3, grid_cols=3, filename='multiplot.png'):
        """
        Takes a list of specs [plot_type, df, col_x, col_y, extra_kwargs]
        or [df, col_x, col_y, extra_kwargs]; auto-detects type if missing.
        Arranges the plots in a grid and saves to filename.
        """
        plot_funcs = []
        for spec in plot_specs:
            # ensure plot_type present
            if not (isinstance(spec[0], str) and spec[0] in ('auto', 'continuous_continuous',
                                                           'continuous_categorical',
                                                           'categorical_categorical')):
                spec = ['auto'] + spec
            # unpack
            if len(spec) == 5:
                plot_type, df, col_x, col_y, extra_kwargs = spec
            else:
                plot_type, df, col_x, col_y = spec
                extra_kwargs = {}
            func = self.get_plot_func(plot_type, df, col_x, col_y, **extra_kwargs)
            plot_funcs.append(func)
        # delegate to the updated create_multiplot
        self.create_multiplot(plot_funcs, grid_rows, grid_cols, filename)

    def scatter_plot(self, df, col_x, col_y,
                     annotate_counts=False, color_intensity=False,
                     group_counts=False, decode_dicts=None):
        """
        Plots a scatter plot of two columns from a DataFrame.
        If group_counts=True, groups by the 'charge' and 'multiplicity' columns.
        """
        plt.figure(figsize=(10, 6))
        if decode_dicts is not None and isinstance(decode_dicts, dict):
            # decoding logic unchanged
            for mapping_dict in [{k: v} for k, v in decode_dicts.items()]:
                for key, mapping in mapping_dict.items():
                    inverted = {v: k for k, v in mapping.items()}
                    if key == col_x:
                        df[col_x] = df[col_x].map(inverted)
                    if key == col_y:
                        df[col_y] = df[col_y].map(inverted)
        if group_counts:
            if 'charge' not in df.columns or 'multiplicity' not in df.columns:
                raise ValueError("DataFrame must contain 'charge' and 'multiplicity' columns for group_counts option.")
            df = self.create_charge_multiplicity_column(df)
            groups = df.groupby('charge_multiplicity')
            for combo, group in groups:
                count = len(group)
                label = f"{combo} (n={count})"
                plt.scatter(group[col_x], group[col_y],
                            s=60, edgecolor='w', alpha=0.8, label=label)
            plt.xlabel(col_x, fontsize=14)
            plt.ylabel(col_y, fontsize=14)
            plt.title(f"{col_y} vs {col_x} (by charge-multiplicity)", fontsize=16)
            plt.legend(title="Charge–Multiplicity", fontsize=12)
        elif color_intensity:
            aggregated = df.groupby([col_x, col_y]).size().reset_index(name='count')
            sc = plt.scatter(aggregated[col_x], aggregated[col_y],
                             c=aggregated['count'], cmap='viridis', s=60,
                             edgecolor='w', alpha=0.8)
            plt.colorbar(sc, label="Count")
            plt.xlabel(col_x, fontsize=14)
            plt.ylabel(col_y, fontsize=14)
        elif annotate_counts:
            aggregated = df.groupby([col_x, col_y]).size().reset_index(name='count')
            plt.scatter(aggregated[col_x], aggregated[col_y],
                        s=60, color='blue', edgecolor='w', alpha=0.8)
            for _, row in aggregated.iterrows():
                plt.text(row[col_x], row[col_y], str(row['count']),
                         fontsize=12, color='black',
                         verticalalignment='bottom', horizontalalignment='right')
            plt.xlabel(col_x, fontsize=14)
            plt.ylabel(col_y, fontsize=14)
        else:
            ax = sns.scatterplot(data=df, x=col_x, y=col_y, s=60,
                                 edgecolor='w', alpha=0.8)
            ax.set_xlabel(col_x, fontsize=14)
            ax.set_ylabel(col_y, fontsize=14)
        plt.tight_layout()
        plt.savefig(f"plots/scatter_{col_x}_{col_y}.png")
        plt.close()

    def plot_value_pair_histogram(self, df, col_x="sum_mulliken_charges.charge_sum", col_y="sum_mulliken_charges.spin_sum"):
        """
        Plots a 2D hexbin histogram of the distribution of value pairs from the specified columns.
        """
        plt.figure(figsize=(8, 6))
        g = sns.jointplot(data=df, x=col_x, y=col_y, kind="hex", height=8, color="skyblue")
        g.set_axis_labels(col_x, col_y)
        #plt.show()

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
        import re
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
                import matplotlib.pyplot as plt
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
                import seaborn as sns
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
            ax = sns.lineplot(
                data=df_melt,
                x="orbital", y="energy",
                hue=hue,
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

    def plot_histograms_by_charge_multiplicity(
        self,
        df,
        col,
        charge_col="charge_multiplicity",
        bin_width=0.1,
        decode_dicts=None,
        figsize=(10, 8),
        xlabel="ΔE (HOMO–LUMO) [a.u.]",
        ylabel="Count",
        tight_layout=True,
        savepath=None,
        axials_as_color=False,
        group_by="axial"   # new parameter: "axial" (default) or "charge"
    ):
        """
        Two overlay modes and two grouping orders:

        1) If axials_as_color=True:
           • Single‐axes overlay of one histogram per axial ligand.  Ignores grouping order.

        Otherwise (axials_as_color=False), there are two grouping orders:
        - group_by="axial" (default):
            • For each axial ligand, create one figure.  In that figure, make a grid of
              subplots—one per charge-multiplicity—and plot the histogram of `col` for
              rows matching (axial, charge-multiplicity).

        - group_by="charge":
            • For each charge-multiplicity, create one figure.  In that figure, make a
              grid of subplots—one per axial ligand—and plot the histogram of `col` for
              rows matching (charge-multiplicity, axial).

        Parameters
        ----------
        df : pandas.DataFrame
            Must contain `col`, `charge`, `multiplicity`, `axial1`, and `axial2`.  `col` must be numeric.
        col : str
            The continuous column to histogram (e.g. "delta_HOMO-LUMO").
        charge_col : str
            Name of the charge–multiplicity column (default "charge_multiplicity").
        bin_width : float
            Width of each bin (global bins = np.arange(min, max+bin_width, bin_width)).
        decode_dicts : dict or None
            If provided, keys = column names to decode via mapping (invert mapping, apply).
        figsize : tuple
            Size (width, height) in inches of each figure.
        xlabel, ylabel : str
            Axis labels.
        tight_layout : bool
            If True, call plt.tight_layout(pad=2) on each figure before showing.
        savepath : str or None
            If not None, figures will be saved.  For group_by="axial", saved as
            "{savepath}_axial_{safe_axial}.png".  For group_by="charge", saved as
            "{savepath}_charge_{safe_charge}.png".  For axials_as_color=True, saved
            exactly to `savepath`.
        axials_as_color : bool
            If True, ignore charge-multiplicity entirely and overlay histograms of `col`
            for each axial ligand in one axes (single‐plot mode).
        group_by : str, either "axial" (default) or "charge"
            If axials_as_color=False, determines whether to loop first over axial ligands
            ("axial") or first over charge-multiplicity values ("charge") when making figures.
        """
        # 1) Decode if needed
        if decode_dicts is not None and not df.attrs.get("decoded", False):
            for c, mapping in decode_dicts.items():
                if c not in df.columns:
                    continue
                inv_map = {v: k for k, v in mapping.items()}
                df[c] = df[c].map(lambda x: inv_map.get(x, x))
            df.attrs["decoded"] = True

        # 2) Drop invalid entries in target column
        df_clean = df[df[col] != -1]

        # 3) Ensure axial columns exist
        if ("axial1" not in df_clean.columns) or ("axial2" not in df_clean.columns):
            raise ValueError("DataFrame must contain 'axial1' and 'axial2' for axial grouping.")

        # 4) Gather all axial ligands
        all_axials = np.unique(df_clean[["axial1", "axial2"]].values.ravel())
        all_axials = [a for a in sorted(all_axials) if pd.notna(a) and a != ""]
        if len(all_axials) == 0:
            raise ValueError("No valid axial-ligand values found in 'axial1'/'axial2'.")

        # 5) Compute global bin edges across the entire df_clean
        global_min = df_clean[col].min()
        global_max = df_clean[col].max()
        bins = np.arange(global_min, global_max + bin_width, bin_width)
        if len(bins) < 2:
            raise ValueError("bin_width is too large or data is too narrow; no bins created.")

        ##########
        # Mode 1: overlay by axial (ignore charge-multiplicity)
        if axials_as_color:
            fig, ax = plt.subplots(figsize=figsize)
            for axial in all_axials:
                mask = (df_clean["axial1"] == axial) | (df_clean["axial2"] == axial)
                sub_axial = df_clean[mask]
                if sub_axial.empty:
                    continue

                if sub_axial[col].nunique() == 1:
                    singleton = sub_axial[col].iloc[0]
                    ax.bar(
                        [singleton],
                        [len(sub_axial)],
                        width=bin_width,
                        alpha=0.5,
                        edgecolor='k',
                        label=str(axial)
                    )
                else:
                    ax.hist(
                        sub_axial[col].dropna(),
                        bins=bins,
                        alpha=0.5,
                        edgecolor='black',
                        label=str(axial)
                    )

            ax.set_xlabel(xlabel, fontsize=14)
            ax.set_ylabel(ylabel, fontsize=14)
            ax.legend(title="Axial", fontsize=12, title_fontsize=14)

            if tight_layout:
                plt.tight_layout(pad=2)
            if savepath:
                fig.savefig(savepath, dpi=300, bbox_inches="tight")

            plt.show()
            plt.close(fig)
            return

        ##########
        # Modes 2 & 3: grid of pairwise (axial, charge) subplots
        if 'charge' not in df_clean.columns or 'multiplicity' not in df_clean.columns:
            raise ValueError("DataFrame must contain 'charge' and 'multiplicity' columns for grouping.")
        
        df_clean = self.create_charge_multiplicity_column(df_clean)
        unique_groups = sorted(df_clean[charge_col].dropna().unique())
        n_groups = len(unique_groups)
        if n_groups == 0:
            raise ValueError(f"No non-null values found in '{charge_col}'.")

        # Compute a “near-square” grid for n subplots
        def _grid_dimensions(n):
            side = int(np.ceil(np.sqrt(n)))
            return side, side

        if group_by.lower() == "axial":
            # Mode 2: first loop over axial → figures per axial, subplots per charge
            nrows, ncols = _grid_dimensions(n_groups)

            for axial in all_axials:
                mask_ax = (df_clean["axial1"] == axial) | (df_clean["axial2"] == axial)
                df_axial = df_clean[mask_ax]
                if df_axial.empty:
                    continue

                fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
                axes_flat = np.array(axes).reshape(-1)

                for idx, group_val in enumerate(unique_groups):
                    if idx >= len(axes_flat):
                        break
                    ax = axes_flat[idx]

                    df_sub = df_axial[df_axial[charge_col] == group_val]
                    charge_str, mult_str = str(group_val).split('_', 1)

                    if df_sub.empty:
                        ax.axis("off")
                        continue

                    if df_sub[col].nunique() == 1:
                        singleton = df_sub[col].iloc[0]
                        ax.bar(
                            [singleton],
                            [len(df_sub)],
                            width=bin_width,
                            alpha=0.7,
                            edgecolor='k'
                        )
                    else:
                        ax.hist(
                            df_sub[col].dropna(),
                            bins=bins,
                            alpha=0.7,
                            edgecolor='black'
                        )

                    total_n = len(df_sub)
                    ax.set_xlim(global_min, global_max)
                    if idx % ncols != 0:
                        ax.set_yticklabels([])

                # Turn off any extra subplots beyond n_groups
                for j in range(n_groups, nrows * ncols):
                    axes_flat[j].axis("off")

                # Shared x/y labels
                for ax in axes_flat[(ncols * (nrows - 1)):]:
                    ax.set_xlabel(xlabel)
                for ax in axes_flat[0::ncols]:
                    ax.set_ylabel(ylabel)

                # Title the figure with the axial ligand name

                if tight_layout:
                    plt.tight_layout(pad=2)
                if savepath:
                    safe_axial = re.sub(r"\W+", "_", str(axial))
                    fig.savefig(f"{savepath}_axial_{safe_axial}.png", dpi=300, bbox_inches="tight")

                plt.show()
                plt.close(fig)

        elif group_by.lower() == "charge":
            # Mode 3: first loop over charge → figures per charge, subplots per axial
            nrows, ncols = _grid_dimensions(len(all_axials))

            for group_val in unique_groups:
                df_charge = df_clean[df_clean[charge_col] == group_val]
                if df_charge.empty:
                    continue

                fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
                axes_flat = np.array(axes).reshape(-1)

                for idx, axial in enumerate(all_axials):
                    if idx >= len(axes_flat):
                        break
                    ax = axes_flat[idx]

                    # Select rows with this charge AND this axial
                    mask = (df_charge["axial1"] == axial) | (df_charge["axial2"] == axial)
                    df_sub = df_charge[mask]

                    charge_str, mult_str = str(group_val).split('_', 1)

                    if df_sub.empty:
                        ax.axis("off")
                        continue

                    if df_sub[col].nunique() == 1:
                        singleton = df_sub[col].iloc[0]
                        ax.bar(
                            [singleton],
                            [len(df_sub)],
                            width=bin_width,
                            alpha=0.7,
                            edgecolor='k'
                        )
                    else:
                        ax.hist(
                            df_sub[col].dropna(),
                            bins=bins,
                            alpha=0.7,
                            edgecolor='black'
                        )

                    total_n = len(df_sub)
                    ax.set_xlim(global_min, global_max)
                    if idx % ncols != 0:
                        ax.set_yticklabels([])

                # Turn off any extra subplots beyond len(all_axials)
                for j in range(len(all_axials), nrows * ncols):
                    axes_flat[j].axis("off")

                # Shared x/y labels
                for ax in axes_flat[(ncols * (nrows - 1)):]:
                    ax.set_xlabel(xlabel)
                for ax in axes_flat[0::ncols]:
                    ax.set_ylabel(ylabel)

                # Title the figure with the charge-multiplicity value
                charge_str, mult_str = str(group_val).split('_', 1)

                if tight_layout:
                    plt.tight_layout(pad=2)
                if savepath:
                    safe_charge = re.sub(r"\W+", "_", str(group_val))
                    fig.savefig(f"{savepath}_charge_{safe_charge}.png", dpi=300, bbox_inches="tight")

                plt.show()
                plt.close(fig)

        else:
            raise ValueError("group_by must be either 'axial' or 'charge'.")

    def plot_publication_histogram(self, df, col, 
                                   group_by_axial=False, 
                                   group_by_charge=False,
                                   decode_dicts=None,
                                   bins=30,
                                   alpha=0.7,
                                   figsize=(10, 8),
                                   save=False,
                                   save_dir='plots',
                                   palette=None,
                                   edgecolor='black',
                                   linewidth=1.0,
                                   replacements=None):
        """
        Create a publication-ready histogram of a column with optional grouping by axial ligands 
        and/or charge-multiplicity combinations.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing the data to plot
        col : str
            Column name to create histogram for
        group_by_axial : bool, default False
            If True, group by combinations of axial1 and axial2 columns
        group_by_charge : bool, default False
            If True, group by combinations of charge and multiplicity columns
        decode_dicts : dict, optional
            Dictionary for decoding categorical variables
        bins : int or sequence, default 30
            Number of bins or bin edges for histogram
        alpha : float, default 0.7
            Transparency of histogram bars
        figsize : tuple, default (10, 8)
            Figure size (width, height) in inches
        save : bool, default False
            If True, save the plot to file
        save_dir : str, default 'plots'
            Directory to save the plot
        palette : list, optional
            Color palette for grouping. If None, uses colorblind-friendly defaults
        edgecolor : str, default 'black'
            Color of histogram bar edges
        linewidth : float, default 1.0
            Width of histogram bar edges
        replacements : dict, optional
            Dictionary for decoding categorical variables in legend labels.
            Format: {column_name: {encoded_value: decoded_string}}
            This is typically the encoding_dict returned by DataPreprocessor.process()

        Returns
        -------
        None
        """
        # Set publication-ready styling
        plt.rcParams.update({
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'font.family': 'serif',
            'font.serif': ['DejaVu Serif'],
            'font.size': 14,
            'axes.titlesize': 18,
            'axes.labelsize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'lines.linewidth': 1.5,
            'axes.linewidth': 1.0,
        })

        # Colorblind-friendly palette
        if palette is None:
            palette = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00']

        # Decode categorical variables if needed
        if decode_dicts is not None and not df.attrs.get("decoded", False):
            df = self.preprocess_dataframe(df, decode_dicts)

        # Filter out invalid values
        df_clean = df[df[col] != -1].copy()

        # Create directory if needed
        os.makedirs(save_dir, exist_ok=True)

        # Helper function to decode legend labels
        def decode_legend_label(label_parts, replacements):
            """
            Decode categorical values in legend labels using the replacements dictionary.
            
            Parameters
            ----------
            label_parts : list
                List of parts that make up the label (e.g., ['axial1_value', 'axial2_value'])
            replacements : dict
                Dictionary for decoding categorical variables
                
            Returns
            -------
            str
                Decoded label string
            """
            if replacements is None:
                return '-'.join(str(part) for part in label_parts)
            
            decoded_parts = []
            for i, part in enumerate(label_parts):
                # Try to decode axial1, axial2, etc. based on position
                if i == 0 and 'axial1' in replacements:
                    try:
                        decoded_parts.append(replacements['axial1'].get(int(part), str(part)))
                    except (ValueError, TypeError):
                        decoded_parts.append(str(part))
                elif i == 1 and 'axial2' in replacements:
                    try:
                        decoded_parts.append(replacements['axial2'].get(int(part), str(part)))
                    except (ValueError, TypeError):
                        decoded_parts.append(str(part))
                else:
                    decoded_parts.append(str(part))
            
            return '-'.join(decoded_parts)

        fig, ax = plt.subplots(figsize=figsize)

        # Handle different grouping scenarios
        if group_by_axial and group_by_charge:
            # Group by both axial and charge combinations
            if not all(c in df_clean.columns for c in ['axial1', 'axial2', 'charge', 'multiplicity']):
                raise ValueError("DataFrame must contain 'axial1', 'axial2', 'charge', and 'multiplicity' columns")
            
            df_clean['axial_combo'] = df_clean['axial1'].astype(str) + '-' + df_clean['axial2'].astype(str)
            df_clean['charge_combo'] = df_clean['charge'].astype(str) + '-' + df_clean['multiplicity'].astype(str)
            df_clean['combined_group'] = df_clean['axial_combo'] + '_' + df_clean['charge_combo']
            
            groups = df_clean.groupby('combined_group')
            title_suffix = "(grouped by axial and charge-multiplicity)"
            
        elif group_by_axial:
            # Group by axial combinations only
            if not all(c in df_clean.columns for c in ['axial1', 'axial2']):
                raise ValueError("DataFrame must contain 'axial1' and 'axial2' columns")
            
            df_clean['axial_combo'] = df_clean['axial1'].astype(str) + '-' + df_clean['axial2'].astype(str)
            groups = df_clean.groupby('axial_combo')
            title_suffix = "(grouped by axial ligands)"
            
        elif group_by_charge:
            # Group by charge-multiplicity combinations only
            if not all(c in df_clean.columns for c in ['charge', 'multiplicity']):
                raise ValueError("DataFrame must contain 'charge' and 'multiplicity' columns")
            
            df_clean['charge_combo'] = df_clean['charge'].astype(str) + '-' + df_clean['multiplicity'].astype(str)
            groups = df_clean.groupby('charge_combo')
            title_suffix = "(grouped by charge-multiplicity)"
            
        else:
            # No grouping - single histogram
            ax.hist(df_clean[col], bins=bins, alpha=alpha, edgecolor=edgecolor, 
                   linewidth=linewidth, color=palette[0])
            ax.set_xlabel(col, fontsize=16)
            ax.set_ylabel('Count', fontsize=16)
            title_suffix = ""

        # Plot grouped histograms if grouping is enabled
        if group_by_axial or group_by_charge:
            legend_handles = []
            
            for i, (group_name, group_data) in enumerate(groups):
                color = palette[i % len(palette)]
                
                # Decode the group name for legend display
                if group_by_axial and not group_by_charge:
                    # For axial grouping: group_name is 'axial1_value-axial2_value'
                    parts = str(group_name).split('-')
                    decoded_label = decode_legend_label(parts, replacements)
                elif group_by_charge and not group_by_axial:
                    # For charge grouping: group_name is 'charge_value-multiplicity_value'
                    decoded_label = str(group_name)  # Keep charge-multiplicity as is
                else:
                    # For combined grouping: group_name is 'axial1-axial2_charge-multiplicity'
                    if '_' in str(group_name):
                        axial_part, charge_part = str(group_name).split('_', 1)
                        axial_parts = axial_part.split('-')
                        decoded_axial = decode_legend_label(axial_parts, replacements)
                        decoded_label = f"{decoded_axial}_{charge_part}"
                    else:
                        decoded_label = str(group_name)
                
                # Plot histogram for this group
                ax.hist(group_data[col], bins=bins, alpha=alpha, 
                       edgecolor=edgecolor, linewidth=linewidth,
                       color=color, label=f'{decoded_label} (n={len(group_data)})')
                
                # Create legend handle with decoded label
                legend_handles.append(mpatches.Patch(color=color, label=f'{decoded_label} (n={len(group_data)})'))
            
            # Add legend
            ax.legend(handles=legend_handles, loc='upper right', frameon=False, fontsize=12)
            
            ax.set_xlabel(col, fontsize=16)
            ax.set_ylabel('Count', fontsize=16)

        # Improve layout
        plt.tight_layout()

        # Save if requested
        if save:
            if group_by_axial and group_by_charge:
                filename = f"histogram_{col}_axial_charge.png"
            elif group_by_axial:
                filename = f"histogram_{col}_axial.png"
            elif group_by_charge:
                filename = f"histogram_{col}_charge.png"
            else:
                filename = f"histogram_{col}.png"
            
            filepath = os.path.join(save_dir, filename)
            fig.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {filepath}")

        plt.show()
        plt.close(fig)

    def plot_histograms(self, df, bins=10, use_axial=False, save=False, save_dir='.',
                        overlay=False, opacity=0.6, color_palette=None):
        """
        Plot publication-ready histograms for each numeric column, grouped by
        either (charge, multiplicity) or (axial1, axial2), decoding categorical
        variables early so decoded (string) columns are excluded from numeric plots.
        Legend is drawn once outside subplots.
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
        COLOR_PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442']
        
        def _decode_df(df, vars_to_decode=None):
            """
            Decode specified variables in the DataFrame to their string labels.
            Returns a copy with decoded columns.
            """
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        # Decode axial and function variables upfront
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2', 'function'])

        # Determine grouping variables
        group_vars = ('axial1', 'axial2') if use_axial else ('charge', 'multiplicity')
        numeric_cols = df_dec.select_dtypes(include=[np.number]).columns
        plot_cols = [c for c in numeric_cols if c not in group_vars]
        combos = df_dec[list(group_vars)].drop_duplicates().to_records(index=False)
        palette = color_palette or COLOR_PALETTE

        if save and not os.path.exists(save_dir):
            os.makedirs(save_dir)

        for col in plot_cols:
            fig, ax = plt.subplots(figsize=(6, 4))
            patches = []
            labels = []
            for i, combo in enumerate(combos):
                val1, val2 = combo
                mask = ((df_dec[group_vars[0]] == val1) &
                        (df_dec[group_vars[1]] == val2) &
                        (df_dec[col] != -1))
                data = df_dec.loc[mask, col]
                label = f"{group_vars[0]}={val1}, {group_vars[1]}={val2}"
                if not data.empty:
                    counts, bins_edges, rects = ax.hist(
                        data, bins=bins, density=True, alpha=opacity,
                        color=palette[i % len(palette)], edgecolor='black'
                    )
                    patch = mpatches.Patch(color=palette[i % len(palette)], label=label)
                    patches.append(patch)
                    labels.append(label)
            ax.set_xlabel(col, fontsize=14)
            ax.set_ylabel('Density', fontsize=14)
            # Draw legend once outside plot
            fig.legend(handles=patches, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=12)
            plt.tight_layout()
            if save:
                fig.savefig(os.path.join(save_dir, f"{col}_overlay.png"), bbox_inches='tight')
            plt.close(fig)

    def plot_nested_histograms(self, df, bins=10, primary='axial', save=False, save_dir='.',
                               opacity=0.6, color_palette=None):
        """
        Plot nested histograms: subplots based on primary grouping and colored by secondary grouping,
        decoding categorical variables early so decoded columns are excluded from numeric plots.
        Legend is drawn once outside the subplots.
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
        COLOR_PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442']
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2', 'function'])

        if primary == 'axial':
            primary_vars = ('axial1', 'axial2')
            secondary_vars = ('charge', 'multiplicity')
        else:
            primary_vars = ('charge', 'multiplicity')
            secondary_vars = ('axial1', 'axial2')

        numeric_cols = df_dec.select_dtypes(include=[np.number]).columns
        plot_cols = [c for c in numeric_cols if c not in primary_vars + secondary_vars]
        combos_primary = df_dec[list(primary_vars)].drop_duplicates().to_records(index=False)
        combos_secondary = df_dec[list(secondary_vars)].drop_duplicates().to_records(index=False)
        palette = color_palette or COLOR_PALETTE

        if save and not os.path.exists(save_dir):
            os.makedirs(save_dir)

        for col in plot_cols:
            n = len(combos_primary)
            ncols = int(np.ceil(np.sqrt(n)))
            nrows = int(np.ceil(n / ncols))
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(4*ncols, 3*nrows),
                                     sharex=True, sharey=True)
            axes_flat = axes.flatten() if isinstance(axes, np.ndarray) else [axes]
            # Prepare legend patches for secondary groups
            legend_patches = []
            for i, secondary_combo in enumerate(combos_secondary):
                sval1, sval2 = secondary_combo
                label = f"{secondary_vars[0]}={sval1}, {secondary_vars[1]}={sval2}"
                patch = mpatches.Patch(color=palette[i % len(palette)], alpha=opacity, label=label)
                legend_patches.append(patch)

            for ax, primary_combo in zip(axes_flat, combos_primary):
                pval1, pval2 = primary_combo
                if col == 'homo_lumo_gap':
                    ax.set_xlabel('LUMO-HOMO [eV]')
                else:
                    ax.set_xlabel(col)
                ax.set_ylabel('Density')
                for i, secondary_combo in enumerate(combos_secondary):
                    sval1, sval2 = secondary_combo
                    mask = (
                        (df_dec[primary_vars[0]] == pval1) &
                        (df_dec[primary_vars[1]] == pval2) &
                        (df_dec[secondary_vars[0]] == sval1) &
                        (df_dec[secondary_vars[1]] == sval2) &
                        (df_dec[col] != -1)
                    )
                    data = df_dec.loc[mask, col]
                    if not data.empty:
                        ax.hist(
                            data, bins=bins, density=True, alpha=opacity,
                            color=palette[i % len(palette)], edgecolor='black'
                        )
            # Hide unused subplots
            for ax in axes_flat[n:]:
                ax.set_visible(False)
            # Draw single legend outside subplots
            fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=12)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            if save:
                fig.savefig(os.path.join(save_dir, f"{col}_nested_{primary}.png"), bbox_inches='tight')
            plt.close(fig)

    def plot_grouped_histogram(self, df, column, group_by='axial', bins=10, 
                               save=False, save_dir='.', opacity=0.6, 
                               color_palette=None, figsize=(12, 8)):
        """
        Plot histogram of a single column grouped by either axial ligand combinations 
        or charge-multiplicity combinations.
        
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
        COLOR_PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442']
        
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
            group_vars = ['axial1', 'axial2']
            df_plot['group_combo'] = df_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
            group_title = 'Axial Ligand Combinations'
        elif group_by == 'charge':
            group_vars = ['charge', 'multiplicity']
            df_plot['group_combo'] = df_plot.apply(lambda r: f"{r.charge}-{r.multiplicity}", axis=1)
            group_title = 'Charge-Multiplicity Combinations'
        else:
            raise ValueError("group_by must be either 'axial' or 'charge'")
        
        # Get unique combinations
        unique_groups = sorted(df_plot['group_combo'].unique())
        palette = color_palette or COLOR_PALETTE
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create histogram for each group
        for i, group in enumerate(unique_groups):
            group_data = df_plot[df_plot['group_combo'] == group][column]
            
            if not group_data.empty:
                ax.hist(group_data, bins=bins, alpha=opacity, 
                       color=palette[i % len(palette)], 
                       label=group, edgecolor='black', linewidth=0.5)
        
        # Set labels and title
        if column == 'homo_lumo_gap':
            ax.set_xlabel('LUMO-HOMO [eV]')
        elif column == 'lft_energy_delta':
            ax.set_xlabel('Δ [a.u.]')
        else:
            ax.set_xlabel(column)
        
        ax.set_ylabel('Frequency')
        
        # Add legend
        ax.legend(title=group_title, loc='upper right', frameon=False, fancybox=False, shadow=False)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save if requested
        if save:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            filename = f"{column}_grouped_by_{group_by}.png"
            fig.savefig(os.path.join(save_dir, filename), bbox_inches='tight')
            print(f"Plot saved as {filename}")
        
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
        COLOR_PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442']
        
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
        unique_groups = sorted(df_plot['group_combo'].unique())
        palette = color_palette or COLOR_PALETTE
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Store histogram data for contour calculation
        histogram_data = {}
        
        # Create histogram for each group
        for i, group in enumerate(unique_groups):
            group_data = df_plot[df_plot['group_combo'] == group][column]
            
            if not group_data.empty:
                color = palette[i % len(palette)]
                
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

    def plot_lft_scatter(self, df,
                          x_col='lft_energy_delta',
                          y_col='iron_natural_charge',
                          color_by='axial',  # 'axial', 'charge_multiplicity', or 'both'
                          save=False,
                          save_dir='.',
                          palette=None,
                          figsize=(10, 8),
                          marker_size=30,
                          draw_dividers=True,
                          light_lines=False,
                          number_quadrants=False):
        """
        Scatterplot of two numeric columns, colored or styled by specified grouping,
        with optional reference lines and quadrant numbering:
          - Vertical line at x=0
          - Horizontal lines at y=0.2, 0.6, 0.85

        Parameters
        ----------
        df : pandas.DataFrame
        x_col, y_col : str
        color_by : str
        save : bool
        save_dir : str
        palette : list
        figsize : tuple
            Figure size (width, height) in inches.
        marker_size : int
            Size of scatter points.
        draw_dividers : bool
            If True, draw vertical and horizontal reference lines.
        light_lines : bool
            If True (and draw_dividers is True), draw lines in a very light color.
        number_quadrants : bool
            If True (and draw_dividers is True), number the regions defined by the lines.
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
        COLOR_PALETTE = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00']
        MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X']
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2'])
        pal = palette or COLOR_PALETTE

        cols = [x_col, y_col, 'axial1', 'axial2', 'charge', 'multiplicity']
        df_plot = df_dec.loc[
            (df_dec[x_col] != -1) & (df_dec[y_col] != -1), cols
        ].dropna()

        df_plot['axial_combo'] = df_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
        df_plot['charge_combo'] = df_plot.apply(lambda r: f"{r.charge}-{r.multiplicity}", axis=1)

        if color_by in ('axial', 'both'):
            groups_axial = df_plot['axial_combo'].unique()
            color_map_axial = {g: pal[i % len(pal)] for i, g in enumerate(groups_axial)}
        if color_by in ('charge_multiplicity', 'both'):
            groups_charge = df_plot['charge_combo'].unique()
            color_map_charge = {g: pal[(i + (len(groups_axial) if 'groups_axial' in locals() else 0)) % len(pal)]
                                for i, g in enumerate(groups_charge)}
            marker_map = {g: MARKERS[i % len(MARKERS)] for i, g in enumerate(groups_charge)}

        fig, ax = plt.subplots(figsize=figsize)

        for _, row in df_plot.iterrows():
            if color_by == 'axial':
                c = color_map_axial[row['axial_combo']]
                m = 'o'
            elif color_by == 'charge_multiplicity':
                c = color_map_charge[row['charge_combo']]
                m = marker_map[row['charge_combo']]
            else:
                c = color_map_axial[row['axial_combo']]
                m = marker_map[row['charge_combo']]
            ax.scatter(row[x_col], row[y_col], color=c, marker=m,
                       alpha=0.8, edgecolor='black', s=marker_size)

        # Optional reference lines
        if draw_dividers:
            line_color = 'lightgrey' if light_lines else 'black'
            # vertical and horizontal lines
            ax.axvline(0, color=line_color, linestyle='--', linewidth=1)
            for y in (0.0, 0.3, 0.6, 0.9):
                ax.axhline(y, color=line_color, linestyle='--', linewidth=1)

            if number_quadrants:
                # get current plot limits
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                x_cuts = [xlim[0], 0, xlim[1]]
                y_cuts = [ylim[0], 0.2, 0.6, 0.85, ylim[1]]
                # annotate each region
                label = 1
                for j in range(len(y_cuts)-1):
                    for i in range(len(x_cuts)-1):
                        x_center = (x_cuts[i] + x_cuts[i+1]) / 2
                        y_center = (y_cuts[j] + y_cuts[j+1]) / 2
                        ax.text(x_center, y_center, str(label),
                                ha='center', va='center', fontsize=14, alpha=0.7)
                        label += 1

        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)

        from matplotlib.lines import Line2D
        legend_handles = []
        if color_by in ('axial', 'both'):
            legend_handles += [mpatches.Patch(color=color_map_axial[g], label=g)
                                for g in color_map_axial]
        if color_by in ('charge_multiplicity', 'both'):
            legend_handles += [Line2D([0], [0], marker=marker_map[g], color='w',
                                       markerfacecolor=color_map_charge[g], markersize=8,
                                       markeredgecolor='black', label=g)
                                for g in marker_map]
        ax.legend(handles=legend_handles, loc='upper left',
                  bbox_to_anchor=(1.02, 1), frameon=False)

        plt.tight_layout(rect=(0, 0, 0.85, 1))
        if save:
            os.makedirs(save_dir, exist_ok=True)
            fname = f"scatter_{x_col}_vs_{y_col}_{color_by}.png"
            fig.savefig(os.path.join(save_dir, fname), bbox_inches='tight')
        plt.close(fig)

    def plot_lft_scatter_with_violin(self, df,
                                   x_col='lft_energy_delta',
                                   y_col='iron_natural_charge',
                                   color_by='axial',  # 'axial', 'charge_multiplicity', or 'both'
                                   save=False,
                                   save_dir='.',
                                   palette=None,
                                   figsize=(14, 8),
                                   marker_size=30,
                                   draw_dividers=True,
                                   light_lines=False,
                                   number_quadrants=False,
                                   violin_alpha=0.7,
                                   violin_width=0.8,
                                   color_grading=True):
        """
        Scatterplot of two numeric columns with violin plots on the side showing
        y-axis density distributions for different groups, plus total density.
        
        Parameters
        ----------
        df : pandas.DataFrame
        x_col, y_col : str
        color_by : str
        save : bool
        save_dir : str
        palette : list
        figsize : tuple
            Figure size (width, height) in inches.
        marker_size : int
            Size of scatter points.
        draw_dividers : bool
            If True, draw vertical and horizontal reference lines.
        light_lines : bool
            If True (and draw_dividers is True), draw lines in a very light color.
        number_quadrants : bool
            If True (and draw_dividers is True), number the regions defined by the lines.
        violin_alpha : float
            Transparency of violin plots (0-1).
        violin_width : float
            Width scaling factor for violin plots.
        color_grading : bool
            If True, apply color intensity grading based on total_energies[0] values.
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
        COLOR_PALETTE = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00']
        MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X']
        
        def _decode_df(df, vars_to_decode=None):
            df_dec = df.copy()
            for var, mapping in REPLACEMENTS.items():
                if vars_to_decode is None or var in vars_to_decode:
                    if var in df_dec.columns:
                        df_dec[var] = df_dec[var].map(mapping)
            return df_dec
        
        df_dec = _decode_df(df, vars_to_decode=['axial1', 'axial2'])
        pal = palette or COLOR_PALETTE

        cols = [x_col, y_col, 'axial1', 'axial2', 'charge', 'multiplicity', 'total_energies[0]']
        df_plot = df_dec.loc[
            (df_dec[x_col] != -1) & (df_dec[y_col] != -1), cols
        ].dropna()

        df_plot['axial_combo'] = df_plot.apply(lambda r: f"{r.axial1}-{r.axial2}", axis=1)
        df_plot['charge_combo'] = df_plot.apply(lambda r: f"{r.charge}-{r.multiplicity}", axis=1)

        if color_by in ('axial', 'both'):
            groups_axial = df_plot['axial_combo'].unique()
            color_map_axial = {g: pal[i % len(pal)] for i, g in enumerate(groups_axial)}
        if color_by in ('charge_multiplicity', 'both'):
            groups_charge = df_plot['charge_combo'].unique()
            color_map_charge = {g: pal[(i + (len(groups_axial) if 'groups_axial' in locals() else 0)) % len(pal)]
                                for i, g in enumerate(groups_charge)}
            marker_map = {g: MARKERS[i % len(MARKERS)] for i, g in enumerate(groups_charge)}

        # Create figure with subplots: scatter plot (left) and violin plot (right)
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.15)
        
        # Main scatter plot
        ax_scatter = fig.add_subplot(gs[0, 0])
        
        # Violin plot
        ax_violin = fig.add_subplot(gs[0, 1], sharey=ax_scatter)

        # Plot scatter points with intensity based on total_energy
        import matplotlib.colors as mcolors
        
        # Calculate energy min/max/range for each axial combination separately
        axial_combinations = ['CYS-HOH', 'HIS-HIS', 'HIS-HOH', 'HIS-MET', 'HIS-OXY']
        energy_stats = {}
        
        # Create marker mapping for axial combinations to distinguish groups
        axial_marker_map = {
            'CYS-HOH': 'o',    # circle
            'HIS-HIS': 's',    # square
            'HIS-HOH': '^',    # triangle up
            'HIS-MET': 'D',    # diamond
            'HIS-OXY': 'v'     # triangle down
        }
        
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
        
        for _, row in df_plot.iterrows():
            if color_by == 'axial':
                base_color = color_map_axial[row['axial_combo']]
                m = axial_marker_map.get(row['axial_combo'], 'o')  # Use axial-specific markers
            elif color_by == 'charge_multiplicity':
                base_color = color_map_charge[row['charge_combo']]
                m = marker_map[row['charge_combo']]
            else:
                base_color = color_map_axial[row['axial_combo']]
                # For 'both' mode, use axial markers to distinguish axial groups
                m = axial_marker_map.get(row['axial_combo'], 'o')
            
            # Calculate intensity based on total_energy using axial-specific ranges (if enabled)
            if color_grading:
                axial_combo = row['axial_combo']
                if axial_combo in energy_stats and energy_stats[axial_combo]['range'] > 0:
                    stats = energy_stats[axial_combo]
                    # Amplify small energy differences (< 1 Hartree) by factor of 10
                    amplified_diff = 10 * (row['total_energies[0]'] - stats['min'])
                    amplified_range = 10 * stats['range']
                    normalized_energy = amplified_diff / amplified_range
                    # Clamp to [0, 1] range in case amplification causes overflow
                    normalized_energy = min(1.0, max(0.0, normalized_energy))
                    # Apply exponential scaling for more dramatic effect (higher energies more pronounced)
                    intensity = 0.1 + 0.9 * (normalized_energy ** 2)
                else:
                    intensity = 1.0
                
                # Convert base color to RGB and apply intensity
                base_rgb = mcolors.to_rgb(base_color)
                # Scale the color intensity (higher energy = more saturated)
                final_color = tuple(intensity * c + (1 - intensity) * 1.0 for c in base_rgb)
            else:
                # Use base color without intensity grading
                final_color = base_color
            
            ax_scatter.scatter(row[x_col], row[y_col], color=final_color, marker=m,
                              alpha=0.8, edgecolor='black', s=marker_size)

        # Optional reference lines on scatter plot
        if draw_dividers:
            line_color = 'lightgrey' if light_lines else 'black'
            # vertical and horizontal lines
            ax_scatter.axvline(0, color=line_color, linestyle='--', linewidth=1)
            for y in (0.0, 0.3, 0.6, 0.9):
                ax_scatter.axhline(y, color=line_color, linestyle='--', linewidth=1)
                ax_violin.axhline(y, color=line_color, linestyle='--', linewidth=1, alpha=0.5)

            if number_quadrants:
                # get current plot limits
                xlim = ax_scatter.get_xlim()
                ylim = ax_scatter.get_ylim()
                x_cuts = [xlim[0], 0, xlim[1]]
                y_cuts = [ylim[0], 0.2, 0.6, 0.85, ylim[1]]
                # annotate each region
                label = 1
                for j in range(len(y_cuts)-1):
                    for i in range(len(x_cuts)-1):
                        x_center = (x_cuts[i] + x_cuts[i+1]) / 2
                        y_center = (y_cuts[j] + y_cuts[j+1]) / 2
                        ax_scatter.text(x_center, y_center, str(label),
                                       ha='center', va='center', fontsize=14, alpha=0.7)
                        label += 1

        # Create violin plots
        # First, create data for violin plots based on grouping
        if color_by == 'axial':
            grouping_col = 'axial_combo'
            groups = groups_axial
            color_map = color_map_axial
        elif color_by == 'charge_multiplicity':
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
            group_data = df_plot[df_plot[grouping_col] == group][y_col].dropna()
            # Need at least 2 data points for violin plot
            if len(group_data) >= 2:
                violin_data.append(group_data)
                violin_labels.append(group)
                violin_colors.append(color_map[group])
        
        # Add total density violin plot (only if we have enough data)
        total_data = df_plot[y_col].dropna()
        if len(total_data) >= 2:
            violin_data.append(total_data)
            violin_labels.append('Total')
            violin_colors.append('gray')

        # Create violin plots
        if len(violin_data) > 0:
            positions = np.arange(len(violin_data))
            
            # Create violin plot parts
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
            
            # Decode violin plot labels using REPLACEMENTS dictionary
            decoded_labels = []
            for label in violin_labels:
                if label == 'Total':
                    decoded_labels.append('Total')
                elif '-' in label:
                    # Handle axial combinations like 'CYS-HOH'
                    parts = label.split('-')
                    if len(parts) == 2:
                        # Already decoded by _decode_df, so keep as is
                        decoded_labels.append(label)
                    else:
                        decoded_labels.append(label)
                else:
                    # Handle charge-multiplicity combinations or other labels
                    decoded_labels.append(label)
            
            # Set violin plot labels
            ax_violin.set_xticks(positions)
            ax_violin.set_xticklabels(decoded_labels, rotation=45, ha='right')
            if color_by == 'axial':
                ax_violin.set_xlabel('Axial Ligands')
            elif color_by == 'charge_multiplicity':
                ax_violin.set_xlabel('Electronic State')
            else:
                ax_violin.set_xlabel('Groups')

        # Set labels and title
        labels=True
        if not labels:
            ax_scatter.set_xlabel(x_col)
            ax_scatter.set_ylabel(y_col)
        else:
            ax_scatter.set_xlabel("Δ [eV]")
            ax_scatter.set_ylabel("Iron Charge (NBO) [a.u.]")
        
        # Ensure y-axis tick labels are visible on scatter plot
        ax_scatter.tick_params(axis='y', which='both', labelleft=True, labelright=False)
        ax_scatter.yaxis.set_tick_params(which='both', labelleft=True)
        # Force matplotlib to show y-axis labels
        plt.setp(ax_scatter.get_yticklabels(), visible=True)
        
        # Remove y-axis labels from violin plot since it shares with scatter plot
        ax_violin.tick_params(axis='y', which='both', labelleft=False, labelright=False)

        # Create legend for scatter plot
        from matplotlib.lines import Line2D
        legend_handles = []
        if color_by in ('axial', 'both'):
            legend_handles += [mpatches.Patch(color=color_map_axial[g], label=g)
                                for g in color_map_axial]
        if color_by in ('charge_multiplicity', 'both'):
            legend_handles += [Line2D([0], [0], marker=marker_map[g], color='w',
                                       markerfacecolor=color_map_charge[g], markersize=8,
                                       markeredgecolor='black', label=g)
                                for g in marker_map]
        # Add legend with appropriate title
        legend_title = None
        if color_by == 'charge_multiplicity':
            legend_title = 'Charge-Multiplicity'
        elif color_by == 'axial':
            legend_title = 'Axial Ligands'
        elif color_by == 'both':
            legend_title = 'Groups'
        
        ax_scatter.legend(handles=legend_handles, loc='upper right',
                         bbox_to_anchor=(0.98, 0.98), frameon=False, fancybox=False, shadow=False,
                         title=legend_title)

        plt.tight_layout()
        if save:
            os.makedirs(save_dir, exist_ok=True)
            fname = f"scatter_violin_{x_col}_vs_{y_col}_{color_by}.png"
            fig.savefig(os.path.join(save_dir, fname), bbox_inches='tight', dpi=300)
        plt.show()
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
                                        titles=None):
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
        COLOR_PALETTE = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#E69F00']
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

        # Combine datasets to get consistent color mapping
        combined_df = pd.concat([df1_plot, df2_plot], ignore_index=True)
        
        # Create color mappings for all schemes
        groups_axial = combined_df['axial_combo'].unique()
        color_map_axial = {g: pal[i % len(pal)] for i, g in enumerate(groups_axial)}
        
        groups_charge = combined_df['charge_combo'].unique()
        color_map_charge = {g: pal[(i + len(groups_axial)) % len(pal)] for i, g in enumerate(groups_charge)}
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
            
            ax.set_xlabel("Δ [eV]", fontsize=16)
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
            fname = f"plots/dual_scatter_violin_{x_col}_vs_{y_col}_{color_by_left}_{color_by_right}.png"
            fig.savefig(os.path.join(save_dir, fname), bbox_inches='tight', dpi=300)
        plt.show()
        plt.close(fig)

    def plot_iron_d_orbital_energy_differences(self, df, group_by='charge_multiplicity', 
                                               save_dir='plots', decode_dicts=None,
                                               figsize=(16, 12), show_groupings='separate', verbose=False):
        """
        Plots iron d-orbital energy differences with comprehensive comparisons.
        
        Creates multiple subplot comparisons showing:
        1. Individual orbital energy differences (dx2y2_dz2_diff, dz2_dxy_diff, etc.)
        2. Orbital distortion type distribution with grouping breakdowns
        3. Grouped comparisons by charge-multiplicity, axial ligands, and PDB-ID
        
        Parameters:
        -----------
        df : pd.DataFrame
            DataFrame containing orbital energy difference columns
        group_by : str
            Primary grouping variable ('charge_multiplicity', 'axials', 'PDB-ID')
        save_dir : str
            Directory to save plots
        decode_dicts : dict
            Dictionary for decoding categorical variables
        figsize : tuple
            Figure size (width, height)
        show_groupings : str
            How to show groupings: 'separate' (two separate barplots), 'combined' (nested stacked bars), 
            'grouped' (side-by-side grouped bars), 'both' (all separate plots)
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
        
        # Create grid layout based on show_groupings parameter
        if show_groupings == 'separate':
            fig = plt.figure(figsize=(16, 10))
            gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.4)
        elif show_groupings == 'combined':
            fig = plt.figure(figsize=(15, 8))
            gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.4)
        elif show_groupings == 'grouped':
            fig = plt.figure(figsize=(15, 8))
            gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.4)
        else:  # 'both'
            fig = plt.figure(figsize=(18, 12))
            gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.4)
        
        colors_distortion = plt.cm.viridis(np.linspace(0, 1, 4))  # 4 distortion types
        
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
                if isinstance(sample_value, (int, float)):
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
        
        # Helper function to create stacked bar plot for groupings
        def create_stacked_barplot(ax, distortion_counts, group_col, colormap):
            """Create stacked barplot showing breakdown by group_col for each distortion type"""
            unique_groups = df_valid[group_col].unique()
            
            print(f"DEBUG: Creating stacked barplot for {group_col}")
            print(f"DEBUG: Unique groups found: {unique_groups}")
            print(f"DEBUG: Distortion counts: {distortion_counts}")
            print(f"DEBUG: Distortion labels: {list(distortion_counts.index)}")
            
            if len(unique_groups) == 0:
                print(f"DEBUG: No unique groups found for {group_col}")
                return
            
            group_colors = colormap(np.linspace(0, 1, len(unique_groups)))
            
            # Create data structure for stacked bars
            stacked_data = {}
            for distortion in distortion_counts.index:
                stacked_data[distortion] = {}
                subset = df_valid[df_valid['distortion_decoded'] == distortion]
                group_counts = subset[group_col].value_counts()
                print(f"DEBUG: For {distortion}, subset size: {len(subset)}")
                print(f"DEBUG: Group counts: {group_counts.to_dict()}")
                for group in unique_groups:
                    stacked_data[distortion][group] = group_counts.get(group, 0)
            
            # FIXED: Use numeric positions for bars
            x_pos = np.arange(len(distortion_counts))
            distortion_labels = list(distortion_counts.index)
            
            print(f"DEBUG: x_pos: {x_pos}")
            print(f"DEBUG: distortion_labels: {distortion_labels}")
            
            # Create stacked bars with proper positioning
            bottom = np.zeros(len(distortion_counts))
            for i, group in enumerate(unique_groups):
                values = [stacked_data[dist][group] for dist in distortion_counts.index]
                print(f"DEBUG: Group {group}, values: {values}, bottom: {bottom}")
                if sum(values) > 0:  # Only plot if there's data
                    bars = ax.bar(x_pos, values, bottom=bottom, width=0.6,
                                 color=group_colors[i], label=str(group), alpha=0.8)
                    print(f"DEBUG: Created bars for {group}: {[bar.get_height() for bar in bars]}")
                    bottom += np.array(values)
            
            # Set proper x-axis labels
            ax.set_xticks(x_pos)
            ax.set_xticklabels(distortion_labels, rotation=45, fontsize=12)
            ax.set_ylabel('Count', fontsize=14)
            ax.tick_params(axis='y', labelsize=12)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, title_fontsize=14)
            
            # Add total count labels on bars
            for i, (distortion, total) in enumerate(distortion_counts.items()):
                ax.annotate(f'{int(total)}',
                           xy=(x_pos[i], total),
                           xytext=(0, 3), textcoords="offset points",
                           ha='center', va='bottom', fontsize=14)
                           
            print(f"DEBUG: Final bottom values: {bottom}")
        
        # Helper function for side-by-side grouped bars (improved combined visualization)
        def create_grouped_barplot(ax, distortion_counts):
            """Create side-by-side grouped bars showing both charge-multiplicity and axial ligands"""
            distortion_names = distortion_counts.index
            x = np.arange(len(distortion_names))
            
            # Get unique groups
            charge_mult_groups = sorted(df_valid['charge_multiplicity'].unique())
            axial_groups = sorted(df_valid['axials'].value_counts().head(5).index.tolist())  # Top 5 axials
            
            print(f"DEBUG: Creating grouped barplot")
            print(f"DEBUG: Distortion names: {list(distortion_names)}")
            print(f"DEBUG: Charge-mult groups: {charge_mult_groups}")
            print(f"DEBUG: Top axial groups: {axial_groups}")
            print(f"DEBUG: x positions: {x}")
            
            # Calculate bar width and positions
            total_bars = len(charge_mult_groups) + len(axial_groups)
            bar_width = 0.8 / total_bars
            
            # Plot charge-multiplicity bars
            charge_colors = plt.cm.Set1(np.linspace(0, 1, len(charge_mult_groups)))
            for i, charge_mult in enumerate(charge_mult_groups):
                values = []
                for distortion in distortion_names:
                    subset = df_valid[
                        (df_valid['distortion_decoded'] == distortion) & 
                        (df_valid['charge_multiplicity'] == charge_mult)
                    ]
                    values.append(len(subset))
                
                print(f"DEBUG: Charge {charge_mult}, values: {values}")
                if sum(values) > 0:
                    offset = (i - len(charge_mult_groups)/2 + 0.5) * bar_width
                    print(f"DEBUG: Plotting charge {charge_mult} with offset {offset}")
                    bars = ax.bar(x + offset, values, bar_width, label=f'Charge: {charge_mult}', 
                                 color=charge_colors[i], alpha=0.8)
                    print(f"DEBUG: Created charge bars: {[bar.get_height() for bar in bars]}")
            
            # Plot axial ligand bars
            axial_colors = plt.cm.Set2(np.linspace(0, 1, len(axial_groups)))
            for i, axial in enumerate(axial_groups):
                values = []
                for distortion in distortion_names:
                    subset = df_valid[
                        (df_valid['distortion_decoded'] == distortion) & 
                        (df_valid['axials'] == axial)
                    ]
                    values.append(len(subset))
                
                print(f"DEBUG: Axial {axial}, values: {values}")
                if sum(values) > 0:
                    offset = (i + len(charge_mult_groups) - total_bars/2 + 0.5) * bar_width
                    print(f"DEBUG: Plotting axial {axial} with offset {offset}")
                    bars = ax.bar(x + offset, values, bar_width, label=f'Axial: {axial}', 
                                 color=axial_colors[i], alpha=0.6, hatch='//')
                    print(f"DEBUG: Created axial bars: {[bar.get_height() for bar in bars]}")
            
            ax.set_xlabel('Distortion Type', fontsize=14)
            ax.set_ylabel('Count', fontsize=14)
            ax.set_xticks(x)
            ax.set_xticklabels(distortion_names, rotation=45, fontsize=12)
            ax.tick_params(axis='y', labelsize=12)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, title_fontsize=14)
        
        # Helper function for combined visualization using patterns (original approach)
        def create_combined_barplot(ax, distortion_counts):
            """Create combined visualization with nested stacking and patterns"""
            charge_mult_groups = df_valid['charge_multiplicity'].unique()
            charge_colors = plt.cm.Set1(np.linspace(0, 1, len(charge_mult_groups)))
            
            print(f"DEBUG: Creating combined barplot")
            print(f"DEBUG: Charge groups: {list(charge_mult_groups)}")
            
            # Get top axial combinations to avoid too many patterns
            axial_counts = df_valid['axials'].value_counts()
            top_axials = axial_counts.head(6).index.tolist()  # Top 6 most common
            patterns = ['', '////', '....', '+++', 'xxx', '|||']
            
            print(f"DEBUG: Top axials: {top_axials}")
            
            x_pos = np.arange(len(distortion_counts))
            bar_width = 0.8
            
            print(f"DEBUG: x_pos: {x_pos}, distortions: {list(distortion_counts.index)}")
            
            # Create nested stacked bars
            for i, distortion in enumerate(distortion_counts.index):
                subset = df_valid[df_valid['distortion_decoded'] == distortion]
                bottom = 0
                print(f"DEBUG: Processing distortion {distortion}, subset size: {len(subset)}")
                
                for j, charge_mult in enumerate(charge_mult_groups):
                    charge_subset = subset[subset['charge_multiplicity'] == charge_mult]
                    if len(charge_subset) == 0:
                        continue
                    
                    print(f"DEBUG: Charge {charge_mult} subset size: {len(charge_subset)}")
                    
                    # Within each charge group, stack by axial ligands
                    axial_counts_in_charge = charge_subset['axials'].value_counts()
                    for k, axial in enumerate(top_axials):
                        if axial in axial_counts_in_charge:
                            count = axial_counts_in_charge[axial]
                            pattern = patterns[k] if k < len(patterns) else ''
                            print(f"DEBUG: Adding bar at x={x_pos[i]}, height={count}, bottom={bottom}")
                            bar = ax.bar(x_pos[i], count, bottom=bottom, 
                                        color=charge_colors[j], alpha=0.8,
                                        hatch=pattern, width=bar_width,
                                        label=f'{charge_mult}_{axial}' if i == 0 else "")
                            print(f"DEBUG: Created bar with height: {bar[0].get_height()}")
                            bottom += count
            
            ax.set_xticks(x_pos)
            ax.set_xticklabels(distortion_counts.index, rotation=45)
            ax.set_ylabel('Count')
            
            # Create custom legend
            from matplotlib.patches import Rectangle
            legend_elements = []
            for j, charge_mult in enumerate(charge_mult_groups):
                legend_elements.append(Rectangle((0,0),1,1, facecolor=charge_colors[j], 
                                               alpha=0.8, label=f'Charge: {charge_mult}'))
            for k, axial in enumerate(top_axials):
                if k < len(patterns):
                    legend_elements.append(Rectangle((0,0),1,1, facecolor='gray', 
                                                   hatch=patterns[k], alpha=0.5, 
                                                   label=f'Axial: {axial}'))
            ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), 
                     loc='upper left', fontsize=12, title_fontsize=14)
        
        # Validate required columns for grouping plots
        required_grouping_cols = ['charge_multiplicity', 'axials']
        missing_grouping_cols = [col for col in required_grouping_cols if col not in df_valid.columns]
        
        if missing_grouping_cols:
            if verbose:
                print(f"DEBUG: Missing required grouping columns: {missing_grouping_cols}")
                print("DEBUG: Cannot create grouping plots, will only create scatter plot")
            # Only create scatter plot
            ax_scatter = fig.add_subplot(gs[0, 0])
        else:
            # We have the required grouping columns, proceed with barplots
            distortion_counts = df_valid['distortion_decoded'].value_counts()
            
            if verbose:
                print(f"DEBUG: Distortion counts for plotting: {distortion_counts.to_dict()}")
                print(f"DEBUG: Charge-multiplicity counts: {df_valid['charge_multiplicity'].value_counts().to_dict()}")
                print(f"DEBUG: Axials counts: {df_valid['axials'].value_counts().to_dict()}")
            
            if show_groupings == 'separate':
                # 1. Charge-multiplicity breakdown
                ax_charge = fig.add_subplot(gs[0, 0])
                create_stacked_barplot(ax_charge, distortion_counts, 'charge_multiplicity', plt.cm.Set1)
                
                # 2. Axial ligand breakdown
                ax_axial = fig.add_subplot(gs[0, 1])
                create_stacked_barplot(ax_axial, distortion_counts, 'axials', plt.cm.Set2)
                
                # 3. Scatter plot
                ax_scatter = fig.add_subplot(gs[1, :])
                
            elif show_groupings == 'combined':
                # 1. Combined visualization (stacked with patterns)
                ax_combined = fig.add_subplot(gs[0, 0])
                create_combined_barplot(ax_combined, distortion_counts)
                
                # 2. Scatter plot
                ax_scatter = fig.add_subplot(gs[0, 1])
                
            elif show_groupings == 'grouped':
                # 1. Grouped side-by-side visualization
                ax_grouped = fig.add_subplot(gs[0, 0])
                create_grouped_barplot(ax_grouped, distortion_counts)
                
                # 2. Scatter plot
                ax_scatter = fig.add_subplot(gs[0, 1])
                
            else:  # 'both'
                # 1. Charge-multiplicity breakdown
                ax_charge = fig.add_subplot(gs[0, 0])
                create_stacked_barplot(ax_charge, distortion_counts, 'charge_multiplicity', plt.cm.Set1)
                
                # 2. Axial ligand breakdown
                ax_axial = fig.add_subplot(gs[0, 1])
                create_stacked_barplot(ax_axial, distortion_counts, 'axials', plt.cm.Set2)
                
                # 3. Combined visualization
                ax_combined = fig.add_subplot(gs[0, 2])
                create_combined_barplot(ax_combined, distortion_counts)
                
                # 4. Scatter plot
                ax_scatter = fig.add_subplot(gs[1, :])
        
        # Scatter plot showing orbital energy relationships
        scatter_cols = ['dx2y2_dz2_diff', 'eg_t2g_diff']
        missing_scatter_cols = [col for col in scatter_cols if col not in df_valid.columns]
        
        if missing_scatter_cols:
            if verbose:
                print(f"DEBUG: Missing scatter plot columns: {missing_scatter_cols}")
                print("DEBUG: Cannot create scatter plot")
            ax_scatter.text(0.5, 0.5, 'Scatter plot data not available', 
                           ha='center', va='center', transform=ax_scatter.transAxes)
        else:
            if verbose:
                print(f"DEBUG: Scatter plot data available, {len(df_valid)} points")
            
            if distortion_col and has_distortion_data:
                # Plot by distortion type
                unique_distortions = df_valid[distortion_col].dropna().unique()
                colors = plt.cm.viridis(np.linspace(0, 1, len(unique_distortions)))
                
                if verbose:
                    print(f"DEBUG: Creating scatter plot with {len(unique_distortions)} distortion types")
                
                for i, distortion_type in enumerate(sorted(unique_distortions)):
                    mask = df_valid[distortion_col] == distortion_type
                    subset = df_valid[mask]
                    if len(subset) > 0:
                        # Get the display label using the appropriate mapping
                        if isinstance(distortion_type, (int, float)):
                            display_label = numeric_distortion_map.get(distortion_type, f'Type {distortion_type}')
                        else:
                            display_label = string_distortion_map.get(str(distortion_type).lower(), f'Type {distortion_type}')
                        
                        ax_scatter.scatter(subset['dx2y2_dz2_diff'], subset['eg_t2g_diff'], 
                                          c=[colors[i]], alpha=0.7, s=50,
                                          label=display_label)
                
                ax_scatter.legend(fontsize=14, loc='upper right', title_fontsize=16)
            else:
                # Plot all points without distortion type coloring
                if verbose:
                    print("DEBUG: Creating scatter plot without distortion type coloring")
                ax_scatter.scatter(df_valid['dx2y2_dz2_diff'], df_valid['eg_t2g_diff'], 
                                  alpha=0.7, s=50, label='All Data')
                ax_scatter.legend(fontsize=14, loc='upper right', title_fontsize=16)
            
            ax_scatter.set_xlabel(r'$d_{x^2-y^2} - d_{z^2}$ Difference', fontsize=16)
            ax_scatter.set_ylabel(r'$e_{g} - t_{2g}$ Difference', fontsize=16)
        
        
        # Save plot
        filename = f'{save_dir}/iron_d_orbital_energy_analysis_{group_by}_{show_groupings}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"Iron d-orbital energy analysis plot saved to: {filename}")
        
        if verbose:
            # Print summary statistics
            print("\nSummary Statistics:")
            print(f"Total structures analyzed: {len(df_valid)}")
            if 'orbital_distortion_type' in df_valid.columns:
                print("\nDistortion type distribution:")
                for distortion, count in df_valid['distortion_decoded'].value_counts().items():
                    print(f"  {distortion}: {count}")
                
                if show_groupings in ['separate', 'both']:
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
            colors_cm = plt.cm.Set1(np.linspace(0, 1, len(merged_df['charge_multiplicity'].unique())))
            colors_axial = plt.cm.Set2(np.linspace(0, 1, len(merged_df['axials'].unique()) if 'axials' in merged_df.columns else 4))
            colors_function = plt.cm.viridis(np.linspace(0, 1, len(merged_df['function'].unique()) if 'function' in merged_df.columns else 6))
            
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
            
            # Scatter plot with error bars
            for i, group in enumerate(merged_df[group_var].unique()):
                group_data = merged_df[merged_df[group_var] == group]
                ax1.errorbar(group_data.index, group_data['E°approx (V)'], 
                            yerr=group_data['std_dev (V)'], fmt='o',
                            color=scatter_colors[i % len(scatter_colors)], 
                            label=group, alpha=0.7, capsize=3, markersize=6)
            
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
                    bars = ax8.bar(x_pos, std_devs, alpha=0.7, color=plt.cm.Reds(np.linspace(0.3, 1, len(pdb_ids))))
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
            filename = f'{save_dir}/reduction_potentials_comprehensive_{primary_grouping}.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            plt.show()

            print(f"Plot saved to: {filename}")
            return merged_df, fig

        except Exception as e:
            print(f"Error in plotting reduction potentials: {e}")
            return None, None


    def plot_high_std_dev_pdbs_detailed(self, reduction_potentials_file='tables/reduction_potentials.csv',
                                      processed_data_file='tables/processed_output.csv',
                                      save_dir='plots', decode_dicts=None,
                                      figsize=(20, 15), top_n=6):
        """
        Creates detailed plots for PDB structures with the highest standard deviation in reduction potentials.
        Shows individual state transition calculations for each high-variance PDB.
        
        Parameters:
        -----------
        reduction_potentials_file : str
            Path to the reduction potentials CSV file
        processed_data_file : str
            Path to processed data CSV file
        save_dir : str
            Directory to save plots
        decode_dicts : dict
            Dictionary for decoding categorical variables
        figsize : tuple
            Figure size (width, height)
        top_n : int
            Number of top high-std-dev PDBs to analyze
        """
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        import ast
        
        os.makedirs(save_dir, exist_ok=True)
        
        try:
            # Load data
            redox_df = pd.read_csv(reduction_potentials_file)
            processed_df = pd.read_csv(processed_data_file)
            
            # Extract PDB-ID from filename in processed data for merging
            processed_df['structure_id'] = processed_df['file_name'].str[:4]
            
            # Merge datasets to get additional metadata
            merged_df = redox_df.merge(
                processed_df[['structure_id', 'charge', 'multiplicity', 'axial1', 'axial2', 'function', 'PDB-ID']].drop_duplicates('structure_id'),
                on='structure_id', how='left'
            )
            
            # Create axials column
            if 'axial1' in merged_df.columns and 'axial2' in merged_df.columns:
                merged_df['axials'] = merged_df['axial1'].astype(str) + '-' + merged_df['axial2'].astype(str)
            
            # Select PDBs with highest standard deviation (must have multiple measurements)
            high_std_pdbs = merged_df[merged_df['std_dev (V)'] > 0].nlargest(top_n, 'std_dev (V)')
            
            if len(high_std_pdbs) == 0:
                print("No PDBs with multiple measurements found")
                return None, None
            
            # Create comprehensive plot
            fig = plt.figure(figsize=figsize)
            
            # Calculate grid layout based on number of PDBs
            n_cols = min(3, top_n)
            n_rows = (top_n + n_cols - 1) // n_cols + 1  # +1 for summary row
            
            gs = fig.add_gridspec(n_rows, n_cols, hspace=0.4, wspace=0.3)
            
            # Color schemes for different transition types
            all_transitions = set()
            for _, row in high_std_pdbs.iterrows():
                try:
                    pair_details = ast.literal_eval(row['pair_details'])
                    for pair in pair_details:
                        transition = f"({pair['red_charge']},{pair['red_mult']})→({pair['ox_charge']},{pair['ox_mult']})"
                        all_transitions.add(transition)
                except:
                    continue
            
            transition_colors = {transition: plt.cm.Set3(i/len(all_transitions)) 
                               for i, transition in enumerate(sorted(all_transitions))}
            
            # Plot individual PDB analyses
            for idx, (_, pdb_row) in enumerate(high_std_pdbs.iterrows()):
                row_idx = idx // n_cols
                col_idx = idx % n_cols
                ax = fig.add_subplot(gs[row_idx, col_idx])
                
                try:
                    pair_details = ast.literal_eval(pdb_row['pair_details'])
                    
                    # Extract individual potentials and transition types
                    potentials = []
                    transitions = []
                    colors = []
                    
                    for pair in pair_details:
                        potential = pair['potential_V']
                        transition = f"({pair['red_charge']},{pair['red_mult']})→({pair['ox_charge']},{pair['ox_mult']})"
                        
                        potentials.append(potential)
                        transitions.append(transition)
                        colors.append(transition_colors[transition])
                    
                    # Create bar plot of individual calculations
                    x_positions = range(len(potentials))
                    bars = ax.bar(x_positions, potentials, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
                    
                    # Add value labels on bars
                    for i, (bar, pot) in enumerate(zip(bars, potentials)):
                        height = bar.get_height()
                        ax.text(bar.get_x() + bar.get_width()/2., height + (0.1 if height > 0 else -0.1),
                               f'{pot:.2f}', ha='center', va='bottom' if height > 0 else 'top', 
                               fontsize=12, fontweight='bold')
                    
                    # Customize plot
                    ax.set_xlabel('Individual Calculations')
                    ax.set_ylabel('E° (V)')
                    ax.set_xticks(x_positions)
                    ax.set_xticklabels([f'#{i+1}' for i in x_positions], fontsize=12)
                    ax.grid(True, alpha=0.3)
                    
                    # Add horizontal line for mean
                    ax.axhline(y=pdb_row['E°approx (V)'], color='red', linestyle='--', 
                             alpha=0.8, linewidth=2, label=f'Mean: {pdb_row["E°approx (V)"]:.2f}V')
                    
                    # Add shaded region for std dev
                    mean_val = pdb_row['E°approx (V)']
                    std_val = pdb_row['std_dev (V)']
                    ax.axhspan(mean_val - std_val, mean_val + std_val, alpha=0.2, color='red')
                    
                    ax.legend(fontsize=12)
                    
                except Exception as e:
                    ax.text(0.5, 0.5, f'Error parsing data\nfor {pdb_row["structure_id"]}', 
                           ha='center', va='center', transform=ax.transAxes, fontsize=14)
            
            # Summary plot showing transition type distributions
            summary_ax = fig.add_subplot(gs[-1, :])
            
            # Collect all transition data
            transition_data = {trans: [] for trans in all_transitions}
            
            for _, row in high_std_pdbs.iterrows():
                try:
                    pair_details = ast.literal_eval(row['pair_details'])
                    for pair in pair_details:
                        transition = f"({pair['red_charge']},{pair['red_mult']})→({pair['ox_charge']},{pair['ox_mult']})"
                        transition_data[transition].append(pair['potential_V'])
                except:
                    continue
            
            # Create box plot of transitions
            box_data = [data for data in transition_data.values() if len(data) > 0]
            box_labels = [trans for trans, data in transition_data.items() if len(data) > 0]
            
            if box_data:
                bp = summary_ax.boxplot(box_data, labels=box_labels, patch_artist=True)
                
                for patch, label in zip(bp['boxes'], box_labels):
                    patch.set_facecolor(transition_colors[label])
                    patch.set_alpha(0.7)
                
                summary_ax.set_ylabel('E° (V)')
                summary_ax.set_xlabel('Transition Type')
                plt.setp(summary_ax.get_xticklabels(), rotation=45, fontsize=12)
                summary_ax.grid(True, alpha=0.3)
            
            
            # Save plot
            filename = f'{save_dir}/high_std_dev_pdbs_detailed_analysis.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            plt.show()
            
            # Print summary
            print(f"High standard deviation PDB analysis saved to: {filename}")
            print(f"\nAnalyzed PDBs:")
            for _, row in high_std_pdbs.iterrows():
                print(f"  {row['structure_id']}: Mean = {row['E°approx (V)']:.3f}V, "
                      f"σ = {row['std_dev (V)']:.3f}V, {row['num_red_ox_pairs']} measurements")
            
            print(f"\nTransition types found: {len(all_transitions)}")
            for transition in sorted(all_transitions):
                count = sum(1 for _, row in high_std_pdbs.iterrows() 
                          if transition in str(row['pair_details']))
                print(f"  {transition}: {count} occurrences")
            
            return high_std_pdbs, fig
            
        except Exception as e:
            print(f"Error in plotting high std dev PDBs: {e}")
            return None, None



#HEATMAP LATER AXIAL OVERVIEW
p = plots_class()
#  (1) Option A: pass decode_dicts directly to get_plot_func
fig, ax = plt.subplots(figsize=(8, 6))
heatmap_func = p.get_plot_func(
    plot_type="auto",
    df=df,
    col_x="axial1",
    col_y="axial2",
    decode_dicts=replacements
)
heatmap_func(ax)
ax.set_title("Decoded heatmap of axial1 vs axial2")
plt.tight_layout()
plt.show()

# PLOT ALL COMBINATIONS OF DIFFERENT VARIABLES AS HISTOGRAMS
do_all = False
if do_all:
    # 1) Instantiate & decode (dropping file_name early)
    p = plots_class()
    p.split_legend_flag = True

    df_decoded = (
        df_plots
        .drop(columns=["file_name"])
        .pipe(lambda df: p.preprocess_dataframe(df, decode_dicts=replacements))
        .dropna()
    )

    # 2) Build a list of columns you actually want to plot
    exclude = {"file_name", "# PDB"}
    selected = [c for c in df_decoded.columns if c not in exclude]

    # 3) Call with selected_columns + histogram mode
    p.plot_all_scatter_combinations(
        df=df_decoded,
        annotate_counts=False,
        color_intensity=False,
        group_counts=True,
        numeric_only=False,
        selected_columns=selected,
        cc_plot="hist"        # ← switch to histograms (instead of KDE)
    )

