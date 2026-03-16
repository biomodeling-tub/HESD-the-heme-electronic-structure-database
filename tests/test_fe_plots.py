#!/usr/bin/env python3
"""
Test script to generate just the Fe coordination plots to see if color coding works.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))

from plot_qm_analysis import QMAnalysisPlotter

def test_fe_plots():
    """Test Fe coordination plotting with color coding."""
    
    print("Testing Fe coordination plot generation...")
    plotter = QMAnalysisPlotter()
    
    # Generate just the Fe coordination analysis plots
    plotter.plot_fe_coordination_analysis(show_charge_mult=True, show_axial=True)
    
    print("Fe coordination plots generated!")
    print(f"Check output directory: {plotter.output_dir}")

if __name__ == "__main__":
    test_fe_plots()
