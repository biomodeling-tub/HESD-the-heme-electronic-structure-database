import plots
import plot_qm_analysis
from scripts import plot_energy_vs_homo_lumo_gap


def test_charge_multiplicity_core_states_use_shared_colors():
    colors = plots.get_charge_multiplicity_colors(['0_1', '1_6', 'q=0,m=1', 'q=1,m=6'])

    assert colors['0_1'] == plots.COLORBLIND_PALETTE[0]
    assert colors['q=0,m=1'] == plots.COLORBLIND_PALETTE[0]
    assert colors['1_6'] == plots.COLORBLIND_PALETTE[3]
    assert colors['q=1,m=6'] == plots.COLORBLIND_PALETTE[3]
    assert plots.CHARGE_MULTIPLICITY_COLOR_MAP['1,5'] == plots.COLORBLIND_PALETTE[4]
    assert plots.CHARGE_MULTIPLICITY_COLOR_MAP['1,6'] == plots.COLORBLIND_PALETTE[3]


def test_charge_multiplicity_helpers_match_across_active_plot_modules():
    combos = ['0_1', '1_6', 'q=0,m=1', 'q=1,m=6']
    expected = plots.get_charge_multiplicity_colors(combos)

    assert plot_qm_analysis.get_charge_multiplicity_colors(combos) == expected
    assert plot_energy_vs_homo_lumo_gap.get_charge_multiplicity_colors(combos) == expected
    assert plot_qm_analysis.CHARGE_MULTIPLICITY_COLOR_MAP['1,6'] == expected['1_6']
