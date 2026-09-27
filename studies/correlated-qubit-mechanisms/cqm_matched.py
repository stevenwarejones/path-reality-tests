"""Matched gamma acquisition audit and conditional acceptance sensitivity.

These are deterministic sensitivity calculations, not fitted full-source models.
"""
from fractions import Fraction as F
import html
from pathlib import Path
import re
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
TABLE_URL = 'https://arxiv.org/html/2503.07354v1'
TABLE_SHA256 = '60d379ae05d606266eecc09ce87eebe4be45d60648a8b7afa806685aa1227dac'


def extract_table(page):
    match = re.search(r'<table id="Sx4.T2.2".*?</table>', page, re.S)
    if not match:
        raise ValueError('paper Table 2 not found')
    return match.group().encode()


def table_cells(table):
    output = []
    for row in re.findall(r'<tr\b.*?</tr>', table.decode(), re.S):
        cells = [re.sub(r'\s+', ' ', html.unescape(re.sub('<[^>]+>', '', cell))).strip()
                 for cell in re.findall(r'<t[dh]\b[^>]*>(.*?)</t[dh]>', row, re.S)]
        if len(cells) == 7 and cells[2].isdigit():
            output.append(cells[2:])
    if len(output) != 12:
        raise ValueError('unexpected Table 2 schema')
    return output


def plusminus(cell):
    m = re.fullmatch(r'(\d+\.\d+)\((\d+)\)', cell)
    if not m:
        raise ValueError('unexpected uncertainty notation')
    return F(m[1]), F(int(m[2]), 10**len(m[1].split('.')[1]))


def odd_signal_bracket(observed, background):
    """Signal contrast 2q when a signal flips ONE readout interval.

    b=P(background has any switch), independent of the signal; background can
    have arbitrary internal dependence. r=P(background is exactly e_j) lies
    in [0,b]. Then p=b+q(1-b-r). This function inverts that relation.
    It assumes faithful resolved readout, not the uncalibrated experimental HMM.
    """
    p, b = F(observed), F(background)
    if not 0 <= b < F(1, 2) or not b <= p <= (1+b)/2:
        raise ValueError('infeasible single-interval signal/background probabilities')
    return 2*(p-b)/(1-b), min(F(1), 2*(p-b)/(1-2*b))


def envelope_bracket(p, pe, b, be, multiplier=4):
    """Union over a deterministic rectangle, clipped to physical feasibility."""
    pl, ph = max(F(0), p-multiplier*pe), min(F(1), p+multiplier*pe)
    bl, bh = max(F(0), b-multiplier*be), min(F(1, 2), b+multiplier*be)
    if bh >= F(1, 2) or ph < bl or pl > (1+bh)/2:
        raise ValueError('invalid/empty diagnostic rectangle')
    # Minimum is at smallest p and largest b, unless q=0 is feasible.
    low = max(F(0), 2*(pl-bh)/(1-bh))
    # For p<=1/2 the inverse upper bound decreases in b. For p>1/2,
    # q=1/2 is feasible somewhere if the above feasibility checks pass.
    high = F(1) if ph >= F(1, 2) else min(F(1), 2*(ph-bl)/(1-2*bl))
    return low, high


def selection_floor(bare_low, copper_high):
    """Common retention rho needed for positive population contrast difference."""
    gap = F(bare_low)-F(copper_high)
    return None if gap <= 0 else 1/(1+gap)


def trace_audit(z):
    result = []
    for panel in 'ab':
        name = f'Figure04/Figure4{panel}_traces.csv'
        with z.open(name) as f:
            a = np.loadtxt(f, delimiter=',', skiprows=1)
        if a.shape != (21001, 6) or not np.isfinite(a).all():
            raise ValueError('unexpected trace schema')
        dt = float(np.median(np.diff(a[:, 0])))
        if not np.allclose(np.diff(a[:, 0]), dt, rtol=1e-7, atol=1e-11):
            raise ValueError('nonuniform figure time coordinate')
        streams = []
        for col, levels in [(3, [-2.5, -1.5]), (5, [-4., -3.])]:
            # Midlevels in these plotted digital traces mark masked intervals.
            valid = any_switch = odd_switch = multiple = 0
            for start in range(0, len(a)-100, 100):
                window = a[start:start+101, col]
                if not np.isin(window, levels).all():
                    continue
                changes = int(np.count_nonzero(np.diff(window)))
                valid += 1
                any_switch += changes > 0
                odd_switch += changes % 2
                multiple += changes > 1
            streams.append({'column_0based': col, 'valid_100_interval_blocks': valid,
                            'any_switch_blocks': any_switch, 'odd_switch_blocks': odd_switch,
                            'multiple_switch_blocks': multiple})
        result.append({'file': name, 'rows': len(a), 'displayed_span_s': float(a[-1, 0]-a[0, 0]),
                       'displayed_sample_period_s': dt, 'streams': streams})
    return result


def matched_analysis(zip_path):
    table = table_cells((zip_path.parent/'gamma-table2.html').read_bytes())
    records = []
    for index, cells in enumerate(table):
        counts, accepted = map(int, cells[:2])
        if not 0 <= counts <= accepted or accepted <= 0:
            raise ValueError('invalid footprint counts')
        p, pe = plusminus(cells[2])
        b, be = plusminus(cells[3])
        reported, reported_error = plusminus(cells[4])
        frequency = F(counts, accepted)
        if abs(frequency-p) > F(1, 200):
            raise ValueError('counts disagree with rounded observation frequency')
        bracket = odd_signal_bracket(frequency, b)
        envelope = envelope_bracket(frequency, pe, b, be)
        paper = 2*(frequency-b)/(1-2*b)
        if abs(paper-reported) > F(1, 200):
            raise ValueError('count/background reconstruction disagrees with reported contrast')
        records.append({'chip': 'nonCu' if index < 6 else '1umCu', 'qubit': index % 6+1,
            'counts': counts, 'accepted_windows': accepted,
            'exposure_hours': 11.78 if index < 6 else 5.35,
            'observed_frequency': float(frequency), 'reported_observed_error': float(pe),
            'background_linear_rate_times_window': float(b), 'reported_background_error': float(be),
            'published_contrast_reconstructed': float(paper),
            'published_contrast_rounded': float(reported), 'reported_contrast_error': float(reported_error),
            'conditional_contrast_bracket': [float(v) for v in bracket],
            'conditional_contrast_four_error_envelope': [float(v) for v in envelope],
            'exact_conditional_contrast_four_error_envelope': [str(v) for v in envelope],
            'same_interval_vs_separate_interval_contrast_difference': float(bracket[1]-bracket[0])})
    # Match the five NON-sensor distances in the two Figure 5c files.
    # This pairs distances; it does NOT establish common impact populations or response.
    pairs = []
    for bare, copper, distance in [(4, 3, 2.03), (1, 6, 4.62), (6, 1, 4.06),
                                   (3, 2, 5.36), (5, 4, 6.66)]:
        x, y = records[bare-1], records[6+copper-1]
        low = F(x['exact_conditional_contrast_four_error_envelope'][0])
        high = F(y['exact_conditional_contrast_four_error_envelope'][1])
        gap = low-high
        floor = selection_floor(low, high)
        pairs.append({'nonCu_qubit': bare, 'Cu_qubit': copper, 'approximate_distance_mm': distance,
                      'conditional_gap_lower_four_error_envelope': float(gap),
                      'exact_conditional_gap_lower': str(gap),
                      'equal_retention_floor_for_positive_population_gap': float(floor) if floor else None,
                      'exact_equal_retention_floor': str(floor) if floor else None})
    with zipfile.ZipFile(zip_path) as z:
        traces = trace_audit(z)
    return {'status': 'CONDITIONAL SENSITIVITY, not a full-record physical separation',
            'source': TABLE_URL+'#Sx4.T2', 'table_sha256': TABLE_SHA256,
            'rows': records, 'distance_pairs': pairs, 'trace_audit': traces,
            'assumptions': ['Four reported errors is a deterministic diagnostic envelope, not simultaneous coverage.',
                'Treat rate times window as background any-switch probability only for this sensitivity calculation.',
                'One signal interval, signal-independent background, faithful digitization and selection are additional assumptions.',
                'Retention must cover charge triggering AND parity masking relative to a specified common impact population.',
                'Distance matching does not establish identical kernels, trigger populations, or charge acceptance.']}
