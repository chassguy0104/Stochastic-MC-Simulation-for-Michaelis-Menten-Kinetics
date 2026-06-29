import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.interpolate import interp1d
import os

plt.style.use('fivethirtyeight')
plt.figure(figsize=(12, 7), dpi=120)

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, 'MCsimulation_data.txt')

df = pd.read_csv(file_path, sep='\t\t', engine='python')

colors = {'S': '#1f77b4', 'E': '#d62728', 'ES': '#9467bd', 'P': '#2ca02c'}
max_time = df['Time'].max()
common_time = np.linspace(0, max_time, 500)

all_S, all_E, all_ES, all_P = [], [], [], []

print("Plotting 100 trials and calculating averages...")
for trial_id, trial_data in df.groupby('Trial'):
    plt.step(trial_data['Time'], trial_data['S'], color=colors['S'], alpha=0.05, where='post')
    plt.step(trial_data['Time'], trial_data['E'], color=colors['E'], alpha=0.05, where='post')
    plt.step(trial_data['Time'], trial_data['ES'], color=colors['ES'], alpha=0.05, where='post')
    plt.step(trial_data['Time'], trial_data['P'], color=colors['P'], alpha=0.05, where='post')

    interp_S = interp1d(trial_data['Time'], trial_data['S'], kind='previous', bounds_error=False, fill_value='extrapolate')
    interp_E = interp1d(trial_data['Time'], trial_data['E'], kind='previous', bounds_error=False, fill_value='extrapolate')
    interp_ES = interp1d(trial_data['Time'], trial_data['ES'], kind='previous', bounds_error=False, fill_value='extrapolate')
    interp_P = interp1d(trial_data['Time'], trial_data['P'], kind='previous', bounds_error=False, fill_value='extrapolate')

    all_S.append(interp_S(common_time))
    all_E.append(interp_E(common_time))
    all_ES.append(interp_ES(common_time))
    all_P.append(interp_P(common_time))

import matplotlib.patheffects as pe
border = [pe.Stroke(linewidth=3, foreground='black'), pe.Normal()]

plt.plot(common_time, np.mean(all_S, axis=0), color=colors['S'], linewidth=1, path_effects=border)
plt.plot(common_time, np.mean(all_E, axis=0), color=colors['E'], linewidth=1, path_effects=border)
plt.plot(common_time, np.mean(all_ES, axis=0), color=colors['ES'], linewidth=1, path_effects=border)
plt.plot(common_time, np.mean(all_P, axis=0), color=colors['P'], linewidth=1, path_effects=border)

custom_legend = [
    Line2D([0], [0], color=colors['S'], lw=3, label='Substrate (S)'),
    Line2D([0], [0], color=colors['E'], lw=3, label='Enzyme (E)'),
    Line2D([0], [0], color=colors['ES'], lw=3, label='Complex (ES)'),
    Line2D([0], [0], color=colors['P'], lw=3, label='Product (P)')
]

plt.legend(handles=custom_legend, loc='center right', fontsize=12, frameon=True, shadow=True)

plt.title('Monte Carlo Ensemble with Mean Trajectory', fontsize=16, pad=15)
plt.xlabel('Time (Seconds)', fontsize=14)
plt.ylabel('Molecule Count', fontsize=14)
plt.tight_layout()

plt.show()