import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from main import *
'''
P_array = df_P_values.iloc[0, :].values #!!!!!!!!!!!!!!!!!!!!!!!
P_matrix = P_array
# print(P_matrix)   #[debugging]
# print(df_task)    #[debugging]
'''
'''for i, j in zip(range(10), range(int(start), int(end)+1)):
    P_matrix[i] = df_task.at[i, 'left bound'] + ( df_task.at[i, 'right bound'] - df_task.at[i, 'left bound'] )* (df_P_values.at[j, f"t{i+1}"])
    print(P_matrix[j])'''



# print(point_A_coords, point_D_coords, angle0, AB_distance, CB_distance, CD_distance, CG_distance, DE_distance)     #[debugging]


fig, ax = plt.subplots()
ax.set_xlim(-150, +60)
ax.set_ylim(-150, +60)
line_AB, = ax.plot([], [], 'k-', lw=2)  # Line AB
line_DC, = ax.plot([], [], 'k-', lw=2)  # Line DC
line_BC, = ax.plot([], [], 'k-', lw=2)  # Line BC
line_CG, = ax.plot([], [], 'k-', lw=2)  # Line for trajectory of G
line_CG_original, = ax.plot([], [], 'k-', lw=2)  # Original line for CG
line_EF, = ax.plot([], [], 'k-', lw=2)  # Line for EF
line_FG, = ax.plot([], [], 'k-', lw=2)  # Line for FG
line_FS, = ax.plot([], [], 'k-', lw=2)  # Line for FP
line_GS, = ax.plot([], [], 'k-', lw=2)  # Line for GP
point_A, = ax.plot(point_A_coords[0], point_A_coords[1], 'ko', markersize=8)  # Point A
point_B, = ax.plot([], [], 'ko', markersize=8)  # Point B
point_C, = ax.plot([], [], 'ko', markersize=8)  # Point C
point_D, = ax.plot(point_D_coords[0], point_D_coords[1], 'ko', markersize=8)  # Point D
point_G, = ax.plot([], [], 'ko', markersize=8)  # Point G
point_E, = ax.plot([], [], 'ko', markersize=8)  # Point E
point_F, = ax.plot([], [], 'ko', markersize=8)  # Point F


