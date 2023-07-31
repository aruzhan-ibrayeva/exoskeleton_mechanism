import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


file_path = r"LPtau-Generator.xlsx"
df_P_values = pd.read_excel(file_path, sheet_name='Sheet1')
df_task = pd.read_excel(file_path, sheet_name = 'Task')


P_array = df_P_values.iloc[0, :].values #!!!!!!!!!!!!!!!!!!!!!!!
P_matrix = P_array
# print(P_matrix)   #[debugging]
# print(df_task)    #[debugging]

start = df_task.at[0, 'start']  # the row 2 in excel = 0 in df_task
end = df_task.at[0, 'end']


N = 20
point_A_coords = (P_matrix[0], P_matrix[1])  
point_D_coords = (P_matrix[2], P_matrix[3])  
angle0 = P_matrix[4]
AB_distance = P_matrix[5]
CB_distance = P_matrix[6]
CD_distance = P_matrix[7]
CG_distance = P_matrix[8]
DE_distance = P_matrix[9]


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


