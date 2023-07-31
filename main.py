import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np
from variables import *

# ACCESSING GENERATOR

file_path = r"LPtau-Generator.xlsx"
df_P_values = pd.read_excel(file_path, sheet_name='generator')
df_task = pd.read_excel(file_path, sheet_name = 'Task')
dyad_break = 0
final_P_matrix = []

start = df_task.at[0, 'start']  # the row 2 in excel and the row 0 in df_task are the same
end = df_task.at[0, 'end']


for j in range(int(start),int(end)+1):

    P_matrix[0] = df_task.at[0, 'left bound'] + ((df_task.at[0, 'right bound'] - df_task.at[0, 'left bound'])* df_P_values.at[j, 't1'])
    P_matrix[1] = df_task.at[1, 'left bound'] + ((df_task.at[1, 'right bound'] - df_task.at[1, 'left bound'])* df_P_values.at[j, 't2'])
    P_matrix[2] = df_task.at[2, 'left bound'] + ((df_task.at[2, 'right bound'] - df_task.at[2, 'left bound'])* df_P_values.at[j, 't3'])
    P_matrix[3] = df_task.at[3, 'left bound'] + ((df_task.at[3, 'right bound'] - df_task.at[3, 'left bound'])* df_P_values.at[j, 't4'])
    P_matrix[4] = df_task.at[4, 'left bound'] + ((df_task.at[4, 'right bound'] - df_task.at[4, 'left bound'])* df_P_values.at[j, 't5'])
    P_matrix[5] = df_task.at[5, 'left bound'] + ((df_task.at[5, 'right bound'] - df_task.at[5, 'left bound'])* df_P_values.at[j, 't6'])
    P_matrix[6] = df_task.at[6, 'left bound'] + ((df_task.at[6, 'right bound'] - df_task.at[6, 'left bound'])* df_P_values.at[j, 't7'])
    P_matrix[7] = df_task.at[7, 'left bound'] + ((df_task.at[7, 'right bound'] - df_task.at[7, 'left bound'])* df_P_values.at[j, 't8'])
    P_matrix[8] = df_task.at[8, 'left bound'] + ((df_task.at[8, 'right bound'] - df_task.at[8, 'left bound'])* df_P_values.at[j, 't9'])
    P_matrix[9] = df_task.at[9, 'left bound'] + ((df_task.at[9, 'right bound'] - df_task.at[9, 'left bound'])* df_P_values.at[j, 't10'])

    point_A_coords = (P_matrix[0], P_matrix[1])  
    point_D_coords = (P_matrix[2], P_matrix[3])  
    angle0 = P_matrix[4]
    AB_distance = P_matrix[5]
    CB_distance = P_matrix[6]
    CD_distance = P_matrix[7]
    CG_distance = P_matrix[8]
    DE_distance = P_matrix[9]

    for i in range(1, N + 1):
        angle_AB_values.append((angle0 * math.pi / 180 + (i - 1) / (N - 1) * 2 * math.pi))

        xB_values.append((point_A_coords[0] + AB_distance * math.cos(angle_AB_values[i])))
        yB_values.append((point_A_coords[1] + AB_distance * math.sin(angle_AB_values[i])))

        DB_distance = math.sqrt((xB_values[i] - point_D_coords[0]) ** 2 + (yB_values[i] - point_D_coords[1]) ** 2)

        if i >= len(cos_psi_values):
            cos_psi_values.append(
                (CD_distance ** 2 + DB_distance ** 2 - CB_distance ** 2) / (2 * CD_distance * DB_distance)
            )

        if (abs(cos_psi_values[i]) <= 1):
            alpha = math.acos(cos_psi_values[i])
        else:
            dyad_break = 1  
            break

        if i >= len(cos_mu_values):
            cos_mu_values.append
            (
                (CD_distance ** 2 + CB_distance ** 2 - DB_distance ** 2) / (2 * CD_distance * CB_distance)
            )

        '''if (abs(cos_mu_values[i]) > 1):
            dyad_break = 1
            break'''

        angle_DC = math.atan2(yB_values[i] - point_D_coords[1], xB_values[i] - point_D_coords[0]) + alpha

        xC_values.append(point_D_coords[0] + CD_distance * math.cos(angle_DC))
        yC_values.append(point_D_coords[1] + CD_distance * math.sin(angle_DC))

        xG_values.append(xC_values[i] + (CG_distance / CB_distance) * (xB_values[i] - xC_values[i]))
        yG_values.append((CG_distance / CB_distance) * (yB_values[i] - yC_values[i]))

        xE_values.append(point_D_coords[0] + (DE_distance / CD_distance) * (xC_values[i] - point_D_coords[0]))
        yE_values.append(point_D_coords[1] + (DE_distance / CD_distance) * (yC_values[i] - point_D_coords[1]))

        xF_values.append(xG_values[i] + xE_values[i] - xC_values[i])
        yF_values.append(yG_values[i] + yE_values[i] - yC_values[i])

        beta_values.append(math.atan2(yF_values[i] - yG_values[i], xF_values[i] - xG_values[i]))

    if(dyad_break==1):
        print("dyad_break condition is not met; table was filled with 0's")
        # FILL THE TABLE WITH 0'S

    if(dyad_break==0):
        for i in range(1, N1+1):
            k += math.cos(beta_values[i])
            m += math.sin(beta_values[i])
            k_alpha += ((i-1)*math.cos(beta_values[i]))/(N1-1)
            k_beta += ((i-1)*math.sin(beta_values[i]))/(N1-1)
            b5 += (- xG_values[i]*math.cos(beta_values[i]) - yG_values[i]*math.sin(beta_values[i]))
            b6 += (xG_values[i]*math.sin(beta_values[i]) - yG_values[i]*math.cos(beta_values[i]))
            b7 += (-xG_values[i])
            b8 += (-yG_values[i])
            k_m += ((i-1)*xG_values[i]/(N1-1))

        a55 = -N1*(2*N1-1)/(6*(N1-1))
        A_matrix = np.array(   [[N1, 0, -k, -m, -k_alpha],
                                [0, N1, m, -k, k_beta],
                                [k, -m, -N1, 0, -(N1/2)],
                                [m, k, 0, -N1, 0],
                                [k_alpha, -k_beta, -(N1/2), 0, a55]])
        # print("A_matrix: ", A_matrix) # debugging

        B_matrix = np.array([b5, b6, b7, b8, -k_m])
        # print ("B matrix: ", B_matrix) # debugging

        X_solution = np.linalg.solve(A_matrix, B_matrix)
        print("Solution:")
        print(X_solution)

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

    xS_local = X_solution[0]
    yS_local = X_solution[1]
    point_S, = ax.plot([], [], 'ko', markersize=8)  # Point S
    trajectory_S, = ax.plot([], [], 'r--', lw=2)  # Trajectory of point S

    # Calculating S
    for i in range(1, N + 1):
        xS_values.append((xG_values[i] + xS_local * math.cos(beta_values[i])) - (yS_local * math.sin(beta_values[i])))
        yS_values.append((yG_values[i] + xS_local * math.sin(beta_values[i])) + (yS_local * math.cos(beta_values[i])))

    # table append
    for i in range(1, N+1):
        angle_AB = angle_AB_values[i]
        xB = xB_values[i]
        yB = yB_values[i]
        xC = xC_values[i]
        yC = yC_values[i]
        xG = xG_values[i]
        yG = yG_values[i]
        xE = xE_values[i]
        yE = yE_values[i]
        xF = xF_values[i]
        yF = yF_values[i]
        xS = xS_values[i]
        yS = yS_values[i]
        table.append((i, angle_AB, xB, yB, xC, yC, xG, yG, xE, yE, xF, yF, xS, yS))

    def update(frame):
        i, angle_AB, xB, yB, xC, yC, xG, yG, xE, yE, xF, yF, xS, yS = table[frame]
        line_AB.set_data([point_A_coords[0], xB], [point_A_coords[1], yB])
        line_DC.set_data([point_D_coords[0], xC], [point_D_coords[1], yC])
        line_BC.set_data([xB, xC], [yB, yC])
        line_CG.set_data([xC, xG], [yC, yG])
        line_CG_original.set_data([xC, xG], [yC, yG])
        line_EF.set_data([xE, xF], [yE, yF])
        line_FG.set_data([xF, xG], [yF, yG])
        line_FS.set_data([xF, xS], [yF, yS])
        line_GS.set_data([xG, xS], [yG, yS])
        point_B.set_data(xB, yB)
        point_C.set_data(xC, yC)
        point_G.set_data(xG, yG)
        point_E.set_data(xE, yE)
        point_F.set_data(xF, yF)
        point_S.set_data(xS, yS)
        trajectory_S.set_data([row[12] for row in table[:frame + 1]], [row[13] for row in table[:frame + 1]])
        return line_AB, line_DC, line_BC, line_CG, line_CG_original, line_EF, line_FG, line_FS, line_GS, point_B, point_C, point_G, point_E, point_F, point_S, trajectory_S

    # Saving the table as an excel file
    df = pd.DataFrame(table, columns = ["i", "φ", "xB", "yB", "xC", "yC", "xG", "yG", "xE", "yE", "xF", "yF", "xS", "yS"])
    df.to_excel("table_data.xlsx", index=False)

    ani = animation.FuncAnimation(fig, update, frames=N, interval=500, blit=True)
    plt.show()
    print(cos_psi_values)
    # end if