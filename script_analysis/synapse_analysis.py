import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import matplotlib.patches as mpatches
import plotly.graph_objects as go

"""
    Here are collection of methods helps with functional analysis of dendritic synapses
"""

def circular_mean(angles):
    """
        Compute average of two angles of orientation data [0 pi]
        problem is that since it's preferred orientation, we need to double the angle,
        otherwise the average of 0 and pi is going to be 2/pi (very wrong, since we assume 0 and pi as same direction)

        Doubles the angles to treat them as directions on a full circle,
        then halves the result to project back to [0, π).

    """
    angles = np.asarray(angles)
    doubled = 2 * angles
    mean_angle = np.arctan2(np.mean(np.sin(doubled)), np.mean(np.cos(doubled))) / 2

    return mean_angle % np.pi  # ensures it lands in [0, π)


def plot_skeleton_nodes_dist(sk_swc, distance_threshold):
    plt.figure(figsize=(8, 8))

    # Mask for nodes within threshold
    mask_near_soma = sk_swc['distance_from_soma'] < distance_threshold
    mask_far_soma = ~mask_near_soma  # Nodes outside threshold

    # Plot nodes based on distance
    plt.scatter(sk_swc.loc[mask_far_soma, 'x'], sk_swc.loc[mask_far_soma, 'y'],
                color='royalblue', label='All Nodes (Far)', s=1.5)
    plt.scatter(sk_swc.loc[mask_near_soma, 'x'], sk_swc.loc[mask_near_soma, 'y'],
                color='dimgray', label=f'Nodes < {distance_threshold}', s=0.1)

    # Plot soma
    mask_soma = sk_swc["parent"] == -1
    soma_df = sk_swc[mask_soma]
    plt.scatter(soma_df["x"], soma_df["y"], color="red", marker="^", s=100, label="Soma")

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'Soma <{distance_threshold} Nodes')
    # plt.legend()

    plt.axis('equal')
    plt.autoscale()
    plt.show()



def assignment_df_mapping(assignment_df, sk_dend_levels, mapping_component, new_column_name):
    """
        Here you basically map assignment_df and add values from sk_dend_levels
        Input:
        mapping_component - from dend_levels, "distance_from_soma"
        new_cloumn_name - name of column you are adding to assignment_df, "distance_from_soma"
    """

    map_component = sk_dend_levels.set_index('id')[mapping_component]
    assignment_df[new_column_name] = assignment_df['node_id'].map(map_component)

    return assignment_df