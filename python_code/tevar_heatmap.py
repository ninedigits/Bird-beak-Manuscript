import math
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pylab as pylab


def center_of_mass(aorta_df):
    centroids = aorta_df.groupby('n')[['x','y','z']].mean()
    centroid_z1 = centroids.iloc[0]['z']
    centroids_crop = centroids[centroids.z >= centroid_z1]
    return centroids_crop



def tevar_individual_fig(pat_id,
                  file_name_root,
                  aortic_contours_fn,
                  arc_curve_diam_cd_fn,
                  aortic_innerline_fn,
                  plp_fn,
                  min_cd,
                  max_cd,
                  bbh_data,
                  pat_index,
                  fig,
                  arrow_direction_vector
                  ):
    aortic_contours, \
    aortic_innerline, \
    curvature_data, \
    plp, \
    bbh = load_and_process(aortic_contours_fn,
                           aortic_innerline_fn,
                           arc_curve_diam_cd_fn,
                           file_name_root,
                           max_cd,
                           min_cd,
                           pat_id,
                           plp_fn,
                           bbh_data,
                           pat_index)

    aortic_innerline, contours_map, index_plp, mapped_centroid_innerpoints, plp = curvature_colormapper(aortic_contours,
                                                                                                        aortic_innerline,
                                                                                                        curvature_data,
                                                                                                        pat_index, plp)

    #fig.subplots_adjust(bottom=0.5)
    fontsize = 20

    ax = fig.add_subplot(1, 1, 1, projection='3d', frame_on=False)

    pad = 1
    ax.tick_params(axis='x', pad=pad)
    ax.tick_params(axis='y', pad=pad)
    ax.tick_params(axis='z', pad=pad)

    ax.set_xlabel('x (mm)', fontsize=fontsize)
    ax.set_ylabel('y (mm)', fontsize=fontsize)
    ax.set_zlabel('z (mm)', fontsize=fontsize)
    ax.set_aspect('equal')
    
    ax.set_facecolor('white')
    
    gs1 = gridspec.GridSpec(4, 4)
    gs1.update(wspace=0.025, hspace=0.05) # set the spacing between axes. 
    
    R, G, B, A = [0, 0, 0, 0]
    
    ax.w_xaxis.set_pane_color((R, G, B, A))
    ax.w_yaxis.set_pane_color((R, G, B, A))
    ax.w_zaxis.set_pane_color((R, G, B, A))
    tick_size = 15
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(tick_size)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(tick_size)
    for tick in ax.zaxis.get_major_ticks():
        tick.label.set_fontsize(tick_size)

    # Contours
    ax.scatter3D(aortic_contours.x,
                 aortic_contours.y,
                 aortic_contours.z, marker='.', s=20,
                 c=contours_map.cd,
                 cmap='jet',
                 vmin=min_cd,
                 vmax=max_cd)
    # Innerline

    # Mapped Centroid Innerline Reference
    ax.scatter3D(mapped_centroid_innerpoints.x,
                 mapped_centroid_innerpoints.y,
                 mapped_centroid_innerpoints.z, marker='*', s=12)
    # Endograft
    # ax.scatter3D(contours_graft_pat_i.x,
    #              contours_graft_pat_i.y,
    #              contours_graft_pat_i.z, c='k', s=1) #, marker='.', s=0.5, c='0.5')
    # Landing Point
    ax.scatter3D(plp.x,
                 plp.y,
                 plp.z, c='k', marker='*', s=5)

    cmap = ax.scatter3D(aortic_innerline.x,
                        aortic_innerline.y,
                        aortic_innerline.z, marker='o', s=6,
                        c=aortic_innerline.cd,
                        cmap='jet',
                        vmin=min_cd,
                        vmax=max_cd)

    u, w, v = arrow_direction_vector
    ax.quiver(plp.x, plp.y, plp.z, u, w, v, arrow_length_ratio=.3, length=30, color='r', pivot='tip')
    ax.view_init(elev=10., azim=150)
    #ax.set_title("Patient P{0} Aorta \nBBH {1} mm".format(pat_id, np.round(bbh, 2)), fontsize=fontsize*1.5, y=.9)
    return cmap



def tevar_fig_gen(pat_id,
                  file_name_root,
                  aortic_contours_fn,
                  arc_curve_diam_cd_fn,
                  aortic_innerline_fn,
                  plp_fn,
                  min_cd,
                  max_cd,
                  bbh_data,
                  pat_index,
                  subplotnum,
                  fig,
                  subplot_dims, arrow_direction_vector,
                  gs):

    aortic_contours, \
    aortic_innerline, \
    curvature_data, \
    plp, \
    bbh = load_and_process(aortic_contours_fn,
                           aortic_innerline_fn,
                           arc_curve_diam_cd_fn,
                           file_name_root,
                           max_cd,
                           min_cd,
                           pat_id,
                           plp_fn,
                           bbh_data,
                           pat_index)

    aortic_innerline, contours_map, index_plp, mapped_centroid_innerpoints, plp = curvature_colormapper(aortic_contours,
                                                                                                        aortic_innerline,
                                                                                                        curvature_data,
                                                                                                        pat_index, plp)
    m, n = subplot_dims

    ax = fig.add_subplot(gs[subplotnum-1], projection='3d', frame_on=False)

    
        # Make background white
    ax.set_facecolor([.5,.5,.5])

    # Makes grid background white
    bgrnd = [.5 for i in [1, 1, 1, 1]]
    ax.w_xaxis.set_pane_color(bgrnd)
    ax.w_yaxis.set_pane_color(bgrnd)
    ax.w_zaxis.set_pane_color(bgrnd)
    # Hide grid lines
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    
    gs1 = gridspec.GridSpec(4, 4)
    gs1.update(wspace=0.025, hspace=0.05) # set the spacing between axes. 
        # Contours
    ax.scatter3D(aortic_contours.x,
                 aortic_contours.y,
                 aortic_contours.z, marker='.', s=3,
                 c=contours_map.cd,
                 cmap='jet',
                 vmin=min_cd,
                 vmax=max_cd)
    # Innerline

    # Mapped Centroid Innerline Reference
    ax.scatter3D(mapped_centroid_innerpoints.x,
                 mapped_centroid_innerpoints.y,
                 mapped_centroid_innerpoints.z, marker='*', s=12)
    # Endograft
    # ax.scatter3D(contours_graft_pat_i.x,
    #              contours_graft_pat_i.y,
    #              contours_graft_pat_i.z, c='k', s=1) #, marker='.', s=0.5, c='0.5')
    # Landing Point
    ax.scatter3D(plp.x,
                 plp.y,
                 plp.z, c='k', marker='*', s=5)


    cmap = ax.scatter3D(aortic_innerline.x,
                 aortic_innerline.y,
                 aortic_innerline.z, marker='o', s=6,
                 c=aortic_innerline.cd,
                 cmap='jet',
                 vmin=min_cd,
                 vmax=max_cd)

    u, w, v = arrow_direction_vector
    ax.quiver(plp.x, plp.y, plp.z, u, w, v, arrow_length_ratio=.3, length=30, color='r', pivot='tip')
    #ax.view_init(elev=10., azim=150)
    ax.view_init(elev=0., azim=170)
    ax.margins(x=-0.4, y=-0.4)
    ax.set_xlim([aortic_contours.x.min(),aortic_contours.x.max()])
    ax.set_ylim([aortic_contours.y.min(),aortic_contours.y.max()])
    ax.set_zlim([aortic_contours.z.min(),aortic_contours.z.max()])
    ax.set_aspect('equal')
    ax.axis('off')
    #plt.subplots_adjust(wspace=None, hspace=None)
    #ax.set_title("Patient T{0} Aorta \nBBH {1} mm".format(pat_id, np.round(bbh, 2)))
    return cmap



def tevar_heatmap(pat_id,
                  file_name_root,
                  aortic_contours_fn,
                  arc_curve_diam_cd_fn,
                  aortic_innerline_fn,
                  plp_fn,
                  min_cd,
                  max_cd,
                  bbh_data,
                  pat_index,
                  ):

    """
    Takes target patient and file structures/paths.
    Returns dataframes of aortic geometry.

    :param pat_id:
    :param file_name_root:
    :param aortic_contours_fn:
    :param arc_curve_diam_cd_fn:
    :param aortic_innerline_fn:
    :param plp_fn:
    :param min_cd:
    :param max_cd:
    :param bbh_data:
    :param pat_index:
    :return:

    """

    aortic_contours, \
    aortic_innerline, \
    curvature_data, \
    plp, \
    bbh = load_and_process(aortic_contours_fn,
                           aortic_innerline_fn,
                           arc_curve_diam_cd_fn,
                           file_name_root,
                           max_cd,
                           min_cd,
                           pat_id,
                           plp_fn,
                           bbh_data,
                           pat_index)

    aortic_innerline, contours_map, index_plp, mapped_centroid_innerpoints, plp = curvature_colormapper(aortic_contours,
                                                                                                        aortic_innerline,
                                                                                                        curvature_data,
                                                                                                        pat_index, plp)

    aorta_and_curvature_subplots(aortic_contours, aortic_innerline, bbh, contours_map, curvature_data, index_plp,
                                 mapped_centroid_innerpoints, max_cd, min_cd, pat_id, plp)

    # plt.savefig("figures/T{}_aortic_heatmap.png".format(pat_id), dpi=500)


def aorta_and_curvature_subplots(aortic_contours, aortic_innerline, bbh, contours_map, curvature_data, index_plp,
                                 mapped_centroid_innerpoints, max_cd, min_cd, pat_id, plp):
    # Subplot 1
    fig = plt.figure(figsize=(30, 7.5))
    plt.tight_layout()
    ax = fig.add_subplot(121, projection='3d', frame_on=False)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_zlabel('z (mm)')
    ax.set_aspect('equal')
    # Contours
    ax.scatter3D(aortic_contours.x,
                 aortic_contours.y,
                 aortic_contours.z, marker='.', s=2,
                 c=contours_map.cd,
                 cmap='jet',
                 vmin=min_cd,
                 vmax=max_cd)
    # Innerline
    ax.scatter3D(aortic_innerline.x,
                 aortic_innerline.y,
                 aortic_innerline.z, marker='o', s=6,
                 c=aortic_innerline.cd,
                 cmap='jet',
                 vmin=min_cd,
                 vmax=max_cd)
    # Mapped Centroid Innerline Reference
    ax.scatter3D(mapped_centroid_innerpoints.x,
                 mapped_centroid_innerpoints.y,
                 mapped_centroid_innerpoints.z, marker='*', s=12)
    # Endograft
    # ax.scatter3D(contours_graft_pat_i.x,
    #              contours_graft_pat_i.y,
    #              contours_graft_pat_i.z, c='k', s=1) #, marker='.', s=0.5, c='0.5')
    # Landing Point
    ax.scatter3D(plp.x,
                 plp.y,
                 plp.z, c='k', marker='*', s=5)
    ax.quiver(plp.x, plp.y, plp.z, 0, 1, 1, arrow_length_ratio=.3, length=20, color='r', pivot='tip')
    ax.view_init(elev=10., azim=330)
    ax.set_title("Patient T{0} Aorta \nBBH {1} mm".format(pat_id, np.round(bbh, 2)))
    # ax.set_facecolor('0.5')
    # Subplot 2
    ax = fig.add_subplot(122)
    ax.set_xlabel('arclength (mm)')
    ax.set_ylabel('curvature-diameter')
    # ax.set_facecolor('0.5')
    ax.scatter(curvature_data.i,
               curvature_data.cd,
               marker='.',
               c=curvature_data.cd,
               cmap='jet',
               vmin=min_cd,
               vmax=max_cd)
    index_plp = index_plp.i.values
    plt.plot([index_plp, index_plp], [curvature_data.cd.min(), curvature_data.cd.max()], color='red', linestyle='--')
    # plt.colorbar(curvature_data.norm)
    ax.set_title("Arclength as a Function of Curvature-Diameter")
    ax.grid()


def curvature_colormapper(aortic_contours, aortic_innerline, curvature_data, pat_index, plp):
    # Match curvature data to inner line
    aortic_innerline = aortic_innerline.merge(curvature_data, on='i', how='inner')
    pat_index += 1
    # Map each centroid to its respective location on the innerline
    mapped_centroid_list = []
    contour_centroids = aortic_contours.groupby(['i']).mean()
    contour_centroids.apply(lambda centroid_i: get_closest_point(centroid_i,
                                                                 aortic_innerline,
                                                                 mapped_centroid_list),
                            axis=1)
    mapped_centroid_innerpoints = pd.concat(mapped_centroid_list)
    aortic_innerline.drop_duplicates(subset=['x', 'y', 'z'], inplace=True)
    mapped_innerline_curve_to_centroid = aortic_innerline.merge(mapped_centroid_innerpoints,
                                                                on=['x', 'y', 'z'],
                                                                how='inner')
    contours_map = pd.merge(mapped_innerline_curve_to_centroid[['norm', 'index_centroid', 'cd']],
                            aortic_contours,
                            left_on='index_centroid',
                            right_on='i',
                            how='inner')
    # Map PLP to innerline curvature
    aortic_innerline = aortic_innerline.round(1)
    plp = plp.round(1)
    index_plp = aortic_innerline[aortic_innerline.apply(lambda x: xyz_coord_matcher(x, plp, 0.5), axis=1)]
    return aortic_innerline, contours_map, index_plp, mapped_centroid_innerpoints, plp


def load_and_process(aortic_contours_fn,
                     aortic_innerline_fn,
                     arc_curve_diam_cd_fn,
                     file_name_root,
                     max_colorscale_val,
                     min_colorscale_val,
                     pat_id, plp_fn, bbh_data, pat_index):
    """
    Reads and processes aortic contour and curvature data

    :param aortic_contours_fn:
    :param aortic_innerline_fn:
    :param arc_curve_diam_cd_fn:
    :param file_name_root:
    :param max_colorscale_val:
    :param min_colorscale_val:
    :param pat_id:
    :param plp_fn:
    :return:
    """
    patient_labels = [1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    patient_labels = ["T{}".format(i) for i in patient_labels]

    aortic_contours = pd.read_csv(file_name_root.format(pat_id, aortic_contours_fn), names=['i', 'x', 'y', 'z'])
    curvature_data = pd.read_csv(file_name_root.format(pat_id, arc_curve_diam_cd_fn),
                                 names=['i', 'curvature', 'diameter', 'cd'])
    aortic_innerline = pd.read_csv(file_name_root.format(pat_id, aortic_innerline_fn), names=['i', 'x', 'y', 'z'])
    # contours_graft_pat_i = pd.read_csv(file_name_root.format(pat_id, graft_fn), names=['i', 'x', 'y', 'z'])
    plp = pd.read_csv(file_name_root.format(pat_id, plp_fn), names=['x', 'y', 'z'])
    # Scaling data
    curvature_data.cd = curvature_data.curvature * curvature_data.diameter
    # max_colorscale_val and min_colorscale_val are the respective min, max of the total population
    curvature_data['norm'] = (curvature_data.cd - min_colorscale_val) / \
                             (max_colorscale_val - min_colorscale_val)
    bbh_data['pat_id'] = pd.Series(patient_labels)
    bbh = bbh_data[bbh_data['pat_id'].values == 'T{}'.format(pat_id)].values[0][0]
    return aortic_contours, aortic_innerline, curvature_data, plp, bbh


def get_closest_point(centroid, innerline, mapped_centroid_list):
    """Map centroids of contours with the innerline indices

    centroid: point
    innerline: many points

    """
    # distance formula
    index = centroid.name
    dx = (centroid.x - innerline.x) ** 2
    dy = (centroid.y - innerline.y) ** 2
    dz = (centroid.z - innerline.z) ** 2

    d = dx + dy + dz
    d = d.apply(lambda x: math.sqrt(x))
    minimum_distance = d.min()
    i_inner = pd.DataFrame({'index_centroid': index,
                            'd': d.values,
                            'x': innerline.x.values,
                            'y': innerline.y.values,
                            'z': innerline.z.values})

    min_coords = i_inner[['index_centroid', 'x', 'y', 'z']][i_inner.d == minimum_distance]
    min_coords = min_coords.iloc[[0]]
    mapped_centroid_list.append(min_coords)


def xyz_coord_matcher(row, target, tolerance):
    # Accounts for rounding errors when matching up cartesian coordinates
    all_true = []
    x_val = row['x']
    y_val = row['y']
    z_val = row['z']
    u_val = target['x'].values
    w_val = target['y'].values
    v_val = target['z'].values
    all_true.append(x_val + tolerance >= u_val >= x_val - tolerance)
    all_true.append(y_val + tolerance >= w_val >= y_val - tolerance)
    all_true.append(z_val + tolerance >= v_val >= z_val - tolerance)
    return all(all_true)


def interpolate(a, b, n):
    """
    Given points a, b, create a line with n equally spaced points between a and b
    :param a: (df)
    :param b: (df)
    :param n: (int)
    :return: (df)
    """
    line_output = []
    ab_dist = math.sqrt((a.x - b.x)**2 + (a.y - b.y)**2 + (a.z - b.z)**2)
    ab_diff = pd.DataFrame(columns=['x', 'y', 'z'])
    ab_diff.x, ab_diff.y, ab_diff.z = b.x - a.x, b.y - a.y, b.z - a.z
    v = ab_diff/ab_dist
    iters = np.linspace(0, ab_dist, n)
    for i in iters:
        increment = i * v
        increment.x, increment.y, increment.z = increment.x + a.x, increment.y + a.y, increment.z + a.z
        line_output.append(increment)
    line_output = pd.concat(line_output)
    return line_output


    return ab_diff
