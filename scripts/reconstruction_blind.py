"""
perform reconstruction on blind collapsed reads
input blind collapsed reads
output reconstrution result
always recon for anchors
"""


import os
import umap
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


# plot the distribution of umi each bead has
def blind_cnt_distribution(bead_all, bead_type):
    plt.figure(figsize=(8,6))
    sns.histplot(np.log10(bead_all['total_cnt']), bins=50)
    plt.xlabel('log10(total count)')
    plt.ylabel('number of '+bead_type)
    plt.title('blind '+bead_type+' total count distribution ({}), median={}'.format(len(bead_all), bead_all['total_cnt'].median()))
    plt.savefig(os.path.join(out_dir,bead_type+'_blind_cnt_distribution.png'),dpi=300)
    plt.close()


# plot the distribution of how many bc each bead covered
def blind_cover_bc_distribution(bead_cover_bc, bead_type):
    plt.figure(figsize=(8,6))
    sns.histplot(bead_cover_bc['cnt'], bins=50)
    plt.xlabel('bead covered')
    plt.ylabel('number of '+bead_type)
    plt.title(bead_type+' bead covered distribution ({}), median={}'.format(len(bead_cover_bc), bead_cover_bc['cnt'].median()))
    plt.savefig(os.path.join(out_dir,bead_type+'_blind_cover_bc_distribution.png'),dpi=300)
    plt.close()


# generate spase matrix from matching, with selection on anchor or target
def get_matrix(match_df,min_a_cnt,max_a_cnt, min_t_cnt,max_t_cnt, anchor, target):
    a_all = match_df.groupby(anchor)['cnt'].size().reset_index(name='bead_cnt')   
    a_sel = a_all.loc[(a_all['bead_cnt']>min_a_cnt) & (a_all['bead_cnt']<max_a_cnt),]
    t_all = match_df.groupby(target)['cnt'].size().reset_index(name='bead_cnt')  
    t_sel = t_all.loc[(t_all['bead_cnt']>min_t_cnt) & (t_all['bead_cnt']<max_t_cnt),]
    match_df = match_df[(match_df[anchor].isin(a_sel[anchor])) & (match_df[target].isin(t_sel[target]))]
    a_list = match_df.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt') 
    t_list = match_df.groupby(target)['cnt'].sum().reset_index(name='total_cnt') 
    print('a: {}'.format(len(a_list)))
    print('t: {}'.format(len(t_list)))
    a_dict = dict()
    t_dict = dict()
    for i in range(len(a_list)):
        a_dict[a_list.iloc[i,0]] = i
    for j in range(len(t_list)):
        t_dict[t_list.iloc[j,0]] = j
    a_coo = []
    t_coo = []
    [a_coo.append(a_dict[a]) for a in match_df[anchor]]
    [t_coo.append(t_dict[t]) for t in match_df[target]]
    counts_coo = sp.coo_matrix((match_df['cnt'], (a_coo, t_coo)))
    counts = counts_coo.tocsr()
    return counts, a_list, t_list

def ot_random_ref(n):
    RADIUS = 1500
    r_vals = np.sqrt(np.random.uniform(size=n)) * RADIUS
    thet_vals = np.random.uniform(size=n) * 2 * np.pi
    a_random = np.column_stack((r_vals * np.cos(thet_vals), r_vals * np.sin(thet_vals)))
    a_random = a_random.astype(float)
    a_random = pd.DataFrame(a_random)
    a_random.columns = ['xcoord','ycoord']
    return a_random


def get_args():
    parser = argparse.ArgumentParser(description='Process recon seq data.')
    parser.add_argument("-d", "--date",
        help="input experiment data.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s", "--sample",
        help="input sample id.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-a", "--anchor",
        help="define anchor bead.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-t", "--target",
        help="define target bead.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-e", "--exptype",
        help="define experiment type (seq or tags).",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c", "--core",
        help="define core type to use.",
        type=str,
        default='CPU',
    )
    args = parser.parse_args()
    return args



args = get_args()
sample = args.sample
date = args.date
anchor = args.anchor
target = args.target 
#sample_folder = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",date,sample)
sample_folder = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",date,sample)


# make output dir
out_dir = os.path.join(sample_folder,'recon')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

blind_raw = pd.read_csv(os.path.join(sample_folder, sample+'_blind_raw_reads_filtered.csv.gz'))
blind_sum = blind_raw.groupby(['R1_bc', 'R2_bc']).size().reset_index(name='cnt')
if args.exptype == 'seq':
    blind_sum.columns = [anchor, target, 'cnt'] #if seq
elif args.exptype == 'tags':
    blind_sum.columns = [target, anchor, 'cnt'] #if tags

# plot blind cnt distribution
a_all = blind_sum.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt')
t_all = blind_sum.groupby(target)['cnt'].sum().reset_index(name='total_cnt')

blind_cnt_distribution(a_all, anchor)
blind_cnt_distribution(t_all, target)

# plot bc covered
a_cover_bc = blind_sum.groupby(anchor).count()
t_cover_bc = blind_sum.groupby(target).count()

blind_cover_bc_distribution(a_cover_bc, anchor)
blind_cover_bc_distribution(t_cover_bc, target)


# generate matrix and reconstruction with CPU
if args.core == 'CPU':
    a_min = 0
    a_max = 1000
    t_min = 0
    t_max = 1000

    counts, a_sel, t_sel = get_matrix(blind_sum, min_a_cnt=a_min, max_a_cnt=a_max, min_t_cnt=t_min, max_t_cnt=t_max, anchor=anchor, target=target)

    reducer = umap.UMAP(metric='cosine',
                        n_neighbors=25, 
                        min_dist=0.99, 
                        low_memory=False, 
                        n_components=2, 
                        # random_state=0, 
                        verbose=True, 
                        n_epochs=50000,
                        # output_dens = True,
                        # local_connectivity = 30,
                        learning_rate = 1)
    embedding = reducer.fit_transform(np.log1p(counts))

    # output reconstruction result
    a_recon = pd.DataFrame(embedding)
    a_recon.columns = ['xcoord','ycoord']
    a_recon.insert(loc=0, column=anchor, value=a_sel[anchor])
    a_recon.to_csv(os.path.join(out_dir,'{}_recon_loc.csv'.format(anchor)), index=False)

    # plot the shape 
    plt.figure(figsize=(6,6))
    plt.scatter(embedding[:,0],embedding[:,1],s=1)
    plt.title(anchor+' umap ({})'.format(len(embedding)))
    plt.savefig(os.path.join(out_dir,'{}_UMAP.png'.format(anchor)),dpi=300)
    plt.close()
    # plt.show()

    # plot uniformness
    a_random = ot_random_ref(len(counts))
    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    im1 = sns.kdeplot(x=a_random.xcoord, y=a_random.ycoord, cmap="Blues", fill=True, levels=30, ax=axes[0], cbar = True)
    axes[0].set_title('random (max={:.7f})'.format(im1.collections[0].colorbar.vmax))
    im2 = sns.kdeplot(x=embedding[:,0],y=embedding[:,1], cmap="Blues", fill=True, levels=30, ax=axes[1], cbar = True)
    axes[1].set_title('recon (max={:.7f})'.format(im2.collections[0].colorbar.vmax))
    plt.savefig(os.path.join(out_dir,'{}_UMAP_density.png'.format(anchor)),dpi=300)
    plt.close()

    # plot covex
    hull = ConvexHull(embedding)
    area = hull.volume
    perimeter = np.sum(np.linalg.norm(embedding[hull.vertices] - np.roll(embedding[hull.vertices], -1, axis=0), axis=1))
    plt.figure(figsize=(6,6))
    plt.scatter(embedding[:,0], embedding[:,1], s=1, alpha=0.8)
    for simplex in hull.simplices:
        plt.plot(embedding[simplex, 0], embedding[simplex, 1], 'r-', linewidth=1.5, alpha = 0.7)
    plt.title(anchor+' umap convex hull ({:.5f})'.format(perimeter**2/area / (4*np.pi)), fontsize=16, pad=20)
    plt.savefig(os.path.join(out_dir,'{}_UMAP_convex.png'.format(anchor)),dpi=300)
    plt.close()


# generate matrix and reconstruction with GPU
elif args.core == 'GPU':
    from cuml.manifold.umap import UMAP as cuUMAP
    min_list = [0,50,100]
    min_dist_list = [0.3,0.99]
    n_epo_list = [100000]
    for m in min_list:
        for md in min_dist_list:
            for nepo in n_epo_list:
                out_gpu = os.path.join(out_dir,'recon_GPU_{}_{}_{}'.format(m,md,nepo))
                if not os.path.exists(out_gpu):
                    os.makedirs(out_gpu)

                a_min = m
                a_max = 1000000
                t_min = m
                t_max = 1000000

                counts, a_sel, t_sel = get_matrix(blind_sum, min_a_cnt=a_min, max_a_cnt=a_max, min_t_cnt=t_min, max_t_cnt=t_max, anchor=anchor, target=target)

                reducer = cuUMAP(metric='cosine',
                                    n_neighbors=25, 
                                    min_dist=md, 
                                    low_memory=False, 
                                    n_components=2, 
                                    # random_state=0, 
                                    verbose=True, 
                                    n_epochs=nepo,
                                    # output_dens = True,
                                    # local_connectivity = 30,
                                    learning_rate = 1)
                embedding = reducer.fit_transform(np.log1p(counts))

                # output reconstruction result
                a_recon = pd.DataFrame(embedding)
                a_recon.columns = ['xcoord','ycoord']
                a_recon.insert(loc=0, column=anchor, value=a_sel[anchor])
                a_recon.to_csv(os.path.join(out_gpu,'{}_recon_loc.csv'.format(anchor)), index=False)

                # plot the shape 
                plt.figure(figsize=(6,6))
                plt.scatter(embedding[:,0],embedding[:,1],s=1)
                plt.title(anchor+' umap ({})'.format(len(embedding)))
                plt.savefig(os.path.join(out_gpu,'{}_UMAP.png'.format(anchor)),dpi=300)
                plt.close()
                # plt.show()

                # plot uniformness
                a_random = ot_random_ref(len(counts))
                fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
                im1 = sns.kdeplot(x=a_random.xcoord, y=a_random.ycoord, cmap="Blues", fill=True, levels=30, ax=axes[0], cbar = True)
                axes[0].set_title('random (max={:.7f})'.format(im1.collections[0].colorbar.vmax))
                im2 = sns.kdeplot(x=embedding[:,0],y=embedding[:,1], cmap="Blues", fill=True, levels=30, ax=axes[1], cbar = True)
                axes[1].set_title('recon (max={:.7f})'.format(im2.collections[0].colorbar.vmax))
                plt.savefig(os.path.join(out_gpu,'{}_UMAP_density.png'.format(anchor)),dpi=300)
                plt.close()

                # plot covex
                hull = ConvexHull(embedding)
                area = hull.volume
                perimeter = np.sum(np.linalg.norm(embedding[hull.vertices] - np.roll(embedding[hull.vertices], -1, axis=0), axis=1))
                plt.figure(figsize=(6,6))
                plt.scatter(embedding[:,0], embedding[:,1], s=1, alpha=0.8)
                for simplex in hull.simplices:
                    plt.plot(embedding[simplex, 0], embedding[simplex, 1], 'r-', linewidth=1.5, alpha = 0.7)
                plt.title(anchor+' umap convex hull ({:.5f})'.format(perimeter**2/area / (4*np.pi)), fontsize=16, pad=20)
                plt.savefig(os.path.join(out_gpu,'{}_UMAP_convex.png'.format(anchor)),dpi=300)
                plt.close()




