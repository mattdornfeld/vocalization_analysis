import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import jump_interface
import cluster_parameter
from IPython import embed

RAT_CLUSTERS = {'V1':[3,6], 'V2':[3,6], 'V3':[3,6], 'V4':[1,2,3,4,5,6,7,8],
        'V5':[3,6], 'V6':[1,2,3,4,5,6,7,8], 'V15':[1,3,6,8], 'V17':[3,6], 
        'V18':[3,6], 'V31':[1,3,6,8], 'V32': [3,6]}

def resample(data, block_size, num_blocks):
    x1_min = floor(np.amin(data[:,0]))
    x1_max = ceil(np.amax(data[:,0]))
    x2_min = floor(np.amin(data[:,1]))
    x2_max = ceil(np.amax(data[:,1]))

    x1_blocks = np.arange(x1_min, x1_max, block_size)
    x2_blocks = np.arange(x2_min, x2_max, block_size)
    
    num_rows = len(x1_blocks)
    num_cols = len(x2_blocks)

    x1_idx = np.random.randint(0, num_rows-1, num_blocks)
    x2_idx = np.random.randint(0, num_cols-1, num_blocks)
    
    verts = []
    codes = []

    for i, j in zip(x1_idx, x2_idx):
    	verts += [(x1_blocks[i], x2_blocks[j])]
    	verts += [(x1_blocks[i], x2_blocks[j+1])]
    	verts += [(x1_blocks[i+1], x2_blocks[j+1])]
    	verts += [(x1_blocks[i+1], x2_blocks[j])]
    	verts += [(x1_blocks[i], x2_blocks[j])]
    	codes += [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO]

    poly = Path(verts, codes)
    
    return data[poly.contains_points(data)]

def find_std(params, num_iterations):
    u = np.mean(params)

    return np.sqrt(sum((params-u)**2)/(num_iterations-1))

def confidence_interval(params, confidence):
    fig1, ax1 = plt.subplots()
    bins, edges, _ = ax1.hist(params, bins=np.sqrt(len(params)), 
        normed=True, cumulative=True)
    fig2, ax2 = plt.subplots()
    ax2.hist(params, bins=np.sqrt(len(params)),
        normed=True, cumulative=False)


    a = sum(bins < confidence)
    a = edges[a]
    b = len(edges) - sum(bins > (1-confidence))
    b = edges[b]
    #embed()

    return [a,b]

def main(rat, estimator):
    if not hasattr(estimator, '__call__'):
        raise Exception("estimator must be a function")
    if not isinstance(rat, str):
        raise Exception("rat must be a string")
    dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
    ji = jump_interface.JumpInterface(dbPath)
    jumps = ji.get_jumps(rat)[:,0:2]
    block_size = 300
    num_blocks = 10000
    num_iterations = 100

    theta = []
    intercept = []
    for n in range(num_iterations):
        if n%10 == 0: print(n)
        resampled_data = resample(jumps, block_size, num_blocks)
        t, i = estimator(resampled_data, RAT_CLUSTERS[rat])
        theta.append(t)
        intercept.append(i)

    [a_theta, b_theta] = confidence_interval(theta, 0.05)
    [a_intercept, b_incercept] = confidence_interval(intercept, 0.05)
    #print(a_theta)
    #print(b_theta)
    #print(a_incercept)
    #print(b_intercept)
    plt.show()

if __name__ == '__main__':
    rat = 'V6'
    main(rat, cluster_parameter.main)





"""
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(jumps[:,0], jumps[:,1], 'o', markersize=2)
plt.show()
"""

"""
patch = patches.PathPatch(poly, facecolor='orange', lw=2)
ax.add_patch(patch)
ax.set_xlim(f1_min, f1_max)
ax.set_ylim(f2_min, f2_max)
plt.show()
"""