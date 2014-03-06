import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
import jump_interface
import minimize_cost
from IPython import embed

DB_PATH = "/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db"
NUM_ITERATIONS = 200

RAT_CLUSTERS = {'V1':[3,6], 'V2':[3,6], 'V3':[3,6], 'V4':[1,2,3,4,5,6,7,8],
        'V5':[3,6], 'V6':[1,2,3,4,5,6,7,8], 'V15':[1,3,6,8], 'V17':[3,6], 
        'V18':[3,6], 'V31':[1,3,6,8], 'V32': [3,6]}

def resample(data):
    idx = np.random.randint(0, len(data), len(data))

    return data[idx,:]

def find_std(measurements, NUM_ITERATIONS):
    u = np.mean(measurements)

    return np.sqrt(sum((measurements-u)**2)/(NUM_ITERATIONS-1))

def confidence_interval(measurements, confidence):
    fig1, ax1 = plt.subplots()
    bins, edges, _ = ax1.hist(measurements, bins=np.sqrt(len(measurements)), 
        normed=True, cumulative=True)
    fig2, ax2 = plt.subplots()
    ax2.hist(measurements, bins=np.sqrt(len(measurements)),
        normed=True, cumulative=False)

    a = sum(bins < confidence)
    a = edges[a]
    b = len(edges) - sum(bins > (1-confidence))
    b = edges[b]
    #embed()

    return [a,b]

def main(data, estimator):
    if not hasattr(estimator, '__call__'):
        raise Exception("estimator must be a function")    

    theta = []
    b = []
    error_count = 0
    
    for n in range(NUM_ITERATIONS):
        if n%10 == 0: print(n)
        try:
            resampled_data = resample(jumps)
            t, i = estimator(resampled_data, RAT_CLUSTERS[rat])
            theta.append(t)
            b.append(i)
        except:
            error_count += 1

    theta_confidence = confidence_interval(theta, 0.05)
    b_confidence = confidence_interval(b, 0.05)
    
    print('error_count=' + str(error_count))
    print('theta_mean=' + str(sum(theta)/len(theta)))
    print('theta_confidence=' + str(theta_confidence))
    print('b_mean=' + str(sum(b)/len(b)))
    print('b_confidence=' + str(b_confidence))


    plt.show()

if __name__ == '__main__':
    rat = 'V6'
    ji = jump_interface.JumpInterface(DB_PATH)
    jumps = ji.get_jumps(rat)[:,0:2]

    main(jumps, minimize_cost.main)





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
    x1_min = floor(np.amin(data[:,0]))
    x1_max = ceil(np.amax(data[:,0]))
    x2_min = floor(np.amin(data[:,1]))
    x2_max = ceil(np.amax(data[:,1]))

    x1_blocks = np.arange(x1_min, x1_max, BLOCK_SIZE)
    x2_blocks = np.arange(x2_min, x2_max, BLOCK_SIZE)

    num_rows = len(x1_blocks)
    num_cols = len(x2_blocks)
      x1_idx = np.random.randint(0, num_rows-1, NUM_BLOCKS)
    x2_idx = np.random.randint(0, num_cols-1, NUM_BLOCKS)
    
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
"""