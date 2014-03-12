import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
import jump_interface
import minimize_cost
from IPython import embed

DB_PATH = "/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db"
NUM_ITERATIONS = 1000

RAT_CLUSTERS = {'V1':[3,6], 'V2':[3,6], 'V3':[3,6], 'V4':[1,2,3,4,5,6,7,8],
        'V5':[3,6], 'V6':[1,2,3,4,5,6,7,8], 'V15':[1,3,6,8], 'V17':[3,6], 
        'V18':[3,6], 'V31':[1,3,6,8], 'V32': [3,6]}

def resample(data):
    idx = np.random.randint(0, len(data), len(data))

    return data[idx,:]

def confidence_interval(measurements, confidence):
    num_bins = int(np.sqrt(len(measurements)))
    plt.hist(measurements, bins=num_bins)
    density, edges = np.histogram(measurements, bins=num_bins)
    density = density / float(sum(density))
    distribution = np.cumsum(density)

    a = sum(distribution < confidence)
    a = edges[a]
    b = len(edges) - sum(distribution > (1-confidence))
    b = edges[b]

    return [a,b]

def main(data, estimator):
    if not hasattr(estimator, '__call__'):
        raise Exception("Variable estimator must be a function.")
    if not isinstance(data, np.ndarray):
        raise Exception("Variable data must be a numpy array.")     

    measurements = np.zeros((NUM_ITERATIONS, data.shape[1]))
    error_count = 0
    
    for n in range(NUM_ITERATIONS):
        if n%10 == 0: print(n)
        try:
            resampled_data = resample(data)
            measurements[n,:] = estimator(resampled_data, RAT_CLUSTERS[rat])
        except:
            error_count += 1

    confidence = [confidence_interval(m, 0.05) for m in measurements.T]

    print('error_count=' + str(error_count))
    print('theta_mean=' + str(np.mean(measurements[:,0])))
    print('theta_confidence=' + str(confidence[0]))
    print('b_mean=' + str(np.mean(measurements[:,1])))
    print('b_confidence=' + str(confidence[1]))

    return measurements, confidence

if __name__ == '__main__':
    ji = jump_interface.JumpInterface(DB_PATH)
    
    rats = ['V5', 'V15', 'V17', 'V18', 'V31', 'V32']
    for rat in rats:
        jumps = ji.get_jumps(rat)[:,0:2]
        measurements, confidence = main(jumps, minimize_cost.main)
        np.save(rat+'_measurements', measurements)





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