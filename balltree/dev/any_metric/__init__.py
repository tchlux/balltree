# Get the version number from the setup file
import os

DIRECTORY = os.path.dirname(os.path.abspath(__file__))
ABOUT_DIR = os.path.join(DIRECTORY, "about")
with open(os.path.join(ABOUT_DIR,"version.txt")) as f:
    __version__ = f.read().strip()

from .balltree import *

# A wrapper around an NMSlib k nearest neighbor search.
# Use `pip install --no-binary :all: nmslib` to build optimized install
# of necessary packages on the current machine. Source code is available
# [here](https://github.com/nmslib/nmslib).
class NMSlib:
    def __init__(self, data, M=15, efC=100, efS=100, num_threads=4, print_progress=False):
        import nmslib
        index_time_params = {'M': M, 'indexThreadQty': num_threads, 'efConstruction': efC}
        self.index = nmslib.init(method='hnsw', space='l2', data_type=nmslib.DataType.DENSE_VECTOR)
        self.index.addDataPointBatch(data)
        self.index.createIndex(index_time_params, print_progress=print_progress)
        self.index.setQueryTimeParams({'efSearch': efS})


    def query(self, points, k=1, return_distance=True, num_threads=4):
        # get all nearest neighbours for all the datapoint
        # using a pool of 4 threads to compute
        neighbours = self.index.knnQueryBatch(points, k=k, num_threads=num_threads)
        indices, distances = [l for l in zip(*neighbours)]
        # Convert indices and distances to NumPy arrays.
        indices = np.asarray(indices)
        distances = np.asarray(distances)
        # Return the results.
        if return_distance: return distances, indices
        else:               return indices

# A wrapper around an NGT tree (uses an indexed method to be faster).
# Use `pip3 install ngt` to install the package, source is available
# [here](https://github.com/yahoojapan/NGT).
class NGT:
    def __init__(self, data):
        import tempfile, ngtpy
        data = np.asarray(data)
        self.index_directory = tempfile.TemporaryDirectory()
        self.index_name = bytes(self.index_directory.name, 'ascii')
        # Can manually delete the index file with:
        #   tree = NGT( <data> )
        #   tree.index_directory.cleanup()
        ngtpy.create(self.index_name, data.shape[1])
        self.index = ngtpy.Index(self.index_name)
        self.index.batch_insert(data)
        self.index.save()

    def query(self, points, k=1, return_distance=True):
        # If only a single point was given, convert it to a matrix.
        if (len(points.shape) == 1): points = np.array([points])
        # `self.index.search` returns a list of the form:
        #   [(idx 1, distance 1), ..., (idx k, distance k)]
        indices = []
        distances = []
        for pt in points:
            output = self.index.search(pt, k)
            ids, dists = [row for row in zip(*output)]
            indices.append( ids )
            distances.append( dists )
        # Convert indices and distances to NumPy arrays.
        indices = np.array(indices)
        distances = np.array(distances)
        # Return the results.
        if return_distance: return distances, indices
        else:               return indices
