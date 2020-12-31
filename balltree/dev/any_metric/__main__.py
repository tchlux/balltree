
# Function for printing the help message.
def print_help_message():
    print("""

balltree -- A very fast python-wrapped fortran implementation of a
    `ball tree`, includes fast `argsort` and `argselect` codes for
    flat arrays as well.

### Python

```python
import numpy as np
import balltree

points = np.random.random(size=(100000,1000))
tree = balltree.BallTree(points, [leaf_size=<leaf size int>])
distances, indices = tree.query(points[:1,:], [k=<k nearest int>])
# 'distances' will have the distances to the `k` nearest points to points[0]
# 'indices' will have the index of the `k` nearest points to points[0]


x = np.random.random(size=(10000,))
x, i = balltree.argsort(x, i=np.arange(len(x)))
# `x` is sorted in place
# `i` is sorted symmetrically with `x` in place, in this case it
#     shows the original indices of the elements of `x` (before sort).


x = np.random.random(size=(10000,))
k = 0
x, i = balltree.argselect(x, k, i=np.arange(len(x)))
# `x` is rearranged such that the element at `x[k]` has rank `k`
# `i` is rearranged symmetrically with `x` in place, in this case it
#     shows the original indices of the elements of `x` (before sort).
```

    """)

# If no arguments were provided, give the usage.
print_help_message()

