import numpy as np


LARGE_TEST   = True
COMPARE_AGAINST_SKLEARN = True
TEST_TREE    = False
TEST_APPROX  = False
TEST_PRUNE   = False
TEST_SORT    = False


# × SUM(SQ) - 3.17 2.37 2.36 2.33 2.33
# ✓ OMP     - 2.51 2.31 2.30 2.26 2.27
# × PRE-ADD - 2.35 2.33 2.31 2.27 2.35


if COMPARE_AGAINST_SKLEARN:
    print()
    print("="*70)

    from util.system import Timer
    t = Timer()

    if LARGE_TEST: train, dim = 1000000, 100
    else:          train, dim = 7, 2
    test = 1
    leaf_size = 10
    k = 5
    print("Initializing data..", flush=True)
    np.random.seed(0)
    x = np.random.random(size=(train,dim))
    z = np.random.random(size=(test,dim))
    print()
    print("x:", x.shape)
    print("z:", z.shape)

    # ----------------------------------------------------------------
    from balltree import BallTree as BT
    print()
    print("Fortran Ball Tree")
    t.start()
    tree = BT(x, leaf_size=leaf_size)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d[0])
    print("i: ",i[0])
    d1, i1 = d[0].copy(), i[0].copy()
    # ----------------------------------------------------------------
    from sklearn.neighbors import BallTree
    print()
    print("Sklearn Ball Tree")
    t.start()
    tree = BallTree(x, leaf_size=leaf_size)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d[0])
    print("i: ",i[0])
    d2, i2 = d[0].copy(), i[0].copy()

    # ----------------------------------------------------------------
    print()
    print("Brute Force")
    # Convert to float64 (for regular test this does nothing, for integer...
    x = x.astype('float64')
    z = z.astype('float64')
    t.start()
    d = np.sqrt(np.sum(x**2, axis=1) + np.sum(z[0]**2) - 2 * np.dot(x, z[0]))
    i = np.argsort(d)
    qt = t.stop()
    print("Query time:", qt)
    i = i[:k]
    d = d[i]
    print("d: ",d)
    print("i: ",i)
    d3, i3 = d.copy(), i.copy()

    # ----------------------------------------------------------------
    max_diff = max(max(abs(d1 - d2)), max(abs(d1-d3)))
    ds_match = max_diff < 2**(-26)
    is_match = np.all(i1 == i2) and np.all(i1 == i3)
    print()
    print(f"Max difference in distance calculations:\n   {max_diff:.3e}")
    try: assert(ds_match and is_match)
    except:
        print()
        print("ERROR")
        print( "  is_match: ",is_match)
        print(f"  ds_match:  {ds_match} {max(abs(d1-d3)):.3e} {max(abs(d1 - d2)):.3e}")


if TEST_TREE:
    from balltree import BallTree
    np.random.seed(0)

    # if LARGE_TEST: size = (20000000,10)
    # if LARGE_TEST: size = (100000,1000)
    if LARGE_TEST: size = (5000000, 100)
    else:          size = (4, 2)
    print()
    print(f"Allocating array.. {size}", flush=True)

    x = np.random.random(size=size)

    # Handle printout based on test size.
    if len(x) < 20:
        print()
        print("x:",x.T.shape,"\n",x)
        print()

    print("Building tree..", flush=True)
    from util.system import Timer
    t = Timer()
    t.start()

    k = 3
    leaf_size = 1
    tree = BallTree(x, leaf_size=leaf_size, reorder=True)

    t.stop()
    print(f"done in {t()} seconds.", flush=True)

    if len(x) < 20:
        print()
        print('-'*70)
        print("Tree:")
        print(tree.tree.T)
        print()
        print("Index mapping:")
        print(tree.index_mapping)
        print()
        print("Order, Radius:")
        print(tree.order)
        print(tree.radii)
        print('-'*70)

    z = np.random.random(size=(1,size[1]))

    if not LARGE_TEST: print("\nz: ",z,"\n")

    
    t.start()
    d,i = tree.query(z, k=k)
    t.stop()
    print("ball tree query in",t(),"seconds")
    d,i = d[0], i[0]
    # i = tree.index_mapping[i]
    t.start()
    true_dists = np.linalg.norm(np.float64(x) - np.float64(z[0,:]), axis=1)
    t.stop()
    print("linear search in",t(),"seconds")
    print()
    print("Tree/Truth")
    print("i: \n",i,"\n",np.argsort(true_dists)[:k])
    print("d: \n",d,"\n",np.sort(true_dists)[:k])


if TEST_APPROX:
    from balltree import BallTree
    size = (15,1)
    x = np.linspace(0,size[0]-1,size[0]).reshape(size)
    print("x: ",x)
    # 
    print("Building tree..", flush=True)
    k = 5
    leaf_size = 1
    tree = BallTree(x, leaf_size=leaf_size)
    print("tree: ",tree)
    print("x: ",x)
    # 
    z = np.array([[1.0,]])
    d, i = tree.nearest(z)
    print("d: ",d)
    print("i: ",i)
    print("tree[i[0]]: ",tree[i[0,0]])
    print()
    print(tree.nearest(z, look_ahead=3))


if TEST_PRUNE:
    from balltree import BallTree, prune
    
    tree_size = 14
    level = 3
    indices = np.zeros(2**level, order='F', dtype=np.int64)
    indices, found = prune.level(tree_size, level, indices)
    assert(found == 7)
    assert(tuple(indices) == (4,5,7,8,11,12,14,0))


    # Function for verifying correctness.
    def small_diff(v1,v2): return abs(v1 - v2) < .01
    # Build a tree over random points in 1 dimension.
    np.random.seed(3)
    pts = np.random.random(size=(100000,1))
    tree = BallTree()
    tree.add(pts)
    tree.build()

    vals = sorted((tree[0], tree[2], tree[1],
                  tree[(tree.tree.shape[1] + 2) // 2 + 1],
                  tree[(tree.tree.shape[1] + 2) // 2]
    ))
    print()
    for v in vals: print(v)
    print()
    assert(all(small_diff(v1,v2) for (v1,v2) in zip(
        vals, np.linspace(0,1,5))))

    tree.prune(4, method="inner")
    print(tree)
    print(tree[:len(tree)][0,:])
    print(sorted(tree[:len(tree)][0,:]))
    print()
    assert(all(small_diff(v1,v2) for (v1,v2) in zip(
        sorted(tree[:len(tree)][0,:]), np.linspace(0,1,len(tree)))))

    tree.prune(3, method="inner")
    print(tree)
    print(tree[:len(tree)][0,:])
    print(sorted(tree[:len(tree)][0,:]))
    print()
    assert(all(small_diff(v1,v2) for (v1,v2) in zip(
        sorted(tree[:len(tree)][0,:]), np.linspace(0,1,len(tree)))))

    tree.prune(2, method="inner")
    print(tree)
    print(tree[:len(tree)][0,:])
    print(sorted(tree[:len(tree)][0,:]))
    print()
    assert(all(small_diff(v1,v2) for (v1,v2) in zip(
        sorted(tree[:len(tree)][0,:]), np.linspace(0,1,len(tree)))))

    tree.prune(1, method="inner")
    print(tree)
    print(tree[:len(tree)][0,:])
    print(sorted(tree[:len(tree)][0,:]))
    print()
    assert(small_diff(tree[0][0], .5))


if TEST_SORT:
    if LARGE_TEST: N = 1000000
    else:          N = 10
    from balltree import fast_sort
    # Import a timer.
    from util.system import Timer
    t = Timer()
    # Generate test numbers.
    print()
    print(f"Generating {N} numbers..", flush=True)
    x = np.random.random(size=N)
    i = np.arange(len(x)) + 1
    print()
    # Test the fortran code.
    pts = x.copy()
    ids = i.copy()
    t.start()
    pts, ids = fast_sort.argsort(pts, ids)
    t.stop()
    # Check for correctness.
    ids_match = np.all(x[ids-1] == pts)
    is_sorted = np.all(np.diff(pts)>=0)
    try: assert(ids_match and is_sorted)
    except:
        print("ERROR")
        print(" ids_match: ",ids_match)
        print(" is_sorted: ",is_sorted)
        print()
        print("pts:", pts)
        print("ids:", ids)
        print()
        print(x)
        print(x[ids-1])
        exit()
    print("argsort: %.6f"%(t.total))
    # Test the NumPy code.
    pts = x.copy()
    ids = i.copy()
    t.start()
    ids = pts.argsort()
    t.stop()
    print("numpy:   %.6f"%(t.total))
