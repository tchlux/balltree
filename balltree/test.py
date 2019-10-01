import numpy as np

TEST_SORT    = False
TREE_TEST    = False
TEST_PRUNE   = True
LARGE_TEST   = False
INTEGER_TEST = False
COMPARE_AGAINST_SKLEARN = False


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
    print("argsort:", t())
    # Test the NumPy code.
    pts = x.copy()
    ids = i.copy()
    t.start()
    ids = pts.argsort()
    t.stop()
    print("numpy:  ", t())


if TEST_PRUNE:
    from balltree import ball_tree, BallTree
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



if TREE_TEST:
    from balltree import BallTree
    np.random.seed(0)

    if LARGE_TEST: size = (100000,1000)
    else:          size = (4,2)
    print()
    print(f"Allocating array.. {size}", flush=True)

    if INTEGER_TEST:
        # Add 128 to convert uint8 into correctly ordered int8
        x = np.asarray((np.random.random(size=size)-1/2)*256, dtype='int8')
    else:
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

    k = 5
    leaf_size = 1
    tree = BallTree(x, leaf_size=leaf_size)

    t.stop()
    print(f"done in {t()} seconds.", flush=True)
    print()

    if len(x) < 20:
        print('-'*70)
        print("Tree:")
        print(tree.tree.T)
        print()
        print("Index mapping:")
        print(tree.index_mapping)
        print()
        print("Order, Radius, Squared Sums:")
        print(tree.order)
        print(tree.radii)
        print(tree.sq_sums)
        print('-'*70)

    if INTEGER_TEST:
        z = np.asarray((np.random.random(size=(1,size[1]))-1/2)*256,
                       dtype='int8')
    else:
        z = np.random.random(size=(1,size[1]))

    if not LARGE_TEST: print("\nz: ",z,"\n")

    d,i = tree.query(z, k=k)
    d,i = d[0], i[0]
    # i = tree.index_mapping[i]
    true_dists = np.linalg.norm(np.float64(x) - np.float64(z[0,:]), axis=1)
    print("Tree/Truth")
    print("i: \n",i,"\n",np.argsort(true_dists)[:k])
    print("d: \n",d,"\n",np.sort(true_dists)[:k])


if COMPARE_AGAINST_SKLEARN:
    print()
    print("="*70)

    from util.system import Timer
    t = Timer()

    if LARGE_TEST: train, dim = 100000, 1000
    else:          train, dim = 7, 2
    test = 1
    leaf_size = 10
    k = 5
    print("Initializing data..", flush=True)
    np.random.seed(0)
    if INTEGER_TEST:
        x = np.asarray((np.random.random(size=(train,dim))-1/2)*128, dtype='int8')
        z = np.asarray((np.random.random(size=(test, dim))-1/2)*128, dtype='int8')
    else:
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

