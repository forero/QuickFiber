import FiberAllocation as FA
import sys, string, os

min_index = int(sys.argv[1])
max_index = int(sys.argv[2])

print range(int(min_index), int(max_index))

#load all the data
objects = FA.ObjectCatalog(filein="/home/forero/work/data/MockDesiRandom/tmp_small_rand00", random=True)
fibers = FA.FiberSetup(filein="../data/desi/fiberpos.txt")
tiles = FA.TilingSetup(filein="../data/desi/plate_centers")

#define the healpix pixels to search for available galaxies
n_side = 256
objects.healpixelize(n_side)
tiles.healpixelize(n_side)

#find the available galaxies for the first 20 plates
for i in range(min_index, max_index):
    x, y = FA.find_available_galaxies(fibers, tiles, objects,i, filename="available_gals_small_rand00")
    print i
