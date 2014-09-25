import FiberAllocation as FA

#load all the data
objects = FA.ObjectCatalog(filein="../data/desi/objects0.rdzipn")
fibers = FA.FiberSetup(filein="../data/desi/fiberpos.txt")
tiles = FA.TilingSetup(filein="../data/desi/plate_centers")

#define the healpix pixels to search for available galaxies
n_side = 256
objects.healpixelize(n_side)
tiles.healpixelize(n_side)

#find the available galaxies for the first 20 plates
for i in range(20):
    x, y = FA.find_available_galaxies(fibers, tiles, objects,i)
    print i
