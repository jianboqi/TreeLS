library(TreeLS)

# open sample plot file
# file = system.file("extdata", "pine_plot.laz", package="TreeLS")
# file = 'D:/DatasetDir/ZCZW/MobileLiDAR/dual_scanners_sub.las'
file = 'D:/DatasetDir/ZCZW/广西数据/bls_uav_top.las'
tls = readTLS(file)

# normalize the point cloud
tls = tlsNormalize(tls, keep_ground = TRUE, force_ground_filtering=TRUE)
x = plot(tls) # colorPalette="White"

# extract the tree map from a thinned point cloud
thin = tlsSample(tls, smp.voxelize(0.01))
map = treeMap(thin, map.hough(min_density = 0.03), 0)
# x = plot(map, size=2)
add_treeMap(x, map, color='yellow', size=2)

# classify tree regions
tls = treePoints(thin, map, trp.autoseg());
add_treePoints(x, tls, size=2)
add_treeIDs(x, tls, cex = 2, col='yellow')

# classify stem points
tls = stemPoints(tls, stm.hough())
add_stemPoints(x, tls, color='red', size=4)

# make the plot's inventory
inv = tlsInventory(tls, d_method=shapeFit(shape='circle', algorithm = 'irls'))
# add_tlsInventory(x, inv)

# extract stem measures
# seg = stemSegmentation(tls, sgt.ransac.circle(n = 20))
# add_stemSegments(x, seg, color='white', fast=T)

# plot everything once
#x = tlsPlot(tls, map, inv, seg, fast=T)

# check out only one tree
#tlsPlot(tls, inv, seg, tree_id = 1)

