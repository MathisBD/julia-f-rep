
# Benchmarking kd-voxelizer (menger_sponge(1) with a sphere cut out) :
[4, 8, 4, 4]     ==> 7002ms
[8, 4, 4, 4]     ==> 7489ms
[32, 4, 4]       ==> 7573ms
[4, 4, 8, 4]     ==> 7927ms
[16, 8, 4]       ==> 8090ms
[8, 16, 4]       ==> 8866ms
[4, 32, 4]       ==> 10466ms
[4, 16, 8]       ==> 13371ms
[16, 4, 8]       ==> 13945ms
[8, 8, 8]        ==> 13989ms
[4, 4, 4, 8]     ==> 14069ms
[4, 8, 16]       ==> 25299ms
[8, 4, 16]       ==> 26346ms
[32, 16]         ==> 26916ms
[4, 4, 32]       ==> 48322ms
[16, 32]         ==> 48422ms

# Benchmarking gpu naive voxelizer (menger_sponge(1) with a sphere cut out) :
[256] ==> 153ms