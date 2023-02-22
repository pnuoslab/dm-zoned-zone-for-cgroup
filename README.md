# dm-zoned-zone-for-cgroup
using cgroup isolated zone each other

# path
/usr/src/{kernel ver}/driver/md/dm-zoned.h

/usr/src/{kernel ver}/driver/md/dm-zoned-metadata.c

/usr/src/{kernel ver}/driver/md/dm-zoned-target.c


# how to use
It working based dm-zoned-cgroup

so if isolate the zone per cgroup 

Need to set the weight as same as the dm-zoned-cgroup(please set same value each Cgroup)
