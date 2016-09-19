# PathORAMSimulator
A detailed simulator for Path ORAM stash overflow and DRAM latency model. Includes basic, recursive (hierarchical) and unified.

Usage: ./BinPathORam_test test_type arg_list

	-s : working set in MB

	-S : ORAM tree size in MB, incompatible with -L

	-L : levels, incompatible with -S

	-B : data ORAM block size in Bytes

	-b : PosMap ORAM block size in Bytes

	-Z : data ORAM blocks per bucket

	-z : PosMap ORAM blocks per bucket

	-m : max Stash size

	-w : simulate worst case?

	-r : simulation length

	-i : integrity verification

Parameter_ORAM.h contains more advanced configurations for power users.
