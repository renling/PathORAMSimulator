
#ifndef PARAMETERS_ORAM_H
#define PARAMETERS_ORAM_H

//#define USE_ORAM							false
//#define USE_UNIFIED_ORAM					false
#define MODEL_STASH							false
#define MODEL_DRAM							false
#define USE_DRAMSIM2						false

// inclusive vs excluse ORAM
#define USE_EXCLUSIVE_ORAM				 	true
#define CACHE_LEAF_ON_CHIP					true

// ORAM interval
#define USE_ORAM_INTERVAL					false
#define O_INT 								100

// capacity, block size and Z
//#define WORKING_SET_SIZE 					(4294967296)
//#define ORAM_SIZE							(8589934592)
//#define FINAL_POS_MAP_SIZE					204800			// 200KB
//#define ORAM_Z								4
//#define ORAM_BLOCK_SIZE						128
// ORAM stash
#define ORAM_LCSZ							200

// DRAM model
#define DRAM_BANDWIDTH						20				// GB/s

// AES and hash
#define AES_SEED_SIZE						8				// Bytes
#define AES_TYPE							128				// AES-128
#define AES_LATENCY							50
#define INTEGRITY							false
#define HASH_LENGTH_IN_BITS					128				// We can always truncate it to 128-bits
#define HASH_LATENCY						20

// PLB
//#define USE_PLB								false
//#define PLB_CAPACITY						(8192)
//#define PLB_WAYS							4

#define UNIFIED_PLB							true
//#define COMPRESSED_POSMAP					false
#define OUTER_COUNTER_WIDTH					64				// in bits
#define INNER_COUNTER_WIDTH					14				// in bits

// early forwarding data
#define ORAM_HALF_LATENCY_TRICK				true

// old HORAM block size and Z
#define POSITION_MAP_Z						oram_z
#define POSITION_MAP_BLOCK_SIZE				32

#endif
