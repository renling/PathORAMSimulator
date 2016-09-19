
#ifndef HIER_BIN_PATH_ORAM
#define HIER_BIN_PATH_ORAM

#include "BinPathORam.h"
using namespace std;

struct ORAMInfoForDRAMSim
{
        int hierarchy;
        short AES_pad_in_bits;
        short integrity;
        float cpu_freq_factor;
        int blockSize[8];
        int blocksPerBucket[8];
        uint64_t blockCount[8];
        uint64_t bucketCount[8];
        int levelCount[8];
        uint64_t leafCount[8];
        int posMapScaleFactor[8];
        char * ini_dir; 
        char * deviceConfFileName; 
        char * systemConfFileName;
};

class HierBinPathORam
{

public:
	HierBinPathORam();
	~HierBinPathORam();
	int Configure(uint64_t workingSet, double * utilization, int * HORamBlockSize, int * HORamBlocksPerBucket, uint32_t maxPosMapSize);
	void DRAMConfigure(short AES_pad_in_bits, float CPU_freq_factor, char * ini_dir, char * deviceConfFileName, char * systemConfFileName);
	void SimulateLatency(int runlength);
	int Initialize();
	int SetPLB(int capacity, int blockSize, int ways);
	void SetMaxStashSize(int maxStashSize);
	int Access(int64_t addr, int hier_to_access, short RWoption);
	void GenerateAddr(int64_t addr);
	int ReadPLB();
	int UpdatePLB(int i);

	inline int GetHierarchy() {return Hierarchy;}
	inline int GetLevelCount(int i) {return LevelCount[i];}
	inline int GetBlockSize(int i) {return BlockSize[i];}
	inline int GetBlocksPerBucket(int i) {return BlocksPerBucket[i];}
	int GetBucketSizeInBits(int i);
	int GetBucketSizeInChunks(int i);
	inline int GetPosMapScaleFactor(int i) {return PosMapScaleFactor[i];}
	inline int GetHitDelay() {return HitDelay[Hierarchy - 1];}
	inline int GetHitDelay(int i) {return HitDelay[i];}
	inline int GetReadyDelay() {return ReadyDelay[Hierarchy - 1];}
	inline int GetReadyDelay(int i) {return ReadyDelay[i];}
	inline int GetFinalPositionMapSize() {return FinalPositionMapSize;}

	bool LocalCacheFull();
	bool LocalCacheAlmostFull();
	void BackgroundEvict();
	int BackgroundEvict_Count();

	uint64_t GetNumDummy() {return NumDummy;}
	uint64_t GetNumAccess() {return NumAccess;}
	uint32_t GetNumAccessInWindow() { return NumAccessInWindow;}
	uint32_t GetNumDummyInWindow() { return NumDummyInWindow;}
	void ResetNumDummyInWindow() { NumAccessInWindow = NumDummyInWindow = 0;}

	BinPathORam ** ORam;

#ifdef USE_DRAMSIM2
	friend class ORAMForDRAMSim;
	ORAMForDRAMSim * Simulator;
#endif
	
private:
	int Hierarchy;
	int MaxHierarchy;

	uint64_t * WorkingSet;
	double * Util;
	int * BlockSize;
	int * BlocksPerBucket;
	uint64_t * BlockCount;
	uint64_t * BucketCount;
	int * LevelCount;
	uint64_t * LeafCount;
	int * PosMapScaleFactor;
	uint64_t FinalPositionMapSize;

	short Integrity;
	float CPU_freq_factor;
	int * HitDelay;
	int * ReadyDelay;

	int64_t *Addr;
	SimpleCache ** PLB;

	uint32_t NumAccessInWindow;
	uint32_t NumDummyInWindow;
	uint64_t NumAccess;
	uint64_t NumDummy;

	char debug;
};

#ifdef USE_DRAMSIM2
class ORAMForDRAMSim
{
	friend class HierBinPathORam;
public: 
	ORAMForDRAMSim(ORAMInfoForDRAMSim info);
	~ORAMForDRAMSim();
	double GetAveHitLatency() {return SumHitLatency * 1.0 / NumAccess;}
	double GetAveReadyLatency() {return SumReadyLatency * 1.0 / NumAccess;}
	void Reset();

	void addRequest(uint64_t addr, bool rw);
	void tick() {
		cycles++;
		mem->update();
	}
	void drain() { while(issued_accesses != finished_reads + finished_writes) tick();}

	void read_complete(unsigned id, uint64_t address, uint64_t clock_cycle);
	void write_complete(unsigned id, uint64_t address, uint64_t clock_cycle);

	void SimulateLatency(int runlength, int hier);
	int SimulateOneAccess(int hier_to_access);
	void SimulateOneORAM(int h, bool rw, bool blocking);

	short Strategy;		// 0 = naive, 1 = subtree, 2 = load-balance, 3 = 1 + 2

private:
	ORAMInfoForDRAMSim Info;
	int Hierarchy;
	DRAMSim::MultiChannelMemorySystem *mem;
	short NumChannels;
	short NumColumns;
	short ChunkSize;
	int RowWidthInBytes;
	AddrShiftObj * AddrMapping;

	int * BlockSize;
	int * BlocksPerBucket;
	uint64_t * BlockCount;
	uint64_t * BucketCount;
	int * LevelCount;
	uint64_t * LeafCount;
	int * PosMapScaleFactor;

	uint64_t * BaseAddr;
	int * BucketSizeInBits;
	int * BucketSizeInChunks;
	short * KPacking;
	int64_t * LonelySubTrees;
	int * LastSubTreeSize;
	uint64_t * Path;
	int * ThisPathLength;
	uint64_t * BucketID;
	short Integrity;

	int total_reads;
	int total_writes;
	int issued_accesses;
	int finished_reads;
	int finished_writes;
	uint64_t cycles;
	int hit_latency;
	int ready_latency;

	uint64_t NumAccess;
	int MinReadyLatency;
	int MaxReadyLatency;
	uint64_t SumReadyLatency;
	int MinHitLatency;
	int MaxHitLatency;
	uint64_t SumHitLatency;

	bool debug_DRAM_access;
	DRAMSim::TransactionCompleteCB *read_cb;
	DRAMSim::TransactionCompleteCB *write_cb;

    float CPU_freq_factor;
};

class AddrShiftObj
{
public:
	// this class shifts addr[begin:begin+len-1] [off] bits ahead. (off < 0 means shifting towards lsb).
	AddrShiftObj(short head, short length, short offset)
	{
		begin = head;
		len = length;
		off = offset;
		if (offset == 0 || head < 0 || len < 0 || begin + len > 64 || head + offset < 0 || begin + len + offset > 64)
			begin = len = off = 0;
	}

	AddrShiftObj(short param[3])
	{
		begin = param[0];
		len = param[1];
		off = param[2];
		if (off == 0 || begin < 0 || len < 0 || begin + len > 64 || begin + off < 0 || begin + len + off > 64)
			begin = len = off = 0;
	}

	uint64_t shiftAddrBits(uint64_t x)
	{
		uint64_t move0, move1;
		move0 = x & (((1 << len) - 1) << begin);
		if (off >= 0)
			move1 = x & (((1 << off) - 1) << (begin + len));
		else
			move1 = x & (((1 << -off) - 1) << (begin + off));

//		printf("%ld, %ld\n", move0, move1);
		x -= move0 + move1;
		if (off >= 0)
			move0 <<= off, move1 >>= len;
		else
			move0 >>= -off, move1 <<= len;
//		printf("%ld, %ld\n", move0, move1);
		return x + move0 + move1;
	}
private:
	short begin;
	short len;
	short off;
};
#endif

#endif

