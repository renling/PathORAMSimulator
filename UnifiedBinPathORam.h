#ifndef UNIFIED_BIN_PATH_ORAM
#define UNIFIED_BIN_PATH_ORAM

#include "BinPathORam.h"
using namespace std;

class UnifiedBinPathORam
{

public:
	UnifiedBinPathORam();
	~UnifiedBinPathORam();
	int Configure(uint64_t workingSet, uint64_t oramSize, int blockSize, int blocksPerBucket, uint32_t maxPosMapSize, bool compressedPosMap, short integrity);
	int Initialize();
	int SetPLB(int capacity, int blockSize, int ways);
	void SetMaxStashSize(int maxStashSize) {ORam->SetMaxStashSize(maxStashSize);}

	int Access(int64_t addr, int dummy, short RWoption);	// dummy is for matching HierBinPathORam
	void GenerateAddr(int64_t addr);
	int ReadPLB();
	int UpdatePLB(int i);
	int UpdateCompressedPosMap(int hier_to_access);

	int BackgroundEvict();
	int BackgroundEvict_Count();

	int GetHierarchy() {return Hierarchy;}
	int GetBlockSize() {return ORam->GetBlockSize();}
	int GetBlocksPerBucket() {return ORam->GetBlocksPerBucket();}
	int64_t GetBlockCount() {return ORam->GetBlockCount();}
	int64_t GetBucketCount() {return ORam->GetBucketCount();}
	int GetLevelCount() {return ORam->GetLevelCount();}
	int GetCurPathLength() {return ORam->GetCurPathLength();}
	int GetValidBlockCount(int i) {return ValidBlockCount[i+1] - ValidBlockCount[i];}
	int GetCurLocalCacheSize() {return ORam->GetCurLocalCacheSize();}
	int GetLastLocalCacheSize() {return ORam->GetLastLocalCacheSize();}
	int GetPeakLocalCacheSize() {return ORam->GetPeakLocalCacheSize();}
	int GetPosMapScaleFactor()	{return PosMapScaleFactor;}
	void ResetPeakLocalCacheSize() { ORam->ResetPeakLocalCacheSize();}
	bool IsPresent(int64_t id)	{return ORam->IsPresent(id);}

	uint64_t GetNumDummy() {return NumDummy;}
	uint64_t GetNumAccess() {return NumAccess;}
	uint32_t GetNumAccessInWindow() { return NumAccessInWindow;}
	uint32_t GetNumDummyInWindow() { return NumDummyInWindow;}
	void ResetNumDummyInWindow() { NumAccessInWindow = NumDummyInWindow = 0;}

	BinPathORam * ORam;
	int64_t *Addr;

#ifdef USE_DRAMSIM2
	friend class ORAMForDRAMSim;
	ORAMForDRAMSim * Simulator;
#endif
	void DRAMConfigure(short AES_pad_in_bits, float cpu_freq_factor, char * ini_dir, char * deviceConfFileName, char * systemConfFileName);
	void SimulateLatency(int runlength);
	int GetHitDelay() {return HitDelay;}
	int GetReadyDelay() {return ReadyDelay;}
private:
	int Hierarchy;
	int MaxHierarchy;
	int PosMapScaleFactor;
	int FinalPositionMapSize;
	uint64_t TotalWorkingSet;
	uint64_t * ValidBlockCount;

	SimpleCache ** PLB;
	bool CompressedPosMap;
	int OuterCounterWidth;
	int InnerCounterWidth;
	int * AccessCounts;

	short Integrity;

	uint32_t NumAccessInWindow;
	uint32_t NumDummyInWindow;
	uint64_t NumAccess;
	uint64_t NumDummy;

	int HitDelay;
	int ReadyDelay;
	float CPU_freq_factor;

	char debug;
};

#endif

