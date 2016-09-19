#ifndef BIN_PATH_ORAM
#define BIN_PATH_ORAM

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <list>
#include <vector>
#include <stdint.h>
#include <map>
#include <memory.h>

#define NDEBUG
#undef NDEBUG	// checking asserts
#include <assert.h>

#include "Parameters_ORAM.h"
#include "SimpleCache.h"

#define USE_DRAMSIM2
#include "DRAMSim.h"

using namespace std;

class LocalCacheLine
{
public:
	LocalCacheLine() {}
	LocalCacheLine(int64_t id, int64_t * PositionMap) {ID = id; Position = PositionMap + id;}
		
	int64_t ID;		// program address
	int64_t * Position;	// leaf mapped to
};

unsigned int HighestBit(uint64_t x);
uint64_t ORamVALookup(uint64_t block_tag, uint64_t page_size_set = 0);

#ifdef USE_DRAMSIM2
class AddrShiftObj;
class ORAMForDRAMSim;
#endif

class BinPathORam 
{
public:

	enum AccessMode {
		read 	= 1,
		write 	= 2,
		erase	= 4,
		write_back = 8,
		fetch 	= 16,
		dummy 	= 32
	};

	BinPathORam();
	~BinPathORam();

	int Configure(uint64_t workingSet, uint64_t oramSize, int blockSize, int blocksPerBucket);
	int Configure(uint64_t workingSet, int levels, int blocksPerBucket);	
	int Initialize();
	int Access(int64_t addr, int64_t position, short RWoption, char* data = NULL);
	int AccessOneBlock(int64_t id, int64_t position, short RWoption);
	int BackgroundEvict(int margin = 0);
	int BackgroundEvict_Count(int margin = 0);
	int InsecureBackgroundEvict();
	int ForegroundEvict1(int64_t position);
	int ForegroundEvict2(int64_t position);
	
	int GetBlockSize() {return BlockSize;}
	int GetBlocksPerBucket() {return BlocksPerBucket;}
	uint64_t GetBlockCount() {return BlockCount;}
	uint64_t GetBucketCount() {return BucketCount;}
	uint64_t GetLeafCount() {return LeafCount;}
	int GetLevelCount() {return LevelCount;}
	int GetCurPathLength() {return CurPathLength;}
	int GetValidBlockCount() {return ValidBlockCount;}
	int GetCurLocalCacheSize() {return LocalCache.size();}
	int GetLastLocalCacheSize() {return LastLocalCacheSize;}
	int GetPeakLocalCacheSize() {return PeakLocalCacheSize;}
	void ResetPeakLocalCacheSize() {PeakLocalCacheSize = 0;}
	bool IsPresent(int64_t id)	{return Present[id];}
	
	uint64_t * GetLCSZHist() {return LCSZHist;}
	uint64_t * GetBucketStatHist() {return BucketStatHist;}
	uint64_t GenerateRandLeaf();
	
	int GetMaxLocalCacheSize() {return MaxLocalCacheSize;}
	void SetMaxStashSize(int maxStashSize);
	bool LocalCacheFull(int margin = 0) {return GetCurLocalCacheSize() >= MaxLocalCacheSize - margin - BlocksPerBucket * LevelCount;}
	bool LocalCacheAlmostFull() {return GetCurLocalCacheSize() >= 0.7 * (MaxLocalCacheSize - BlocksPerBucket * LevelCount);}
	uint64_t GetNumAccess() {return NumAccess;}
	uint64_t GetNumDummy() {return NumDummy;}
	int GetuselessDummy() {return uselessDummy;}

	void EnableHistogram(int histSize, char recordPeak);
	void RecordHistogram();
	void DumpHistogram(const char * filename);

	void PrintLocalCache();

	std::vector<int64_t> ServerSees;

public:
	uint64_t WorkingSet;	
	uint64_t ORAMSize;
	int BlockSize;
	int BlocksPerBucket;

	uint64_t ValidBlockCount;
	uint64_t BlockCount;
	uint64_t BucketCount;
	uint64_t LeafCount;
	int LevelCount;

	bool * Present;
	int64_t * PositionMap;			// block - leaf mapping
	int64_t * ProgAddr;			// program address of a block, used for locating interesting block
	char* ProgData;				// program data, for functional simulation
	
	std::list<LocalCacheLine> LocalCache; 	// Local cache
	std::list<LocalCacheLine>::iterator iter;
	int PeakLocalCacheSize;			// peak local cache size ever seen
	int LastLocalCacheSize;			// local cache size after in last access after path read
	int MaxLocalCacheSize;			// max local cache size allowed	

	short RWOption;
	int CurPathLength;

	int64_t * EvictQueue;			// Evict Queue
	short* EvictQueueCount;
	short* EvictQueueCountUB;
	int64_t * CurPathBuffer;		// store the path on this access
	int EvictDepth;

	int ReadPath(int64_t interest, uint64_t leaf);
	int ReadPathForOneBlock(int64_t interest, uint64_t leaf);
	void Remap(int64_t interest, uint64_t newLeaf);
	int ScanCurPath(int64_t interest, uint64_t leaf);
	int ScanStash(int64_t interest, uint64_t leaf);
	int WritePath(int64_t interest, uint64_t leaf);
	
	int ReadBucket(uint64_t bucketIndex);
	int WriteBucket(uint64_t bucketIndex, uint64_t leafToEvict);
	bool LegalReside(uint64_t position, uint64_t bucketIndex);
	
	inline void ResetEvictQueue();
	//inline int PosDiff(int64_t blockPosition, int64_t oldPosition);
	int PosDiff(int64_t blockPosition, int64_t oldPosition);
	inline int FindSpaceOnPath(LocalCacheLine Block, int64_t oldPosition);
	inline void EvictBlock(LocalCacheLine Block, int height);

	uint64_t * LCSZHist;
	int HistSize;
	char RecordOption;
	uint64_t * BucketStatHist;

	uint64_t NumAccess;
	uint64_t NumDummy;

	char debug;
	int * randPool;
	int randPoolIdx;
	int uselessDummy;
};

#endif

