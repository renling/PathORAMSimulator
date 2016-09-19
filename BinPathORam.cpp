
#include "BinPathORam.h"

using namespace std;

unsigned int HighestBit(uint64_t x)
{
	unsigned int hb = 0;
	while (x)
	{
		hb++;
		x>>=1;
	}
	return hb;
}

uint64_t ORamVALookup(uint64_t vaddr, uint64_t page_size_set)
{
	static map<uint64_t, uint64_t> ORamVA;
	static uint64_t page_size = 0;
	// vaddr in the same page (e.g. 4K) are mapped to consecutive oaddrs
	// for UORAM, set it large such that an initial access misses all the PLBs.

	if (page_size == 0)
	{
		page_size = page_size_set;
        printf("page_size: %ld\n", page_size_set);
        return 0;
	}

	uint64_t page_idx = vaddr / page_size;
	uint64_t page_offset = vaddr % page_size;
	uint64_t oaddr;
	map<uint64_t, uint64_t>::iterator ORamVAIt = ORamVA.find(page_idx);
	if (ORamVAIt == ORamVA.end())
	{
		ORamVA.insert(pair<uint64_t, uint64_t>(page_idx, ORamVA.size()));
		oaddr = ORamVA.size() - 1;
	}
	else
		oaddr = ORamVAIt->second;
	return oaddr * page_size + page_offset;
}

BinPathORam::BinPathORam() {}

BinPathORam::~BinPathORam() {}

int BinPathORam::Configure(uint64_t workingSet, uint64_t oramSize, int blockSize, int blocksPerBucket)
{
	debug = 0;
	WorkingSet = workingSet;
	ORAMSize = oramSize;
	BlockSize = blockSize;
	BlocksPerBucket = blocksPerBucket;

	ValidBlockCount = WorkingSet / blockSize;
	BlockCount = ORAMSize / BlockSize;
	BucketCount = BlockCount / BlocksPerBucket;		// we may lost several blocks here, but who cares
	BlockCount = BucketCount * BlocksPerBucket;
	LevelCount = HighestBit(BucketCount);
	LeafCount = (BucketCount + 1) / 2;

	assert (ValidBlockCount <= BlockCount);
	return 0;
}

int BinPathORam::Configure(uint64_t workingSet, int levels, int blocksPerBucket)
{
	WorkingSet = workingSet;
	LevelCount = levels;
	BlocksPerBucket = blocksPerBucket;

	BlockSize = 1;
	ValidBlockCount = WorkingSet / BlockSize;
	BucketCount = (1 << levels) - 1;
	LeafCount = (BucketCount + 1) / 2;
	BlockCount = BucketCount * BlocksPerBucket;

	return 0;
}

int BinPathORam::Initialize()
{
//	cout<<"LevelCount: "<<LevelCount<<"\tBucketCount: "<<BucketCount<<"\tLeafCount: "<<LeafCount<<endl;

	Present = new bool [ValidBlockCount + 1];	// one extra address reserved for dummy
	PositionMap = new int64_t [ValidBlockCount + 1];
	ProgAddr = new int64_t [BlockCount];
	EvictQueue = new int64_t [LevelCount * BlocksPerBucket];
	EvictQueueCount = new short [LevelCount];
	CurPathBuffer = new int64_t [LevelCount * BlocksPerBucket];
	assert(Present && PositionMap && ProgAddr && EvictQueue && EvictQueueCount && CurPathBuffer);

	memset(Present, 0, sizeof(bool) * (ValidBlockCount + 1));	// all blocks initiated as "not present"
	for (uint64_t i = 0; i < ValidBlockCount + 1; i += 1)
	{
		int64_t newPosition = GenerateRandLeaf();			// all position maps initialized to random
		PositionMap[i] = newPosition;					// blocks in the same superblock have the same leaf
	}
	memset(ProgAddr, -1, sizeof(int64_t) * BlockCount);  // all blocks in the tree initiated as dummy

	MaxLocalCacheSize = 1024 * 1024;					// max local cache size allowed, initialized to almost inf
	LastLocalCacheSize = PeakLocalCacheSize = 0;
	NumAccess = NumDummy = 0;

	HistSize = 0;
	RecordOption = 0;
	LCSZHist = NULL;
	BucketStatHist = NULL;

	return 0;
}


int BinPathORam::Access(int64_t id, int64_t position, short RWoption, char* data)
{
	NumAccess ++;
	if (id < 0)
		id = ValidBlockCount;
	assert(id < (int64_t) ValidBlockCount + 1);
	assert(!LocalCacheFull() || (RWoption & dummy));
	RWOption = RWoption;

	if (RWOption & write_back)			// This is a block evicted from LLC, append it to stash without accessing any path
	{
		assert(!Present[id]);			// This block should not exist at this point
		LocalCache.push_back(LocalCacheLine(id, PositionMap));
		Present[id] = true;
		return 0;
	}

	ResetEvictQueue();
	int Traffic = 0;

	if ((RWoption & dummy) && position >= 0)
		PositionMap[id] = position;
	// Step 1: Lookup Position map
	int64_t oldPosition	= PositionMap[id];
	int64_t newPosition = GenerateRandLeaf();
	assert(position < 0 || position == oldPosition);	// Leaf, if provided, should match position map
#ifdef BE_SECURE
	assert(0);
	ServerSees.push_back(oldPosition);					// the path accessed is revealed
#endif
	if (debug == 1)
		cout<<"Access Block "<<id<<", mapped to Leaf "<<oldPosition<<", remapped to Leaf "<<newPosition
		<<", this path length: "<<HighestBit(oldPosition + 1) << endl;

	// Step 2: Remap
	Remap(id, newPosition);

	// Step 3: Read path
	Traffic += ReadPath(id, oldPosition);

	if (!Present[id])					// New block, not currently in ORAM
	{
		if (RWOption & read)			// Should not read a non-existent block unless RWoption == dummy
		{
			cout << "Read non existent block!!!! " << id << endl;
			assert(RWOption & dummy);
		}
		else if (RWOption & write)		// create the block if it's a write
		{
			Present[id] = true;
			LocalCache.push_back(LocalCacheLine(id, PositionMap));
		}
	}

	LastLocalCacheSize = LocalCache.size();		// record local cache size
	PeakLocalCacheSize = PeakLocalCacheSize > LastLocalCacheSize? PeakLocalCacheSize: LastLocalCacheSize;
	if (debug == 1)
		PrintLocalCache();

	// Step 4: Scan Stash
	ScanCurPath(id, oldPosition);
	ScanStash(id, oldPosition);

	if (debug == 1)
		PrintLocalCache(),	getchar();

	// Step 5: Write path back
	Traffic += WritePath(id, oldPosition);

	return Traffic;
}

int BinPathORam::AccessOneBlock(int64_t id, int64_t position, short RWoption)
{
	NumAccess ++;
	assert(id >= 0 && id < (int64_t) ValidBlockCount + 1);
	assert(!LocalCacheFull() || (RWoption & dummy));
	RWOption = RWoption;

	if (RWOption & write_back)			// This is a block evicted from LLC, append it to stash without accessing any path
	{
		assert(!Present[id]);			// This block should not exist at this point
		LocalCache.push_back(LocalCacheLine(id, PositionMap));
		PositionMap[id] = PositionMap[id];	// mapped to the 'leader' of this superblock
		Present[id] = true;
		return 0;
	}

	// ResetEvictQueue();
	int Traffic = 0;

	// Step 1: Lookup Position map
	if ((RWoption & dummy) && position >= 0)
		PositionMap[id] = position;
	int64_t oldPosition	= PositionMap[id];
	int64_t newPosition = GenerateRandLeaf();
	assert(position < 0 || position == oldPosition);	// Leaf, if provided, should match position map

	if (debug == 1)
		cout<<"Access Block "<<id<<", mapped to Leaf "<<oldPosition<<", remapped to Leaf "<<newPosition
		<<", this path length: "<<HighestBit(oldPosition + 1) << endl;

	// Step 2: Remap
	Remap(id, newPosition);

	// Step 3: Read path
	// Traffic += ReadPath(id, oldPosition);
	Traffic += ReadPathForOneBlock(id, oldPosition);

	if (!Present[id])					// New block, not currently in ORAM
	{
		if (RWOption & read)			// Should not read a non-existent block unless RWoption == dummy
		{
			cout << "Read non existent block!!!! " << id << endl;
			// assert(RWOption & dummy);
		}
		else if (RWOption & write)		// create the block if it's a write
		{
			Present[id] = true;
			LocalCache.push_back(LocalCacheLine(id, PositionMap));
		}
	}

	LastLocalCacheSize = LocalCache.size();		// record local cache size
	PeakLocalCacheSize = max(PeakLocalCacheSize, LastLocalCacheSize);
	if (debug == 1)
		PrintLocalCache();

	// Step 4: Scan Stash
	if (debug == 1)
		PrintLocalCache(),	getchar();

	// Step 5: Write path back
	return Traffic;
}

int BinPathORam::ReadPath(int64_t interest, uint64_t leaf)
{
	uint64_t bucketIndex = leaf;
	CurPathLength = HighestBit(leaf + 1);
	for (int height = 0; height < CurPathLength; height++)
	{
		assert(bucketIndex || (bucketIndex == 0 && height == CurPathLength - 1));
		assert(bucketIndex * BlocksPerBucket < BlockCount);
		memcpy(CurPathBuffer + height * BlocksPerBucket, ProgAddr + bucketIndex * BlocksPerBucket, sizeof(int64_t) * BlocksPerBucket);
		bucketIndex = (bucketIndex - 1) / 2;
	}
	return CurPathLength * BlocksPerBucket;
}

int BinPathORam::ReadPathForOneBlock(int64_t interest, uint64_t leaf)
{
	uint64_t bucketIndex = leaf;
	CurPathLength = HighestBit(leaf + 1);
	bool BlockFound = false;
	for (int height = 0; height < CurPathLength; height++)
	{
		assert(bucketIndex || (bucketIndex == 0 && height == CurPathLength - 1));
		assert(bucketIndex * BlocksPerBucket < BlockCount);
		for (int j = 0; j < BlocksPerBucket; j++)
			if (ProgAddr[bucketIndex * BlocksPerBucket + j] == interest)
			{
				LocalCacheLine Block(ProgAddr[bucketIndex * BlocksPerBucket + j], PositionMap);
				LocalCache.push_back(Block);
				BlockFound = true;
				ProgAddr[bucketIndex * BlocksPerBucket + j] = -1;
			}
		bucketIndex = (bucketIndex - 1) / 2;
	}
	assert(BlockFound == !(RWOption & dummy));
	return CurPathLength * BlocksPerBucket;
}

void BinPathORam::Remap(int64_t interest, uint64_t newLeaf)
{
	PositionMap[interest] = newLeaf;
}

int BinPathORam::ScanCurPath(int64_t interest, uint64_t leaf)
{
	for (int j = 0; j < BlocksPerBucket * CurPathLength; j++)
	//for (int j = BlocksPerBucket * CurPathLength-1; j >= 0; j--)	// a discovery by Xiao. RAW ORAM is sensitive to this
	{
		LocalCacheLine Block(CurPathBuffer[j], PositionMap);
		if (Block.ID == -1)
			continue;
		else if (Block.ID == interest)		// block of interest goes into stash
		{
			CurPathBuffer[j] = -1;
			LocalCache.push_back(Block);
			continue;											// cannot call FindSpaceOnPathOnPath twice for this blok
		}
		int newHeight = FindSpaceOnPath(Block, leaf);
		if (newHeight > -1)
		{
			EvictBlock(Block, newHeight);
			CurPathBuffer[j] = -1;
		}
	}

	return 0;
}

int BinPathORam::ScanStash(int64_t interest, uint64_t leaf)
{
	bool BlockFound = false;
	iter = LocalCache.begin();
	while (iter != LocalCache.end())
	{
		if ((*iter).ID == interest)
		{
			BlockFound = true;
			if (RWOption & erase)
			{
				Present[(*iter).ID] = false;
				iter = LocalCache.erase(iter);
				continue;
			}
		}

		int newHeight = FindSpaceOnPath(*iter, leaf);
		if (newHeight != -1)
		{
			EvictBlock(*iter, newHeight);
			iter = LocalCache.erase(iter);
		}
		else
			iter++;
	}
	assert(BlockFound == !(RWOption & dummy));
	return 0;
}

int BinPathORam::WritePath(int64_t interest, uint64_t leaf)
{
	int Traffic = 0;
	uint64_t bucketIndex = leaf;
	CurPathLength = HighestBit(leaf + 1);
	for (int height = 0; height < CurPathLength; height++)
	{
		assert(bucketIndex || (bucketIndex == 0 && height == CurPathLength - 1));
		assert(bucketIndex < BucketCount);
		Traffic ++;
		int64_t * src = EvictQueue + height * BlocksPerBucket;
		memcpy(ProgAddr + bucketIndex * BlocksPerBucket, src, sizeof(int64_t) * BlocksPerBucket); // write bucket
		bucketIndex = (bucketIndex - 1) / 2;
	}
	return Traffic * BlocksPerBucket;
}


int BinPathORam::BackgroundEvict(int margin)
{
	int Traffic = 0;
	while (LocalCacheFull(margin))
	{
		NumDummy++;
		Traffic += Access(ValidBlockCount, -1, dummy, NULL);
	}
	return Traffic;
}

int BinPathORam::BackgroundEvict_Count(int margin)
{
	int count = 0;
	while (LocalCacheFull(margin))
	{
		NumDummy++;
		count++;
		Access(ValidBlockCount, -1, dummy, NULL);
	}
	return count;
}

int BinPathORam::InsecureBackgroundEvict()
{
	int Traffic = 0;
	int64_t insecure_position = -1;
	while (LocalCacheFull())
	{
		NumDummy++;
//		Traffic += Access(LocalCache.front().ID, -1, dummy, NULL);
		if (insecure_position != *(LocalCache.front().Position))
		 	insecure_position = *(LocalCache.front().Position);
	 	else
	 		insecure_position = GenerateRandLeaf();
 //		*(LocalCache.front().Position) = GenerateRandLeaf();
		//*(LocalCache.front().Position) = GenerateRandLeaf();
		Traffic += Access(ValidBlockCount, insecure_position, dummy, NULL);
//		*(LocalCache.front().Position) = insecure_position;
	}
	return Traffic;
}

int BinPathORam::ReadBucket(uint64_t bucketIndex)
{
	assert(bucketIndex * BlocksPerBucket < BlockCount);
	for (int j = 0; j < BlocksPerBucket; j++)
	{
		LocalCacheLine Block(ProgAddr[bucketIndex * BlocksPerBucket + j], PositionMap);
		if (Block.ID >= 0)
			LocalCache.push_back(Block);
		ProgAddr[bucketIndex * BlocksPerBucket + j] = -1;				
	}
	return BlocksPerBucket;
}

bool SortToLeaf(LocalCacheLine x, LocalCacheLine y, uint64_t leaf)
{
	// first align to the same level and then xor
	
	uint64_t l[2];
	l[0] = *(x.Position) + 1,	l[1] = *(y.Position) + 1;	
	for (int j = 0; j < 2; j++)
		switch (HighestBit(l[j]) - HighestBit(leaf+1))
		{
			//case 1: l[j] >>= 1; break;
			case 0: break;
			//case -1: l[j] <<= 1; break;
			default: cout << l[j] << '\t'<< leaf << endl; assert(0);
		}
	return (l[0] ^ leaf) < (l[1] ^ leaf);
}

class SorterToLeaf {
      uint64_t leaf_to_sort_;
public:
      SorterToLeaf(uint64_t leaf) {	leaf_to_sort_ = leaf;}
      bool operator()(LocalCacheLine const x, LocalCacheLine const y) const {
            return SortToLeaf(x, y, leaf_to_sort_ );
      }
};


int BinPathORam::WriteBucket(uint64_t bucketIndex, uint64_t leafToEvict)
{
	assert(bucketIndex * BlocksPerBucket < BlockCount);
	
	// leafToEvict = 0 skips the sorting phase and evicts arbitrary blocks
	if (leafToEvict > 0)
	{
		// sort the stash based on the leaf we are evicting to
		LocalCache.sort(SorterToLeaf(leafToEvict));			
	}		
	
	// evict up Z blocks to that bucket
	iter = LocalCache.begin();
	int legal_writes = 0;
	while (iter != LocalCache.end() && legal_writes < BlocksPerBucket)
	{
		if (LegalReside((*(*iter).Position), bucketIndex))
		{
	//		printf("writing block %ld to bucket %ld[%d]\n", (*iter).ID, bucketIndex, legal_writes);	
			assert(ProgAddr[bucketIndex * BlocksPerBucket + legal_writes] == -1);
			ProgAddr[bucketIndex * BlocksPerBucket + legal_writes] = (*iter).ID;
			iter = LocalCache.erase(iter);
			legal_writes ++;
		}
		else
			iter ++;
	}
	
	// fill the rest with dummy
	for (int j = legal_writes; j < BlocksPerBucket; j++)
		ProgAddr[bucketIndex * BlocksPerBucket + j] = -1;
				
	return BlocksPerBucket;
}

bool BinPathORam::LegalReside(uint64_t position, uint64_t bucketIndex)
{
	position ++;
	bucketIndex ++;
	
	position >>= (HighestBit(position) - HighestBit(bucketIndex));	
	return (position == bucketIndex);
}

int BinPathORam::ForegroundEvict2(int64_t leaf)		// read/write path + siblings
{
	CurPathLength = HighestBit(leaf + 1);
	
	// read the path and the siblings
	uint64_t bucketIndex = leaf;
	for (int height = 0; height < CurPathLength-1; height++)
	{
		assert(bucketIndex);
		assert(bucketIndex * BlocksPerBucket < BlockCount);
		
		ReadBucket(bucketIndex);
		ReadBucket(bucketIndex - 1 + 2 * (bucketIndex % 2));

		bucketIndex = (bucketIndex - 1) / 2;
	}
	ReadBucket(0);

	// write the path and the siblings
	bucketIndex = leaf;	
	for (int height = 0; height < CurPathLength-1; height++)
	{
		assert(bucketIndex);
		assert(bucketIndex * BlocksPerBucket < BlockCount);
		
		WriteBucket(bucketIndex, leaf);
		WriteBucket(bucketIndex - 1 + 2 * (bucketIndex % 2), leaf);

		bucketIndex = (bucketIndex - 1) / 2;
	}
	WriteBucket(0, leaf);
	
	return CurPathLength * BlocksPerBucket;	
}

int BinPathORam::ForegroundEvict1(int64_t leaf)		// Shi et al. swap buckets + Gentry or random
{
	CurPathLength = HighestBit(leaf + 1);
	
	uint64_t bucketIndex = 0;
	for (int height = 0; height < CurPathLength-1; height++)
	{
	//	cout << bucketIndex << endl;
		assert(bucketIndex * BlocksPerBucket < BlockCount);
		
		// read bucket and two children
		ReadBucket(bucketIndex);
		ReadBucket(2 * bucketIndex + 1);
		ReadBucket(2 * bucketIndex + 2);
		
		//	move to children
		bucketIndex = 2 * bucketIndex + 1 + ((leaf + 1) >> (CurPathLength-2-height)) % 2;
					
		// write bucket and parent
		WriteBucket(bucketIndex, leaf);			
		WriteBucket((bucketIndex-1)/2, leaf);	
		
		// write sibling
		WriteBucket(bucketIndex - 1 + 2 * (bucketIndex % 2), 0);				
	}
	
	return CurPathLength * BlocksPerBucket;	
}

inline void BinPathORam::ResetEvictQueue()
{
	memset(EvictQueue, -1, sizeof(int64_t) * LevelCount * BlocksPerBucket);
	memset(EvictQueueCount, 0, sizeof(short) * LevelCount);
}


uint64_t BinPathORam::GenerateRandLeaf() { return (rand() % LeafCount) + BucketCount - LeafCount;}

int BinPathORam::PosDiff(int64_t blockPosition, int64_t oldPosition)	// the highest intersection (towards leaf) of two paths
{

	assert(blockPosition > 0 && oldPosition > 0);
	blockPosition++, oldPosition++;
	switch (HighestBit(blockPosition) - CurPathLength)
	{
		case 1: blockPosition >>= 1; break;
		case 0: break;
		case -1: blockPosition <<= 1; break;
		default: assert(0);
	}
	return HighestBit(blockPosition ^ oldPosition);
}

inline int BinPathORam::FindSpaceOnPath(LocalCacheLine Block, int64_t oldPosition)
{
	assert (Block.ID >= 0);
	int64_t blockPosition = *(Block.Position);

	for (int height = PosDiff(blockPosition, oldPosition); height < CurPathLength; height++)	// try to add this block to EvictQueue
		if (EvictQueueCount[height] < BlocksPerBucket)
			return height;
	return -1;
}

inline void BinPathORam::EvictBlock(LocalCacheLine Block, int height)
{
	EvictQueue[height * BlocksPerBucket + EvictQueueCount[height]] = Block.ID;
	EvictQueueCount[height]++;
	if (debug == 1)
		cout<<"Evicting Block " << Block.ID << " to height "<< height <<endl;
}

void BinPathORam::SetMaxStashSize(int maxLocalCacheSize)
{
	assert(maxLocalCacheSize > BlocksPerBucket * LevelCount);
	MaxLocalCacheSize = maxLocalCacheSize;
}

void BinPathORam::PrintLocalCache()
{
	cout<<"Printing Local Cache: ";
	for (iter = LocalCache.begin(); iter != LocalCache.end(); iter++)
		cout<<"("<<iter->ID<<", "<<*(iter->Position)<<"), ";
	cout<<'\n';
}

void BinPathORam::EnableHistogram(int histSize, char recordOption)
// recordOption: 0 disable; 1 after path read; 2 after path write; 3 bucket status
{
	HistSize = histSize;
	LCSZHist = new uint64_t [HistSize];
	if (RecordOption < 3)
	{
		memset(LCSZHist, 0, sizeof(uint64_t) * HistSize);
		RecordOption = recordOption;
	}
	if (RecordOption == 3 || RecordOption == 4)
	{
		BucketStatHist = new uint64_t [LevelCount * (BlocksPerBucket+1)];
		memset(BucketStatHist, 0, sizeof(uint64_t) * LevelCount * (BlocksPerBucket+1));
	}
}

void BinPathORam::RecordHistogram()
{
	int HistIndex;
	if (RecordOption < 3)
	{
		HistIndex = (RecordOption == 1)? GetLastLocalCacheSize() : GetCurLocalCacheSize();
		HistIndex = HistIndex < HistSize ? HistIndex : HistSize - 1;
		LCSZHist[HistIndex]++;
	}
	if (RecordOption == 3)
	{
		for (int height = 0; height < LevelCount; height++)
			BucketStatHist[height * (BlocksPerBucket+1) + EvictQueueCount[height]]++;
	}
	else if (RecordOption == 4)
	{
		int ValidBlocksOnPath = 0;
		for (int height = 0; height < LevelCount; height++)
			ValidBlocksOnPath += EvictQueueCount[height];
		BucketStatHist[ValidBlocksOnPath]++;
	}
}

void BinPathORam::DumpHistogram(const char * filename)
{
	ofstream fout;
	fout.open(filename, ios_base::out);
	if (RecordOption < 3)
		for (int j = 0; j < HistSize; j++)
			fout << LCSZHist[j] << '\n';
	if (RecordOption == 3)
		for (int height = 0; height < LevelCount; height++)
		{
			for (int j = 0; j < (BlocksPerBucket+1); j++)
				fout<<BucketStatHist[height * (BlocksPerBucket+1) + j]<<'\t';
			fout << '\n';
		}
	else if (RecordOption == 4)
		for (int num = 0; num < LevelCount * BlocksPerBucket + 1; num++)
			fout << BucketStatHist[num] << '\n';
	fout.close();
}

