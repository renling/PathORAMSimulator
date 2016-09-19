#include "UnifiedBinPathORam.h"
#include "HierBinPathORam.h"

UnifiedBinPathORam::UnifiedBinPathORam()
{
	MaxHierarchy = 20;
	Hierarchy = 0;
	debug = 0;
	CompressedPosMap = 0;

	ValidBlockCount = new uint64_t [MaxHierarchy + 1];
#ifdef USE_DRAMSIM2
	Simulator = NULL;
#endif
	PLB = NULL;
	AccessCounts = NULL;
}

UnifiedBinPathORam::~UnifiedBinPathORam()
{
	if (PLB)
		delete PLB;
	if (AccessCounts)
		delete []AccessCounts;
	delete ORam;
	delete []ValidBlockCount;
	delete Addr;
}

int UnifiedBinPathORam::Configure(uint64_t workingSet, uint64_t oramSize, int blockSize, int blocksPerBucket, uint32_t maxPosMapSize, bool compressedPosMap, short integrity)
{
	assert(!Hierarchy);					// should not be called twice
	assert(workingSet > maxPosMapSize);	// otherwise, do not need external memory
	uint64_t BlockCount = oramSize / blockSize;
	uint64_t BucketCount = (BlockCount + blocksPerBucket - 1) / blocksPerBucket;
	int LevelCount = HighestBit(BucketCount);

	CompressedPosMap = compressedPosMap;
	Integrity = integrity;
	
	if (CompressedPosMap)
	{
		OuterCounterWidth = OUTER_COUNTER_WIDTH;
		InnerCounterWidth = INNER_COUNTER_WIDTH;
		PosMapScaleFactor = (blockSize * 8 - OuterCounterWidth) / InnerCounterWidth;
	}
	else if (!CompressedPosMap && Integrity)
		PosMapScaleFactor = blockSize * 8 / OUTER_COUNTER_WIDTH;	// has to be large enough
	else
		PosMapScaleFactor = blockSize * 8 / LevelCount;	// leaf labels in each block
		
	int logPosMapScaleFactor = log2(PosMapScaleFactor);
	PosMapScaleFactor = 1 << logPosMapScaleFactor;	// PosMapScaleFactor may have to a power of 2

	int i = 1;
	TotalWorkingSet = 0;
	ValidBlockCount[1] = workingSet / blockSize;
	while (ValidBlockCount[i] * blockSize > maxPosMapSize)	// iterate until fit into the specified on-chip storage
	{
		TotalWorkingSet += ValidBlockCount[i] * blockSize;
		ValidBlockCount[i+1] = (ValidBlockCount[i] + PosMapScaleFactor - 1) / PosMapScaleFactor;
		i++;
	}
	Hierarchy = i - 1;
	assert(Hierarchy < MaxHierarchy);
	FinalPositionMapSize = ValidBlockCount[Hierarchy+1] * blockSize;

	ValidBlockCount[0] = 0;
	for (int i = 0; i < Hierarchy; i++)
		ValidBlockCount[i+1] += ValidBlockCount[i];		// starting virtual address for each logical ORAM
	Addr = new int64_t [Hierarchy];

	// if (debug)
	{
		printf("Working set = %ld MB, ORAM tree is %ld MB\n", (int64_t) round(workingSet / 1024.0 / 1024), (int64_t) round(oramSize / 1024.0 / 1024));
                printf("Using ORAM Int: %d. Interval: %d\n", USE_ORAM_INTERVAL, O_INT);
		printf("Unified ORAM tree has %d levels. Each block contains %d leaf labels\n", LevelCount, PosMapScaleFactor);
		for (int i = 0; i < Hierarchy; i++)
			printf("ORAM %d has %ld blocks\n", i, ValidBlockCount[i+1] - ValidBlockCount[i]);		// ValidBlockCount is CSR format to enalbe this, notice!!
		printf("Total hierachy = %d, On-chip PosMap = %d Bytes\n", Hierarchy, FinalPositionMapSize);
	}

	if (CompressedPosMap)
	{
		AccessCounts = new int [ValidBlockCount[Hierarchy]];
		memset(AccessCounts, 0, sizeof(int) * ValidBlockCount[Hierarchy]);
	}

	ORam = new BinPathORam;
	ORam->Configure(TotalWorkingSet, oramSize, blockSize, blocksPerBucket);

	return 0;
}

int UnifiedBinPathORam::Initialize()
{
	NumAccess = NumDummy = 0;
	return ORam->Initialize();
}

int UnifiedBinPathORam::SetPLB(int capacity, int blockSize, int ways)
{
	assert(!PLB);
	PLB = new SimpleCache * [Hierarchy];
	PLB[0] = new SimpleCache(capacity, blockSize, ways);

	for (int i = 1; i < Hierarchy; i++)
		PLB[i] = UNIFIED_PLB ? PLB[0] : new SimpleCache(capacity, blockSize, ways);

	return 0;
}

int UnifiedBinPathORam::Access(int64_t addr, int hier_to_access, short RWoption)
{
	if (addr < 0)
		return ORam->Access(ORam->GetValidBlockCount(), -1, BinPathORam::dummy);

	assert(addr < (int64_t) ValidBlockCount[1]);
	GenerateAddr(addr);
	bool isWriteBack = BinPathORam::write_back & RWoption;
	assert(hier_to_access >= 0 && hier_to_access <= Hierarchy);
	if (isWriteBack)
		assert(hier_to_access == 0);

	int Traffic = 0;
	for (int i = hier_to_access - 1; i > 0; i--)
	{
		Traffic += ORam->Access(Addr[i], -1, BinPathORam::read | BinPathORam::erase);
		if (UpdatePLB(i))
			ORam->Access(PLB[i]->evicted->tag, -1, BinPathORam::write | BinPathORam::write_back);
	}
	Traffic += ORam->Access(Addr[0], -1, RWoption, NULL);

	return Traffic;
}

int UnifiedBinPathORam::UpdateCompressedPosMap(int hier_to_access)
{
	assert(CompressedPosMap);

	int Count = 0;
	for (int i = hier_to_access - 1; i >= 0; i--)
	{
		int64_t addr = Addr[i];
		assert(AccessCounts[addr] < (1 << InnerCounterWidth));
		AccessCounts[addr]++;

		if (AccessCounts[addr] < (1 << InnerCounterWidth))
			continue;			// inner counter not overflow, do nothing

		// otherwise, need to bump the outer counter and reset all inner counters
		int64_t Base = addr / PosMapScaleFactor * PosMapScaleFactor;
		for (int j = Base; j < Base + PosMapScaleFactor; j++)
		{
			AccessCounts[j] = 0;
			if (j == addr)
				continue;

			Count ++;
			if (!MODEL_STASH)	// TODO: too hacky ...
				continue;

			if (ORam->IsPresent(j))
				ORam->Access(j, -1, BinPathORam::read);
			else
				Count --;
			Count += ORam->BackgroundEvict_Count();
		}
	}
	return Count;
}

void UnifiedBinPathORam::GenerateAddr(int64_t addr)
{
	if (addr < 0 || addr >= (int64_t) ValidBlockCount[1])
		cout << addr << '\t' << ValidBlockCount[1] << endl;
	assert(addr >= 0 && addr < (int64_t) ValidBlockCount[1]);		// the region of virtual memory space
	Addr[0] = addr;
	for (int i = 1; i < Hierarchy; i++)
		Addr[i] = Addr[i-1] / PosMapScaleFactor;

	if (debug)
	{
		printf("Virtual block ID = { ");
		for (int i = 0; i < Hierarchy; i++)
			printf("%ld ", Addr[i]);
		printf("}\n");
	}

	for (int i = 0; i < Hierarchy; i++)
		Addr[i] += ValidBlockCount[i];
}

int UnifiedBinPathORam::ReadPLB()
{
	assert(PLB);
	int hier_to_access = 0;
	for (int i = 1; i < Hierarchy; i++)
	{
		if (PLB[i]->Read(Addr[i]) >= 0)
			break;
		else
			hier_to_access++;
	}

	return hier_to_access;
}

int UnifiedBinPathORam::UpdatePLB(int i)
{
	assert(PLB);
	assert(i > 0 && i < Hierarchy);
	return PLB[i]->Write(Addr[i]);
}

int UnifiedBinPathORam::BackgroundEvict()
{
	return ORam->BackgroundEvict(Hierarchy);
}

int UnifiedBinPathORam::BackgroundEvict_Count()
{
	return ORam->BackgroundEvict_Count(Hierarchy);
}


#ifdef USE_DRAMSIM2
void UnifiedBinPathORam::DRAMConfigure(short AES_pad_in_bits, float cpu_freq_factor, char * ini_dir, char * deviceConfFileName, char * systemConfFileName)
{
	struct ORAMInfoForDRAMSim info;
	info.hierarchy = 1;
	info.integrity = Integrity;	// PMMAC integrity only involves adding a flit to each bucket
	info.AES_pad_in_bits = AES_pad_in_bits + Integrity * HASH_LENGTH_IN_BITS * GetBlocksPerBucket();	// AES_pad is now really AES + SHA storage overhead
	info.cpu_freq_factor = CPU_freq_factor = cpu_freq_factor;
	info.blockSize[0] = GetBlockSize();
	info.blocksPerBucket[0] = GetBlocksPerBucket();
	info.blockCount[0] = GetBlockCount();
	info.bucketCount[0] = GetBucketCount();
	info.levelCount[0] = GetLevelCount();
	info.leafCount[0] = ORam->GetLeafCount();
	info.posMapScaleFactor[1] = PosMapScaleFactor;

	info.ini_dir = ini_dir;
	info.deviceConfFileName = deviceConfFileName; 
    info.systemConfFileName = systemConfFileName;

	Simulator = new ORAMForDRAMSim(info);
}

void UnifiedBinPathORam::SimulateLatency(int runlength)
{
	HitDelay = 0;
	ReadyDelay = 0;
	Simulator->Reset();
	Simulator->SimulateLatency(runlength, 0);
	HitDelay = Simulator->GetAveHitLatency();
	ReadyDelay = Simulator->GetAveReadyLatency();

	// convert into CPU cycles (adding AES latency too)
	int PerLevelLatency = AES_LATENCY;

	HitDelay = CPU_freq_factor * HitDelay + PerLevelLatency;
	ReadyDelay = CPU_freq_factor * ReadyDelay + PerLevelLatency;

	printf("Latency in CPU cycles: hit_delay = %d,  ready_delay = %d\n", HitDelay, ReadyDelay);
}
#endif

