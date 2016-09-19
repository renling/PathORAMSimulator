
#include "HierBinPathORam.h"

HierBinPathORam::HierBinPathORam()
{
	MaxHierarchy = 20;
	Hierarchy = 0;
	debug = 0;

	ORam = new BinPathORam * [MaxHierarchy];
	WorkingSet = new uint64_t [MaxHierarchy];
	Util = new double [MaxHierarchy];
	BlockSize = new int [MaxHierarchy];
	BlocksPerBucket = new int [MaxHierarchy];
	BlockCount = new uint64_t [MaxHierarchy];
	BucketCount = new uint64_t [MaxHierarchy];
	LevelCount = new int [MaxHierarchy];
	LeafCount = new uint64_t [MaxHierarchy];

	PosMapScaleFactor = new int [MaxHierarchy];
#ifdef USE_DRAMSIM2
	Simulator = NULL;
#endif
	PLB = NULL;
}

HierBinPathORam::~HierBinPathORam() {}

int HierBinPathORam::Configure(uint64_t workingSet, double * utilization, int * HORamBlockSize, int * HORamBlocksPerBucket, uint32_t maxPosMapSize)
{
	assert(!Hierarchy);			// should not be called twice
	int i = 0;
	WorkingSet[0] = workingSet;
	assert(workingSet > maxPosMapSize);	// otherwise, do not need external memory

	uint64_t ValidBlockCount = 0;
	FinalPositionMapSize = WorkingSet[0];
	while (FinalPositionMapSize > (uint64_t) maxPosMapSize && i < MaxHierarchy)	// iterate until fit into the specified on-chip storage
	{
		ORam[i] = NULL;
		BlockSize[i] = HORamBlockSize[i];
		BlocksPerBucket[i] = HORamBlocksPerBucket[i];
		Util[i] = utilization[i];
		ValidBlockCount = WorkingSet[i] / BlockSize[i];
		BlockCount[i] = ceil(WorkingSet[i] / Util[i] / BlockSize[i]);
		BucketCount[i] = BlockCount[i] / BlocksPerBucket[i];
		BlockCount[i] = BucketCount[i] * BlocksPerBucket[i];
		LevelCount[i] = HighestBit(BucketCount[i]);
		LeafCount[i] = (BucketCount[i] + 1) / 2;
		
		PosMapScaleFactor[i+1] = (HORamBlockSize[i+1] * 8 / LevelCount[i]);
		int logPosMapScaleFactor = log2(PosMapScaleFactor[i+1]);
		PosMapScaleFactor[i+1] = 1 << logPosMapScaleFactor;

		// next level working set
		WorkingSet[i+1] = (ValidBlockCount + PosMapScaleFactor[i+1]- 1) / PosMapScaleFactor[i+1] * HORamBlockSize[i+1];
		FinalPositionMapSize = ValidBlockCount * LevelCount[i] / 8;
		i++;
	}
	Hierarchy = i;
	cout << BlockCount[Hierarchy-1] << endl;
	cout << LevelCount[Hierarchy-1] << endl;
	Addr = new int64_t [Hierarchy];

	// if (debug)
	{
		printf("Working set = %ld MB, ORAM tree is %ld MB\n", (int64_t) round(workingSet / 1024.0 / 1024), (int64_t) round(workingSet / Util[0] / 1024 / 1024));
		for (int i = 0; i < Hierarchy; i++)
			printf("ORAM %d has %d levels\n", i, LevelCount[i]);
		printf("Total hierachy = %d, On-chip PosMap = %ld Bytes\n", Hierarchy, FinalPositionMapSize);
	}

	for (int i = 0; i < Hierarchy; i++)
	{
		ORam[i] = new BinPathORam;
		ORam[i]->Configure(WorkingSet[i], WorkingSet[i] / Util[i], BlockSize[i], BlocksPerBucket[i]);
	}
	return 0;
}

int HierBinPathORam::GetBucketSizeInChunks(int i)
{
#ifdef USE_DRAMSIM2
	if (Simulator) return Simulator->BucketSizeInChunks[i];
	else
#endif
	return 0;
}

int HierBinPathORam::GetBucketSizeInBits(int i)
{
#ifdef USE_DRAMSIM2
	if (Simulator) return Simulator->BucketSizeInBits[i];
	else
#endif
	return 0;
}

int HierBinPathORam::Initialize()
{
	for (int i = 0; i < Hierarchy; i++)
	{
		ORam[i]->Initialize();
		assert(BlockCount[i] == ORam[i]->GetBlockCount());
		assert(BucketCount[i] == ORam[i]->GetBucketCount());
		assert(LevelCount[i] == ORam[i]->GetLevelCount());
	}
	NumAccess = NumDummy = 0;
	return 0;
}

int HierBinPathORam::SetPLB(int capacity, int blockSize, int ways)
{
	assert(!PLB);
	PLB = new SimpleCache * [Hierarchy];

	for (int i = 0; i < Hierarchy; i++)
		PLB[i] = new SimpleCache(capacity, blockSize, ways);

	return 0;
}

void HierBinPathORam::SetMaxStashSize(int maxStashSize)
{
	for (int i = 0; i < Hierarchy; i++)
		ORam[i]->SetMaxStashSize(maxStashSize);

}

int HierBinPathORam::Access(int64_t id, int hier_to_access, short RWoption)
{
	if(hier_to_access < 0 && hier_to_access > Hierarchy)
		cout << "hier to access: " << hier_to_access << endl;
	assert(hier_to_access >= 0 && hier_to_access <= Hierarchy);	// hierarchies to access does not exceed total

	if (debug == 1)
		cout << "block id: "<< id << "\t RWoption: " << (int) RWoption << endl;
	NumAccess ++;
	NumDummy += id < 0;

	int Traffic = 0;
	if (id < 0) 		// background eviction
	{
		for (int i = hier_to_access-1; i >= 0; i--)
			Traffic += ORam[i]->Access(ORam[i]->GetValidBlockCount(), -1, BinPathORam::dummy, NULL);
		return Traffic;
	}
	else for (int i = hier_to_access-1; i > 0; i--)
		assert (ORam[i] && !ORam[i]->LocalCacheFull());		// check stash not full

	assert(ORam[0]->GetCurLocalCacheSize() <= ORam[0]->GetMaxLocalCacheSize() + 1);

	// normal access
	GenerateAddr(id);
	for (int i = hier_to_access-1; i > 0; i--)
		Traffic += ORam[i]->Access(Addr[i], -1, BinPathORam::write, NULL);
	Traffic = ORam[0]->Access(Addr[0], -1, RWoption, NULL);

	return Traffic;
}

void HierBinPathORam::GenerateAddr(int64_t addr)
{
	Addr[0] = addr;
	for (int i = 1; i < Hierarchy; i++)
	{
		Addr[i] = Addr[i-1] / PosMapScaleFactor[i];
		if (debug)	cout << Addr[i] << endl;
	}
}

int HierBinPathORam::ReadPLB()
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

int HierBinPathORam::UpdatePLB(int i)
{
	assert(PLB);
	assert(i > 0 && i < Hierarchy);
	return PLB[i]->Write(Addr[i]);
}

bool HierBinPathORam::LocalCacheFull()
{
	for (int i = Hierarchy - 1; i >= 0; i--)
		if (ORam[i])
			if (ORam[i]->LocalCacheFull())
				return true;
	return false;
}

bool HierBinPathORam::LocalCacheAlmostFull()
{
	for (int i = Hierarchy - 1; i >= 0; i--)
		if (ORam[i]->LocalCacheAlmostFull())
			return true;
	return false;
}

void HierBinPathORam::BackgroundEvict()
{
	while (LocalCacheFull())
	{
		Access(-1, Hierarchy, BinPathORam::dummy);
	}
}

int HierBinPathORam::BackgroundEvict_Count()
{
	int count = 0;
	while (LocalCacheFull())
	{
		Access(-1, Hierarchy, BinPathORam::dummy);
		count ++;
	}
	return count;
}

#ifdef USE_DRAMSIM2
void HierBinPathORam::DRAMConfigure(short AES_pad_in_bits, float cpu_freq_factor, char * ini_dir, char * deviceConfFileName, char * systemConfFileName)
{
	struct ORAMInfoForDRAMSim info;
	assert(Hierarchy);		// must call configure first
	info.hierarchy = Hierarchy;
	info.AES_pad_in_bits = AES_pad_in_bits;
	info.integrity = Integrity = false;			// Recursive ORAM does not have integrity right now
	info.cpu_freq_factor = CPU_freq_factor = cpu_freq_factor;
	for (int i = 0; i < 8; i++)
	{
		info.blockSize[i] = BlockSize[i];
		info.blocksPerBucket[i] = BlocksPerBucket[i];
		info.blockCount[i] = BlockCount[i];
		info.bucketCount[i] = BucketCount[i];
		info.levelCount[i] = LevelCount[i];
		info.leafCount[i] = LeafCount[i];
		info.posMapScaleFactor[i] = PosMapScaleFactor[i];
	}
	info.ini_dir = ini_dir;
	info.deviceConfFileName = deviceConfFileName; 
    info.systemConfFileName = systemConfFileName;
	
	Simulator = new ORAMForDRAMSim(info);
}

void HierBinPathORam::SimulateLatency(int runlength)
{
	HitDelay = new int [Hierarchy];
	ReadyDelay = new int [Hierarchy];
	for (int hier = 0; hier < Hierarchy; )
	{
		Simulator->Reset();
		Simulator->SimulateLatency(runlength, hier);
		HitDelay[hier] = Simulator->GetAveHitLatency();
		ReadyDelay[hier] = Simulator->GetAveReadyLatency();

		// convert into CPU cycles (adding AES latency too)
		int direction = 2 - (bool) ((Simulator->Strategy & 4) && hier);
		assert(direction == 1 || direction == 2);
		int PerLevelLatency = AES_LATENCY;
		PerLevelLatency *= direction;

		HitDelay[hier] = CPU_freq_factor * HitDelay[hier] + PerLevelLatency * hier;
		ReadyDelay[hier] = CPU_freq_factor * ReadyDelay[hier] + PerLevelLatency * hier;

		printf("Latency in CPU cycles: hier = %d, hit_delay = %d,  ready_delay = %d\n", hier + 1, HitDelay[hier], ReadyDelay[hier]);
		
		if (hier == Hierarchy-1)
			break;
		else
			hier = Hierarchy-1;
	}
	cout << "PosMap portion = " << 1.0 * (ReadyDelay[0]) / ReadyDelay[Hierarchy-1]  << endl;

}
#endif

int findParam(const char * filename, const char * keyword)
{
	char tmp[4096];
	ifstream file(filename);
	assert(file.good());
	file.getline(tmp, 4096);
	while(strncmp(tmp, keyword, strlen(keyword)))
		file.getline(tmp, 4096);
	file.close();
	cout << keyword <<atoi(tmp + strlen(keyword)) << endl;
	return atoi(tmp + strlen(keyword));
}

#ifdef USE_DRAMSIM2
ORAMForDRAMSim::ORAMForDRAMSim(ORAMInfoForDRAMSim info)
{
	read_cb = new DRAMSim::Callback<ORAMForDRAMSim, void, unsigned, uint64_t, uint64_t>(this, &ORAMForDRAMSim::read_complete);
	write_cb = new DRAMSim::Callback<ORAMForDRAMSim, void, unsigned, uint64_t, uint64_t>(this, &ORAMForDRAMSim::write_complete);
	string deviceConf(info.ini_dir), systemConf(info.ini_dir);
	deviceConf.append(info.deviceConfFileName);
	systemConf.append(info.systemConfFileName);		
	//deviceConf.append("DDR3_micron_16M_8B_x8_sg15_ORAM.ini");
	//systemConf.append("system_ORAM.ini");		

	NumChannels = findParam(systemConf.c_str(), "NUM_CHANS=");
	NumColumns = findParam(deviceConf.c_str(), "NUM_COLS=") / findParam(deviceConf.c_str(), "NUM_BANKS=");
	ChunkSize = findParam(deviceConf.c_str(), "BL=") * findParam(systemConf.c_str(), "JEDEC_DATA_BUS_BITS=") / 8;
	RowWidthInBytes = NumChannels * NumColumns * ChunkSize;
	char ch = HighestBit(NumChannels) - 1;
	AddrMapping = new AddrShiftObj(6 + ch, 7, 7 - ch);

	Info = info;
	Hierarchy = Info.hierarchy;
	Integrity = Info.integrity;
	CPU_freq_factor = Info.cpu_freq_factor;
	BlockSize = Info.blockSize;
	BlocksPerBucket = Info.blocksPerBucket;
	BlockCount = Info.blockCount;
	BucketCount = Info.bucketCount;
	LevelCount = Info.levelCount;
	LeafCount = Info.leafCount;
	PosMapScaleFactor = Info.posMapScaleFactor;
	
	BaseAddr = new uint64_t [Hierarchy + 1];
	BucketSizeInBits = new int [Hierarchy];
	BucketSizeInChunks = new int [Hierarchy];
	KPacking = new short [Hierarchy];
	LonelySubTrees = new int64_t [Hierarchy];
	LastSubTreeSize = new int [Hierarchy];
	Path = new uint64_t  [Hierarchy];
	ThisPathLength = new int [Hierarchy];

	BucketID = new uint64_t [LevelCount[0]];

	
	BaseAddr[0] = 0;
	for (int i = 0; i < Hierarchy; i++)
	{
		BucketSizeInBits[i] = BlocksPerBucket[i] * (BlockSize[i] * 8 + 2 * LevelCount[i] + 1) +  Info.AES_pad_in_bits;
											// AES_pad_in_bits = 64 is the counter trick

		BucketSizeInChunks[i] = ceil(BucketSizeInBits[i] / 8.0 / ChunkSize);
		printf("L = %d, bucket size = %d bits\n", LevelCount[i], BucketSizeInBits[i]);

		int BucketsInRow = RowWidthInBytes / ChunkSize / BucketSizeInChunks[i];
		if (BucketsInRow > 0)
			KPacking[i] = floor(log2(BucketsInRow+1));
		else
			KPacking[i] = 1;
		if (KPacking[i] >= LevelCount[i] / 2.0)
			KPacking[i] = ceil(LevelCount[i] / 2.0);
		int NumComSubTree = LevelCount[i] / KPacking[i];					// # levels of complete subtrees
		LonelySubTrees[i] = 1 << (NumComSubTree * KPacking[i]);
		int alpha = 1 << (KPacking[i]);
		int64_t alphaSum = (pow(alpha, NumComSubTree) - 1) / (alpha - 1);	// 1 + alpha + alpha^2 + .... + alpha^{L-1}
		BaseAddr[i+1] = BaseAddr[i] + (alphaSum + 0) * RowWidthInBytes;		// complete subtrees

		LastSubTreeSize[i] = (1 << (LevelCount[i] % KPacking[i])) - 1;		// # of incomplete subtrees
		LastSubTreeSize[i] *= BucketSizeInChunks[i] * ChunkSize;
		BaseAddr[i+1] += LastSubTreeSize[i] * (LonelySubTrees[i] + 1);		// incomplete subtrees

/*
		cout <<  alpha << endl;
		cout <<  alphaSum << endl;
		cout <<  BaseAddr[i+1] << endl;
		cout <<  LastSubTreeSize[i] << endl;
		cout <<  LonelySubTrees[i] << endl;
		cout <<  BaseAddr[i+1] << endl;
*/
	}

	uint64_t DRAM_Capacity = 1 << ((int) ceil(log2(BaseAddr[Hierarchy])) - 20);
	mem = DRAMSim::getMemorySystemInstance(deviceConf.c_str(), systemConf.c_str(), "./", "Ascend", DRAM_Capacity);
	mem->RegisterCallbacks(read_cb, write_cb, NULL);

	debug_DRAM_access = 0;
	Strategy = 7;

	Reset();
}

ORAMForDRAMSim::~ORAMForDRAMSim()
{
	delete read_cb;
	delete write_cb;
	delete AddrMapping;
	delete BaseAddr;
	delete []BucketSizeInChunks;
	delete []KPacking;
	delete []Path;
	delete []ThisPathLength;
	delete []BucketID;
}

void ORAMForDRAMSim::Reset()
{
	NumAccess = 0;
	MaxReadyLatency = SumReadyLatency = MaxHitLatency = SumHitLatency = 0;
	MinReadyLatency = MinHitLatency = 1024 * 1024 * 1024;
	for (int i = 0; i < 100; i++)
		mem->update();
}

void ORAMForDRAMSim::SimulateLatency(int runlength, int hier)
{
//	srand (time(NULL));
	srand (0);
	for (int iter = 0; iter < runlength; iter++)
		SimulateOneAccess(hier);
		printf("# read = %d, # write = %d. Latency in DRAM cycles: hit_latency = %f [%d, %d],  ready_latency = %f [%d, %d]\n", total_reads, total_writes, GetAveHitLatency(), MinHitLatency, MaxHitLatency, GetAveReadyLatency(), MinReadyLatency, MaxReadyLatency);

}

int ORAMForDRAMSim::SimulateOneAccess(int hier_to_access)
{
	NumAccess++;
	cycles = 0;
	total_reads = total_writes = issued_accesses = finished_reads = finished_writes = 0;
	for (int i = 0; i <= hier_to_access; i++)
		Path[i] = rand() % LeafCount[i] + BucketCount[i] - LeafCount[i];

	for (int i = hier_to_access; i >= 0; i--)
	{
		ThisPathLength[i] = HighestBit(Path[i] + 1);
		assert(LevelCount[i] - ThisPathLength[i] <= 1);
		total_reads += ThisPathLength[i] * BucketSizeInChunks[i];
		total_writes += ThisPathLength[i] * BucketSizeInChunks[i];
	}


	if (Strategy & 4)
		for (char rw = 0; rw < 2; rw ++)
			for (int i = hier_to_access; i >= 0; i--)
			{
				bool blocking = !rw;
				if (i == 0 && rw == 1)
				{
					while ((int64_t) (cycles - hit_latency) * CPU_freq_factor <= AES_LATENCY)
						tick();
				}
				SimulateOneORAM(i, (bool) rw, blocking);
			}
	else
		for (int i = hier_to_access; i >= 0; i--)
			for (char rw = 0; rw < 2; rw ++)
			{
				assert(false);
				bool blocking = rw;
				SimulateOneORAM(i, (bool) rw, blocking);
			}
	drain();
	ready_latency = cycles;

	SumHitLatency += hit_latency;
	SumReadyLatency += ready_latency;
	assert(hit_latency > 0);
	MinHitLatency = min(hit_latency, MinHitLatency);
	MaxHitLatency = max(hit_latency, MaxHitLatency);
	MinReadyLatency = min(ready_latency, MinReadyLatency);
	MaxReadyLatency = max(ready_latency, MaxReadyLatency);

	return ready_latency;
}

void ORAMForDRAMSim::SimulateOneORAM(int h, bool rw, bool blocking)
{
	if (blocking)
		drain();
	int BucketSizeInCBytes = BucketSizeInChunks[h] * ChunkSize;		// buckets should be multiple of chunks
	int CurPathLength = ThisPathLength[h];
	int KPack = KPacking[h];

	BucketID[CurPathLength-1] = Path[h];
	for (int l = CurPathLength - 1; l > 0; l--)
		BucketID[l-1] = (BucketID[l] - 1) / 2;

	uint64_t BaseAddr1 = BaseAddr[h], BaseAddr2 = 0, ThisSubTreeSize = RowWidthInBytes;
	int64_t offset = 0, subBucketID = 0;
	for (int l = 0; l < CurPathLength; l++)
	{
		if (l % KPack == 0)
		{
			if (l > 0)		BaseAddr1 += ( 1 << (l - KPack) ) * (uint64_t) RowWidthInBytes;					// next level of subtrees

			bool IsIncomplete = l == (CurPathLength/KPack*KPack) && LastSubTreeSize[h];				// is last level subtree && incomplete tree?
			ThisSubTreeSize = IsIncomplete ? LastSubTreeSize[h] : RowWidthInBytes;
			int64_t tmp = 1;
			tmp <<= l;
			offset = (BucketID[l] + 1 - tmp);
			assert(offset < 2 * tmp);
			BaseAddr2 = BaseAddr1 + offset * ThisSubTreeSize;										// locate the subtree
			subBucketID = 0;
		}
		else
			subBucketID = 2 * subBucketID + BucketID[l] - 2 * BucketID[l-1];						// locate the bucket in the subtree

		for (int j = 0; j < BucketSizeInChunks[h]; j++)												// read bucket
		{
			uint64_t addr = BaseAddr2 + subBucketID * BucketSizeInCBytes + j * ChunkSize;
			if (addr >= BaseAddr[h+1])
				printf("%d, %d, %d, %ld, %ld, %ld, %ld, %ld, %ld, %ld\n", h, l, KPack, BaseAddr[h], BaseAddr1, BaseAddr2, offset, addr, BaseAddr[h+1], ThisSubTreeSize);
			assert(addr < BaseAddr[h+1]);
			addRequest(addr, rw);
		}
	}
}

void ORAMForDRAMSim::addRequest(uint64_t addr, bool rw)
{
	addr = AddrMapping->shiftAddrBits(addr);
	if (debug_DRAM_access)
		printf("\t\taccessing DRAM address %ld\n", addr);
	mem->addTransaction(rw, addr);
	issued_accesses ++;
	while (issued_accesses - finished_reads - finished_writes > 16 * NumChannels)
		tick();
}

void ORAMForDRAMSim::read_complete(unsigned id, uint64_t address, uint64_t clock_cycle)
{
	finished_reads ++;
	if (finished_reads == total_reads)
	{
		hit_latency = cycles;		// record hit delay;
		assert(cycles > 0);
		if (debug_DRAM_access)
			printf("[DRAM Response] entire read complete: %d reads at cycle=%lu\n", issued_accesses, clock_cycle);
	}
	if (debug_DRAM_access)
		printf("[DRAM Response] read complete: address 0x%16lx at cycle=%lu\n", address, clock_cycle);
}

void ORAMForDRAMSim::write_complete(unsigned id, uint64_t address, uint64_t clock_cycle)
{
	finished_writes ++;
	if (debug_DRAM_access)
		printf("[DRAM Response] write complete: address 0x%16lx at cycle=%lu\n", address, clock_cycle);
	if (debug_DRAM_access)
	if (finished_writes == total_writes)
		printf("[DRAM Response] entire write complete: %d writes at cycle=%lu\n", issued_accesses / 2, clock_cycle);
}
#endif

