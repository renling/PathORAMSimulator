#include "BinPathORam.h"
#include "HierBinPathORam.h"
#include "UnifiedBinPathORam.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <memory.h>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>
using namespace std;
//#define NDEBUG
#include <assert.h>

int workingSetInMB = 8;
int oramSizeInMB = 16;
int levels = 0;
int dataBlockSize = 64;
int posBlockSize = 32;
int dataZ = 3;
int posZ = 3;
int superBlockSize = 1;
int maxLCSZ = ORAM_LCSZ;
int onchipPosMapSize = 512 * 1024;
int runlength = 4;
bool worstCase = 0;
int histSize = 1000;
short integrity = 0;

double RAWA = 1.0;
int ExtraDummy = 4;		// this is X
double RAWEarlyReshuffleRate = 1.0;
int LogBurstLength = 0;

ofstream fout;
string outFileName("./exp_res/result_");

void SingleORAMTest();
void RAWORAMTest();
void RAWPPORAMTest();
void RAWBurstORAMTest(double);
void HierORAMTest();
void UnifiedORAMTest();

void BESecureTest();
int LocalCacheHist(const char * filename);

void SimDRAMLatency();

void Miscellany();

int main(int argc, char* argv[])
{	
	bool illegal_arg = false;
	for (int n = 2; n < argc; n++) 
	{
		assert(argv[n][0] == '-');
		switch (argv[n][1]) 
		{
		case 's' :
			workingSetInMB = atoi(&argv[n][2]);
			break;
		case 'S' :
			oramSizeInMB = atoi(&argv[n][2]);
			break;
		case 'L' :
			levels = atoi(&argv[n][2]);
			break;
		case 'B' : 
			dataBlockSize = atoi(&argv[n][2]);
			break; 
		case 'b' : 
			posBlockSize = atoi(&argv[n][2]);
			break; 
		case 'Z' : 
			dataZ = atoi(&argv[n][2]);
			break;   
		case 'z' : 
			posZ = atoi(&argv[n][2]);
			break;    
		case 'p' : 
			superBlockSize = atoi(&argv[n][2]);
			break;         
		case 'm' :
			maxLCSZ = atoi(&argv[n][2]);
			break;
		case 'w' :
			assert(argv[n][2] == '0' || argv[n][2] == '1');
			worstCase = argv[n][2] - '0';
			break;
		case 'r' :
			runlength = atoi(&argv[n][2]);
			break;
		case 'h' :
			histSize = atoi(&argv[n][2]);
			break;
		case 'i':
			assert(atoi(&argv[n][2]) < 3);
			integrity = atoi(&argv[n][2]);
			break;
		case 'A' :
			RAWA = atof(&argv[n][2]);
			break;
		case 'X' :
			ExtraDummy = atoi(&argv[n][2]);
			break;		
		case 'x' :
			RAWEarlyReshuffleRate = atof(&argv[n][2]);
			break;		
		case 't' :
			LogBurstLength = atoi(&argv[n][2]);
			break;		
		default :
			illegal_arg = true;			
		}
	}

	if (argc < 2)
		illegal_arg = true;
	
	if (illegal_arg)
	{
		cout<< "Usage: ./BinPathORam_test test_type arg_list\n"
			<< "\t-s : working set in MB\n"
			<< "\t-S : ORAM tree size in MB, incompatible with -L\n"
			<< "\t-L : levels, incompatible with -S\n"
			<< "\t-B : data ORAM block size in Bytes\n"
			<< "\t-b : PosMap ORAM block size in Bytes\n"
			<< "\t-Z : data ORAM blocks per bucket\n"	
			<< "\t-z : PosMap ORAM blocks per bucket\n"
			<< "\t-m : max Stash size\n" 
			<< "\t-w : simulate worst case?\n"
			<< "\t-r : simulation length\n"
			<< "\t-i : integrity verification\n"	
		exit(0);
	}

	outFileName.append(argv[1]);
	for (int i = 2; i < argc; i++) 
		outFileName.append(argv[i]);

	ifstream try_read(outFileName.c_str());
	if (try_read.good())
	{
//	  	cout<<outFileName.c_str()<<" exists"<<endl;
		try_read.close();
//		exit(0);
	}
	fout.open(outFileName.c_str(), ios_base::out);	
	fout.close();	

	if (! strcmp(argv[1], "single"))
		SingleORAMTest();
	else if (! strcmp(argv[1], "raw"))
		RAWORAMTest();
	else if (! strcmp(argv[1], "rawpp"))
		RAWPPORAMTest();	
	else if (! strcmp(argv[1], "rawburst"))
	{
		RAWBurstORAMTest(0.053687091);
		RAWBurstORAMTest(0.107374182);
		RAWBurstORAMTest(0.161061274);
		RAWBurstORAMTest(0.214748365);
		RAWBurstORAMTest(0.268435456);
		RAWBurstORAMTest(0.536870912);
	}
	else if (! strcmp(argv[1], "hier"))
		HierORAMTest();
	else if (! strcmp(argv[1], "BEsecure"))
		BESecureTest();
	else if (! strcmp(argv[1], "localcache"))
	{
		int sum = 0;
		for (int i = 0; i < 1; i++)
			sum += LocalCacheHist(outFileName.c_str());
		cout << sum / 1 << endl;
	}
	else if (! strcmp(argv[1], "dramsim"))
		SimDRAMLatency();
	else if (! strcmp(argv[1], "miscellany"))
		Miscellany();

	return 0;
} 

int64_t divceil(int64_t x, int64_t y) { return x < 0? 0 : (x + y - 1) / y; }
int64_t roundto(int64_t x, int64_t y) { return divceil(x, y) * y; }

int ModelEarlyReshuffle(double prob, int level)
{
	int res = 0;
	const int rand_max = RAND_MAX / 2;
	for (int j = 0; j < level; j++)		
		if (rand() % rand_max < rand_max * prob)
			res++;
	return res;		
}

void RAWPPORAMTest()
{
	const int PinBW = 16; // Bytes per cycle
	
	bool unified_oram = false;
	bool use_plb = unified_oram;
	
	BinPathORam PathORam, RAWORam;
	
	int PathORamHalfLat = 0, ROR_Lat = 0, ROR_XOR_Lat = 0, ROW_Lat = 0, RW_Lat = 0, RW_Bkt_Lat = 0, RW_Pos_Bkt_Lat = 0;
	int RAWLevelCount = 0, RAWPosLevelCount = 0;
	int hierarchy = 0;
	
	int ZplusA = dataZ + RAWA;
	int ZplusAplusX = ZplusA + ExtraDummy;
	
	if (unified_oram)
	{
		// Path ORAM latency
		int PathORamZ = 4;
		PathORam.Configure(workingSetInMB * (int64_t) 1024 * 1024, workingSetInMB * 2 * (int64_t) 1024 * 1024, dataBlockSize, PathORamZ);
		int PathORamLevelCount = PathORam.GetLevelCount();
		int PathH_Bkt = 64 + PathORamZ * 2 * roundto(PathORam.GetLevelCount(), 32);
		PathORamHalfLat = divceil(PathH_Bkt, 8); + PathORamZ * dataBlockSize;	// encSeed + Z(U+L+B)
		PathORamHalfLat = divceil(PathORamHalfLat, PinBW) * PathORamLevelCount;
					
		// RAW ORAM latency
		RAWORam.Configure(workingSetInMB * (int64_t) 1024 * 1024, oramSizeInMB * (int64_t) 1024 * 1024, dataBlockSize, dataZ);
		RAWLevelCount = RAWORam.GetLevelCount();
	
			// RO overhead per level in Bytes
		int RORH_Bkt = 64 + ZplusAplusX + dataZ * roundto(RAWLevelCount, 32) + dataZ * (HighestBit(ZplusAplusX) + 1);	// one encseed, valid, Z * U + Z * log(Z+A+X)
		int RORD_Bkt = dataBlockSize;
		int ROWH_Bkt = 1; //  update valid bits in the clear
		int RWH_Bkt = RORH_Bkt + dataZ * roundto(RAWLevelCount, 32);	// extra: Z*L
	
		RWH_Bkt  = divceil(RWH_Bkt, 8);
		ROWH_Bkt = divceil(ROWH_Bkt, 8);
		RORH_Bkt = divceil(RORH_Bkt, 8);
		
			// RW overhead per level in Bytes	
		int RWR_Bkt = RWH_Bkt + dataZ * dataBlockSize;
		int RWW_Bkt = RWH_Bkt + ZplusAplusX * dataBlockSize;	// header, Z leaf labels, Z+A+X blocks
		
		// latency per level
		int ROR_Bkt_Lat = divceil(RORH_Bkt, PinBW) + divceil(RORD_Bkt, PinBW);
		int ROW_Bkt_Lat = divceil(ROWH_Bkt, PinBW);
		
		ROR_Lat = ROR_Bkt_Lat * RAWLevelCount;
		ROW_Lat = ROW_Bkt_Lat * RAWLevelCount;
		ROR_XOR_Lat = divceil(RORH_Bkt, PinBW) * RAWLevelCount + divceil(RORD_Bkt, PinBW);	// latency with XOR
				
		RW_Bkt_Lat = divceil(RWR_Bkt, PinBW) + divceil(RWW_Bkt, PinBW);	
		RW_Lat = RW_Bkt_Lat * RAWLevelCount;
	}
	
	else
	{
		// recursive Path ORAM
		int PathORamZ = 4;
		uint64_t workingSet = workingSetInMB * (int64_t) 1024 * 1024;
		double ORAMUtil[10] = {workingSetInMB * 1.0 / oramSizeInMB, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
		int HORamBlockSize[10] = {dataBlockSize};
		int HORamBlocksPerBucket[10] = {PathORamZ};	 
		for (int i = 1; i < 10; i++)
			HORamBlockSize[i] = posBlockSize, HORamBlocksPerBucket[i] = PathORamZ;
	
		HierBinPathORam PathORam;
		PathORam.Configure(workingSet, ORAMUtil, HORamBlockSize, HORamBlocksPerBucket, onchipPosMapSize);
	
		hierarchy = PathORam.GetHierarchy();
		for (int h = 0; h < hierarchy; h++)
		{
			int PathH_Bkt = 64 + PathORamZ * 2 * roundto(PathORam.GetLevelCount(h), 32);
			PathORamHalfLat += (divceil(PathH_Bkt, 8) + PathORamZ * PathORam.GetBlockSize(h)) * PathORam.GetLevelCount(h);
		}	
		PathORamHalfLat = divceil(PathORamHalfLat, PinBW);	
		PathORamHalfLat = 5028 / 2;	
	
		// recursive RAW ORAM
		HORamBlocksPerBucket[0] = dataZ;
		for (int i = 1; i < 10; i++)
			HORamBlocksPerBucket[i] = dataZ;
			
		HierBinPathORam RAWORam;
		RAWORam.Configure(workingSet, ORAMUtil, HORamBlockSize, HORamBlocksPerBucket, onchipPosMapSize);
		hierarchy = RAWORam.GetHierarchy();
		
		RAWLevelCount = RAWORam.GetLevelCount(0);
		for (int h = 1; h < hierarchy; h++)
			RAWPosLevelCount += RAWORam.GetLevelCount(h);
			
		for (int h = 0; h < hierarchy; h++)
		{
			int this_level_count = RAWORam.GetLevelCount(h);
			int this_block_size = RAWORam.GetBlockSize(h);
			
				// RO overhead per level in Bytes
			int RORH_Bkt = 64 + ZplusAplusX + dataZ * roundto(this_level_count, 32) + dataZ * (HighestBit(ZplusAplusX) + 1);	// one encseed, valid, Z * U + Z * log(Z+A+X)	
			int ROWH_Bkt = 1;
			int RWH_Bkt = RORH_Bkt + dataZ * roundto(this_level_count, 32);	// extra: Z*L
					
			RWH_Bkt  = divceil(RWH_Bkt, 8);
			ROWH_Bkt = divceil(ROWH_Bkt, 8);
			RORH_Bkt = divceil(RORH_Bkt, 8);
			
			int RWR_Bkt = RWH_Bkt + dataZ * this_block_size;
			int RWW_Bkt = RWH_Bkt + ZplusAplusX * this_block_size;	// header, Z leaf labels, Z+A+X blocks
						
			// latency per level
			int ROR_Bkt_Lat = divceil(RORH_Bkt, PinBW) + divceil(this_block_size, PinBW);
			int ROW_Bkt_Lat = divceil(ROWH_Bkt, PinBW);
					
			// RW overhead per level in Bytes			
			ROR_Lat += ROR_Bkt_Lat * this_level_count;
			ROW_Lat += ROW_Bkt_Lat * this_level_count;
			ROR_XOR_Lat = divceil(RORH_Bkt, PinBW) * this_level_count + divceil(this_block_size, PinBW);	// latency with XOR
					
			int this_bkt_lat =  divceil(RWR_Bkt, PinBW) + divceil(RWW_Bkt, PinBW);	
			
			RW_Lat += this_bkt_lat * this_level_count;			
			if (h==0)
				RW_Bkt_Lat = this_bkt_lat;
			else
				RW_Pos_Bkt_Lat = this_bkt_lat;
		}		
	}
		
	int DRAM_Lat = 50;
	cout << "Path ORAM half latency = " << PathORamHalfLat << endl;	
	cout << "RAW ORAM RO_R latency = " << ROR_Lat << endl;
	cout << "RAW ORAM RO_W latency = " << ROW_Lat << endl;
	cout << "RAW ORAM RWBk latency = " << RW_Bkt_Lat << endl;
	cout << "RAW ORAM RW.. latency = " << RW_Lat << endl;	

	uint64_t TraceSegSize = 100000;
	int * read_addr = new int [TraceSegSize];
	int * idle_time = new int [TraceSegSize];
	ifstream traceFile;

	UnifiedBinPathORam * UORam = new UnifiedBinPathORam;
	if (unified_oram)
	{	
		UORam->Configure(workingSetInMB * (int64_t) 1024 * 1024, oramSizeInMB * (int64_t) 1024 * 1024, dataBlockSize, dataZ, onchipPosMapSize, 1, 0);
	
		hierarchy = UORam->GetHierarchy();
		if (use_plb)
		{
			UORam->SetPLB(8192, dataBlockSize, 1);
			uint64_t page_size = dataBlockSize;
			for (int i = 1; i < hierarchy; i++)
				page_size *= UORam->GetPosMapScaleFactor();
			ORamVALookup(0, page_size);
		}
	}
	
	ofstream fout;
	fout.open("trace_sim.out");
	fout << "bench runtime:DRAM Path RAW RAW_wd RAW_XOR resptime:Path RAW RAW_wd RAW_XOR\n";
	cout << "bench runtime:DRAM Path RAW RAW_wd RAW_XOR resptime:Path RAW RAW_wd RAW_XOR\n";
	fout.close();
	
	string benches[12] = {"tpcc", "ycsb", "astar", "bzip2", "gcc", "gobmk", "h264ref", "libquantum", "mcf", "omnetpp", "perlbench", "sjeng"};
	for (int bench_id = 0; bench_id < 12; bench_id++)
	{
		// read trace
		string traceFileName = string("input/");
		traceFileName.append(benches[bench_id]).append(".trace");
		traceFile.open(traceFileName.c_str());
					
		// statistics to gather	
		uint64_t TotalTime_RAW = 0, TotalTime_RAW_wd = 0, TotalTime_RAW_XOR = 0, TotalTime_RAW_wd_XOR = 0, TotalTime_Path = 0, TotalTime_DRAM = 0;
		uint64_t NewTotalTime_RAW = 0, NewTotalTime_RAW_wd = 0, NewTotalTime_RAW_XOR = 0, NewTotalTime_RAW_wd_XOR = 0, NewTotalTime_Path = 0, NewTotalTime_DRAM = 0;
		uint64_t TotalRespTime_RAW = 0, TotalRespTime_RAW_wd = 0, TotalRespTime_RAW_XOR = 0, TotalRespTime_RAW_wd_XOR = 0, TotalRespTime_Path = 0, TotalRespTime_DRAM = 0;
		uint64_t RWAccessCount_RAW = 0, RWAccessCount_RAW_wd = 0, RWAccessCount_RAW_XOR = 0, RWAccessCount_RAW_wd_XOR = 0;
		uint64_t ROAccessCount = 0, ORAMAccessCount = 0;
		int RW_Lat_Done_wd = 0, RW_Lat_Done_wd_XOR = 0;
	    
	    int thisTraceSegSize = TraceSegSize;
	    while (thisTraceSegSize == TraceSegSize)
	    {
	    	// read in some traces
	    	thisTraceSegSize = 0;
	    	assert(traceFile.good());
	    	for (uint64_t i = 0; i < TraceSegSize && !traceFile.eof(); i++, thisTraceSegSize++)
	    	{
	    		// if (i == 0) cout << "reading more traces from " << traceFile.tellg() << endl;
	    		uint64_t addr_in, idle_time_in;
				traceFile >> addr_in >> idle_time_in;
				if (bench_id < 2)
					addr_in = ORamVALookup(addr_in, 4096) / 64;
				assert(addr_in < 2147483646 && idle_time_in < 10000000);
				read_addr[i] = addr_in;
				idle_time[i] = idle_time_in;	
			}
			if (traceFile.eof( ))
				thisTraceSegSize--;
			
			for (uint64_t i = 0; i < thisTraceSegSize; i++, ORAMAccessCount++)
			{				
    			int	hier_to_access = 1;
				if (unified_oram)
				{	
					assert(use_plb);
					// access plb to figure out how many uoram accesses should be made
					UORam->GenerateAddr(read_addr[i]);
					hier_to_access = UORam->ReadPLB() + 1;
					for (int i = hier_to_access - 1; i > 0; i--)
           				UORam->UpdatePLB(i);
				}

    			
				for (int h = 0; h < hier_to_access; h++)
				{
					int this_idle_time;
					if (h == 0)
						this_idle_time = idle_time[i];
					else
						this_idle_time = 0;	
					
					// RAW ORAM without de-amortization	
					if (ROAccessCount - RWAccessCount_RAW * RAWA >= 0)
					{
						NewTotalTime_RAW += max(this_idle_time, ROW_Lat + RW_Lat);
						RWAccessCount_RAW++;
					}
					else
						NewTotalTime_RAW += max(this_idle_time, ROW_Lat);				
					NewTotalTime_RAW += DRAM_Lat + ROR_Lat;	// online latency
					
					// RAW ORAM with de-amortization
					if (ROAccessCount - RWAccessCount_RAW_wd * RAWA >= 0)
					{
						NewTotalTime_RAW_wd += max(this_idle_time, ROW_Lat + RW_Lat - RW_Lat_Done_wd); 
						RWAccessCount_RAW_wd++;
						RW_Lat_Done_wd = 0;
					}
					else 
					{
						NewTotalTime_RAW_wd += max(this_idle_time, ROW_Lat);
						RW_Lat_Done_wd += max(0, this_idle_time - ROW_Lat);
					}						
					NewTotalTime_RAW_wd += DRAM_Lat + ROR_Lat;	// online latency
					
					// RAW ORAM without de-amortization but with XOR	
					if (ROAccessCount - RWAccessCount_RAW_XOR * RAWA >= 0)
					{
						NewTotalTime_RAW_XOR += max(this_idle_time, ROW_Lat + RW_Lat);
						RWAccessCount_RAW_XOR++;
					}
					else
						NewTotalTime_RAW_XOR += max(this_idle_time, ROW_Lat);
					NewTotalTime_RAW_XOR += DRAM_Lat + ROR_XOR_Lat;	// online latency
										
					// RAW ORAM with de-amortization and XOR
					if (ROAccessCount - RWAccessCount_RAW_wd_XOR * RAWA >= 0)
					{
						NewTotalTime_RAW_wd_XOR += max(this_idle_time, ROW_Lat + RW_Lat - RW_Lat_Done_wd_XOR); 
						RWAccessCount_RAW_wd_XOR++;
						RW_Lat_Done_wd_XOR = 0;
					}
					else 
					{
						NewTotalTime_RAW_wd_XOR += max(this_idle_time, ROW_Lat);
						RW_Lat_Done_wd_XOR += max(0, this_idle_time - ROW_Lat);
					}						
					NewTotalTime_RAW_wd_XOR += DRAM_Lat + ROR_XOR_Lat;	// online latency
									
					// Path ORAM
					NewTotalTime_Path += max(this_idle_time, PathORamHalfLat);
					NewTotalTime_Path += DRAM_Lat + PathORamHalfLat;	
					
					// early reshuffle for RAW
					if (unified_oram)
					{
						NewTotalTime_RAW += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);	
						NewTotalTime_RAW_wd += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);
						NewTotalTime_RAW_XOR += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);
						NewTotalTime_RAW_wd_XOR += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);	
					}
					else
					{
						NewTotalTime_RAW += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);	
						NewTotalTime_RAW_wd += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);
						NewTotalTime_RAW_XOR += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);
						NewTotalTime_RAW_wd_XOR += (RW_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);
					
						NewTotalTime_RAW += (RW_Pos_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWPosLevelCount);	
						NewTotalTime_RAW_wd += (RW_Pos_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWPosLevelCount);	
						NewTotalTime_RAW_XOR += (RW_Pos_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWPosLevelCount);
						NewTotalTime_RAW_wd_XOR += (RW_Pos_Bkt_Lat) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWPosLevelCount);		
					}								
				}
											
				// DRAM
				NewTotalTime_DRAM += idle_time[i];
				NewTotalTime_DRAM += DRAM_Lat;	// a flat DRAM latency	
				
				// record response time
				TotalRespTime_RAW += NewTotalTime_RAW - (TotalTime_RAW + idle_time[i]);
				TotalRespTime_RAW_wd += NewTotalTime_RAW_wd - (TotalTime_RAW_wd + idle_time[i]);
				TotalRespTime_RAW_XOR += NewTotalTime_RAW_XOR - (TotalTime_RAW_XOR + idle_time[i]);
				TotalRespTime_RAW_wd_XOR += NewTotalTime_RAW_wd_XOR - (TotalTime_RAW_wd_XOR + idle_time[i]);
				TotalRespTime_Path += NewTotalTime_Path - (TotalTime_Path + idle_time[i]);
											
				// update time
				TotalTime_RAW = NewTotalTime_RAW;
				TotalTime_RAW_wd = NewTotalTime_RAW_wd;
				TotalTime_RAW_XOR = NewTotalTime_RAW_XOR;
				TotalTime_RAW_wd_XOR = NewTotalTime_RAW_wd_XOR;
				TotalTime_Path = NewTotalTime_Path;
				TotalTime_DRAM = NewTotalTime_DRAM;
				
				ROAccessCount += hier_to_access;
			}
		}
		traceFile.close();
	
		TotalRespTime_Path /= ORAMAccessCount;
		TotalRespTime_RAW /= ORAMAccessCount;
		TotalRespTime_RAW_wd /= ORAMAccessCount;
		TotalRespTime_RAW_XOR /= ORAMAccessCount;
		TotalRespTime_RAW_wd_XOR /= ORAMAccessCount;
	
		fout.open("trace_sim.out", ios::out | ios::app);
		fout << benches[bench_id] << ' ' << TotalTime_DRAM << ' ' << TotalTime_Path << ' ' << TotalTime_RAW << ' ' << TotalTime_RAW_wd << ' ' << TotalTime_RAW_XOR << ' ' << TotalTime_RAW_wd_XOR << ' ' << TotalRespTime_Path << ' ' << TotalRespTime_RAW << ' ' << TotalRespTime_RAW_wd << ' ' << TotalRespTime_RAW_XOR << ' ' << TotalRespTime_RAW_wd_XOR << endl;
		cout << benches[bench_id] << ' ' << TotalTime_DRAM << ' ' << TotalTime_Path << ' ' << TotalTime_RAW << ' ' << TotalTime_RAW_wd << ' ' << TotalTime_RAW_XOR << ' ' << TotalTime_RAW_wd_XOR << ' ' << TotalRespTime_Path << ' ' << TotalRespTime_RAW << ' ' << TotalRespTime_RAW_wd << ' ' << TotalRespTime_RAW_XOR << ' ' << TotalRespTime_RAW_wd_XOR << endl;
		fout.close();
	}
}	


void RAWBurstORAMTest(double PinBW)
{
	// RAW ORAM configuration
	BinPathORam RAWORam;
	
	workingSetInMB = 32 * 1024 * 1024; // 32TB
	oramSizeInMB = 2 * workingSetInMB;
	dataBlockSize = 4096;
		
	RAWORam.Configure(workingSetInMB * (int64_t) 1024 * 1024, oramSizeInMB * (int64_t) 1024 * 1024, dataBlockSize, dataZ);
	int64_t RAWLevelCount = RAWORam.GetLevelCount();
	int64_t ZplusA = RAWORam.GetBlocksPerBucket() + RAWA;
	int64_t ZplusAplusX = ZplusA + ExtraDummy;
	
	RAWLevelCount = RAWLevelCount - LogBurstLength;		// tree top caching	
	
	int64_t  Blk_Lat = dataBlockSize * 8 / PinBW;
	cout << "Blk_Lat = " << Blk_Lat << endl;
	
	uint64_t TraceSegSize = 32768;
	uint64_t TotalAccessCount = 88322452;
	int64_t * req_size = new int64_t [TraceSegSize];
	int64_t * req_time = new int64_t [TraceSegSize];
	
	int64_t * resp_time_raw = new int64_t [TotalAccessCount];
	int64_t * resp_time_opt = new int64_t [TotalAccessCount];
	ifstream traceFile;

	ofstream fout;
	fout.open("netapp_trace_sim.out");	
	
	// netapp trace simulation
	
	// read trace
	string traceFileName = string("input/data_merged_0925_1009");
	traceFile.open(traceFileName.c_str());
	char tmp[255];
	traceFile.getline(tmp, 255);
				
	// statistics to gather	
	int64_t TotalTime_RAW = 0, TotalTime_Opt = 0;
	int64_t TotalRespTime_RAW = 0, TotalRespTime_Opt = 0;
	int64_t ROAccessCount = 0, RWAccessCount = 0;
	int64_t RWBktCount = RAWLevelCount * (dataZ + ZplusAplusX), RWBktCountMax = RAWLevelCount * (dataZ + ZplusAplusX); 
    bool outstandingRWAccess = false;
    
    int64_t thisTraceSegSize = TraceSegSize;
    while (thisTraceSegSize == TraceSegSize)
    {
    	// read in some traces
    	thisTraceSegSize = 0;
    	assert(traceFile.good());
    	for (uint64_t i = 0; i < TraceSegSize && !traceFile.eof(); i++, thisTraceSegSize++)
    	{
    		//if (i == 0) cout << "reading more traces from line " << ROAccessCount << endl;
    		assert(traceFile.good());
    		int64_t offset; // useless
			traceFile >> tmp[0] >> tmp[1] >> req_time[i] >> tmp[2] >> req_size[i] >> tmp[3] >> offset;
			if (req_size[i] <= 0)
				i--, thisTraceSegSize--;
			// cout  << ROAccessCount + i << ' ' << traceFile.eof() << tmp[0] << ' '<< req_time[i] << ' ' << req_size[i] << ' ' << offset << endl;					
		}
		if (traceFile.eof())
			thisTraceSegSize--;	
	
		for (uint64_t i = 0; i < thisTraceSegSize; i++)
		{
			// RAW ORAM with delayed eviction
			int64_t idle_time = req_time[i] - TotalTime_RAW;	// use this idle period to do RW access
			while (idle_time > 0 && outstandingRWAccess)
			{ 
			
				int64_t rw_blk;
				rw_blk = divceil(idle_time, Blk_Lat);
				//rw_blk = idle_time / Blk_Lat;
				rw_blk = min(rw_blk, RWBktCountMax - RWBktCount);
				RWBktCount += rw_blk;
				idle_time -= rw_blk * Blk_Lat;
				
				if (RWBktCount == RWBktCountMax)	// finished 1 RW
				{
					RWAccessCount++;
					RWBktCount = 0;
					outstandingRWAccess = ROAccessCount - RWAccessCount * RAWA > 0;
				}
				
			}
			TotalTime_RAW = max(req_time[i], req_time[i] - idle_time);					
			if (ROAccessCount - RWAccessCount * RAWA >= ((int64_t) RAWA) << LogBurstLength)	// too many delayed RWs, need to perform 1 RW
			{
				TotalTime_RAW += Blk_Lat * (RWBktCountMax-RWBktCount);
				RWAccessCount++;
				RWBktCount = 0;
				outstandingRWAccess = ROAccessCount - RWAccessCount * RAWA > 0;
			}
			
			TotalTime_Opt = max(req_time[i], TotalTime_Opt);
				
			// now we can serve this request, at t = TotalTime		
			int64_t this_req_size = divceil(req_size[i], 8 * dataBlockSize);	// req_size is in bits ! 
			assert(this_req_size > 0);	
			
						
			for (int j = 0; j < this_req_size; j++)
			{
				// RAW ORAM
				TotalTime_RAW += Blk_Lat;	// online latency
				resp_time_raw[ROAccessCount] = TotalTime_RAW - req_time[i];	
				TotalRespTime_RAW += resp_time_raw[i];	
				TotalTime_RAW += Blk_Lat * (dataZ + ZplusAplusX) * ModelEarlyReshuffle(RAWEarlyReshuffleRate / RAWA, RAWLevelCount);	// early reshuffle
				
				// Optimal, no ORAM	
				TotalTime_Opt += Blk_Lat;	// online latency			
				resp_time_opt[ROAccessCount] = TotalTime_Opt - req_time[i];	
				TotalRespTime_Opt += resp_time_opt[i];			
				ROAccessCount++;		
			}
			
			outstandingRWAccess = ROAccessCount - RWAccessCount * RAWA > 0;
		}		
	}
	
	double resp_time_90, resp_time_95, resp_time_99, resp_time_999;
	vector<int64_t> resp_time_sort_vec;
	
	resp_time_sort_vec.assign(resp_time_raw, resp_time_raw + ROAccessCount);
	sort (resp_time_sort_vec.begin(), resp_time_sort_vec.end());  
	  
	resp_time_90 = 50 + resp_time_sort_vec[ROAccessCount * .90] / 1000000.0;
	resp_time_95 = 50 + resp_time_sort_vec[ROAccessCount * .95] / 1000000.0;   
	resp_time_99 = 50 + resp_time_sort_vec[ROAccessCount * .99] / 1000000.0;
	resp_time_999 = 50 + resp_time_sort_vec[ROAccessCount * .999] / 1000000.0;
	cout << PinBW << ' ' << resp_time_90 << ' ' << resp_time_95 << ' ' << resp_time_99 << ' ' << resp_time_999 << endl;
/*
	resp_time_sort_vec.assign(resp_time_opt, resp_time_opt + ROAccessCount);
	sort (resp_time_sort_vec.begin(), resp_time_sort_vec.end());  
	resp_time_90 = 50 + resp_time_sort_vec[ROAccessCount * .90] / 1000000.0;
	resp_time_95 = 50 + resp_time_sort_vec[ROAccessCount * .95] / 1000000.0;   
	resp_time_99 = 50 + resp_time_sort_vec[ROAccessCount * .99] / 1000000.0;
	resp_time_999 = 50 + resp_time_sort_vec[ROAccessCount * .999] / 1000000.0;
	cout << PinBW << ' ' << resp_time_90 << ' ' << resp_time_95 << ' ' << resp_time_99 << ' ' << resp_time_999 << endl;
*/	
	fout.close();
	
	resp_time_sort_vec.clear();
	delete []req_size;
	delete []req_time;
	delete []resp_time_raw;
	delete []resp_time_opt;
}	


void HierORAMTest()
{
	uint64_t workingSet = workingSetInMB * ((uint64_t) 1024) * 1024;
	double ORAMUtil[10] = {workingSetInMB * 1.0 / oramSizeInMB, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
	int HORamBlockSize[10] = {dataBlockSize};
	int HORamBlocksPerBucket[10] = {dataZ};	 
	for (int i = 1; i < 10; i++)
		HORamBlockSize[i] = posBlockSize, HORamBlocksPerBucket[i] = posZ;

	int AES_pad_in_bits = 8 * 8;		//
	int chunk_size_in_bytes = 64;		// should be multiple of DRAM access width

	int bucket_size_in_chunks[10];	

	HierBinPathORam HORam;
	HORam.Configure(workingSet, ORAMUtil, HORamBlockSize, HORamBlocksPerBucket, onchipPosMapSize);

	fout.open(outFileName.c_str(), ios_base::out);
	cout << "Final position map size = " << HORam.GetFinalPositionMapSize() / 1024 << " KB" << endl;
	fout<<"Progress	Useful / Total"<<endl;
	fout.close();	

	HORam.Initialize();
	for (int i = 0; i < HORam.GetHierarchy(); i++)
		HORam.ORam[i]->SetMaxStashSize(maxLCSZ);

	int64_t BlockID;
	uint64_t validBlockCount = HORam.ORam[0]->GetValidBlockCount();
	assert(RAND_MAX > validBlockCount);
	uint64_t RecordInterval = validBlockCount >> 2;

	for (uint64_t k = 0; k < validBlockCount; k++)
	{
		HORam.Access(k, -1, BinPathORam::write);
		HORam.BackgroundEvict();
		if (k % RecordInterval == RecordInterval - 1)			
		{
			fout.open(outFileName.c_str(), ios_base::app);				
			fout<< k / RecordInterval << '\t' << HORam.GetNumAccess() << '\t' << HORam.GetNumDummy() << '\t' << k * 1.0 / HORam.GetNumAccess() <<endl;
			fout.close();
		}
	}
	
	fout.open(outFileName.c_str(), ios_base::app);	
	fout << endl;	
	fout.close();
	
	uint64_t Total = runlength * validBlockCount;
	for (uint64_t k = validBlockCount; k < Total; k++)
	{
		BlockID = rand() % validBlockCount;
		if (worstCase)
			BlockID = k % validBlockCount;
		HORam.Access(BlockID, -1, BinPathORam::read);
		HORam.BackgroundEvict();
		if (k % RecordInterval == RecordInterval - 1)			
		{
			fout.open(outFileName.c_str(), ios_base::app);				
			fout<< k / RecordInterval << '\t' << HORam.GetNumAccess() << '\t' << HORam.GetNumDummy() << '\t' << k * 1.0 / HORam.GetNumAccess() <<endl;
			fout.close();
		}
	}
	fout.open(outFileName.c_str(), ios_base::app);	
	fout << endl << Total * 1.0 / HORam.GetNumAccess() <<endl;
	fout.close();
}

void SingleORAMTest()
{
	BinPathORam ORam;
	ORam.Configure(workingSetInMB * (int64_t) 1024 * 1024, oramSizeInMB * (int64_t) 1024 * 1024, dataBlockSize, dataZ);
	ORam.Initialize();

	int64_t BlockID;
	int validBlockCount = ORam.GetValidBlockCount();
	assert(RAND_MAX > validBlockCount);

	ORam.SetMaxStashSize(maxLCSZ);
	int RecordInterval = validBlockCount >> 2;

	fout.open(outFileName.c_str(), ios_base::out);	
	fout<<"Progress	Rate / Total"<<endl;
	fout.close();	

	for (uint64_t k = 0; k < validBlockCount; k++)
	{
		ORam.Access(k, -1, BinPathORam::write, NULL);
		ORam.BackgroundEvict();
		if (k % RecordInterval == RecordInterval - 1)
		{	
			fout.open(outFileName.c_str(), ios_base::app);	
			fout<< k / RecordInterval << '\t' << ORam.GetNumAccess() << '\t' << ORam.GetNumDummy() << '\t' << k * 1.0 / ORam.GetNumAccess() <<endl;
			fout.close();
		}
	}

	fout.open(outFileName.c_str(), ios_base::app);	
	fout << endl;	
	fout.close();

	int64_t total_traffic = 0, total_levels = 0;

	ORam.EnableHistogram(1000, 2);

	uint64_t Total = runlength * validBlockCount;
	for (uint64_t k = validBlockCount; k < Total; k++)
	{
		BlockID = rand() % validBlockCount;
		if (worstCase)
			BlockID = k % validBlockCount;
	
		total_traffic += ORam.BackgroundEvict();			
		total_traffic += ORam.Access(BlockID, -1, BinPathORam::read, NULL);
		ORam.RecordHistogram();	
		
		total_levels += ORam.GetCurPathLength();	
		
		if (k % RecordInterval == RecordInterval - 1)
		{
			fout.open(outFileName.c_str(), ios_base::app);				
			fout<< (k - validBlockCount) / RecordInterval << '\t' <<  total_traffic * 1.0 / (k - validBlockCount) <<endl;
			fout.close();
		}
	}
	ORam.DumpHistogram("PathLoad.txt");

	double ave_traffic = total_traffic * 1.0 / (Total - validBlockCount);
	fout.open(outFileName.c_str(), ios_base::app);	
	fout << "Access overhead: \t" << ave_traffic << "\t\t"<< ave_traffic / total_levels * (Total - validBlockCount) << endl;
	fout.close();
	cout << "Access overhead: \t" << ave_traffic << "\t\t"<< ave_traffic / total_levels * (Total - validBlockCount) << endl;
	
	cout<< Total / RecordInterval << '\t' << ORam.GetNumAccess() << '\t' << ORam.GetNumDummy() << '\t' << Total * 1.0 / ORam.GetNumAccess() <<endl;
	
}

inline uint64_t ReverseBit(uint64_t input, int n)
{
	uint64_t output;
	for (int i = 0; i < n; i++)
	{
		output <<= 1;
		output += input & 1;
		input >>= 1;
	}
	return output;
}

void RAWORAMTest()
{
	BinPathORam ORam;
	// ORam.Configure(workingSetInMB * (int64_t) 1024 * 1024, oramSizeInMB * (int64_t) 1024 * 1024, dataBlockSize, dataZ);
	int ORamLevels = levels;
	ORam.Configure((1 << ORamLevels) / 4 * RAWA, ORamLevels, dataZ);	// (1 << ORamLevels) / 4 * RAWA, 0r 25% util is what's proven.
	ORam.Initialize();

	int64_t BlockID;
	int validBlockCount = ORam.GetValidBlockCount();
	assert(RAND_MAX > validBlockCount);

	// ORam.SetMaxStashSize(maxLCSZ);
	int RecordInterval = validBlockCount >> 2;

	double EvictRate = (double) RAWA;
	uint64_t EvictCounter = 0;
	uint64_t ROCounter = 0;
	int LeafBase = 1 << (ORamLevels - 1);
	LeafBase --;
	double AveStashOccupancy = 0;
	int64_t total_traffic = 0, total_levels = 0;
	uint64_t Total = runlength * validBlockCount;

	ORam.EnableHistogram(1000, 2);
	for (uint64_t k = 0; k < 2 * validBlockCount; k++)
	{
		ORam.AccessOneBlock(k % validBlockCount, -1, BinPathORam::write);
		ROCounter++;
		if (ROCounter > EvictRate * EvictCounter)
		{		
			uint64_t EvictPath = ReverseBit(EvictCounter, ORamLevels - 1);
			assert(EvictPath < LeafBase + 1);
			EvictPath += LeafBase;
			EvictCounter++;
			ORam.Access(-1, EvictPath, BinPathORam::dummy, NULL);
		}
		
	}

	for (uint64_t k = validBlockCount; k < Total; k++)
	{
		BlockID = rand() % validBlockCount;
		if (worstCase)
			BlockID = k % validBlockCount;

		total_traffic += ORam.AccessOneBlock(BlockID, -1, BinPathORam::read);
		total_levels += ORam.GetCurPathLength();	
		
		ROCounter++;
		while (ROCounter > EvictRate * EvictCounter)
		{						
			uint64_t EvictPath = ReverseBit(EvictCounter, ORamLevels - 1);
			assert(EvictPath < LeafBase + 1);
			EvictPath += LeafBase;
			EvictCounter++;
			
			//total_traffic += ORam.Access(-1, EvictPath, BinPathORam::dummy, NULL);
			total_traffic += ORam.ForegroundEvict1(EvictPath);
			
			AveStashOccupancy += ORam.GetCurLocalCacheSize();
			ORam.RecordHistogram();		
		}	
		
	}

	ORam.DumpHistogram("PathLoad.txt");

	double ave_traffic = total_traffic * 1.0 / (Total - validBlockCount);
	AveStashOccupancy /= (Total - validBlockCount) / EvictRate;
	cout << "Access overhead: " << ave_traffic << "\t"<< ave_traffic / total_levels * (Total - validBlockCount) << "\nMaximum/Ave Stash: " << ORam.GetPeakLocalCacheSize() << '\t' << AveStashOccupancy << endl;
}


void BESecureTest()
{
/*	int workingSetInMB = 16;
	int oramSizeInMB = 63;
	int dataBlockSize = 1;	
	int dataZ = 1;
	int superBlockSize = 1;
	int maxLCSZ = 2;
	int runlength = 100;
	bool worstCase = 1;
*/
	srand (time(NULL));
	BinPathORam ORam;
	ORam.Configure(workingSetInMB, oramSizeInMB, dataBlockSize, dataZ);
	ORam.Initialize();

	int BlockID;
	int blockCount = ORam.GetValidBlockCount();
	assert(RAND_MAX > blockCount);

	fout.open("./common_path.txt", ios_base::app);

	ORam.SetMaxStashSize(maxLCSZ);
	int RecordInterval = blockCount >> 2;
	for (int k = 0; k < blockCount; k++)
	{
		ORam.Access(k, -1, BinPathORam::write, NULL);
		ORam.InsecureBackgroundEvict();
//		ORam.BackgroundEvict();
	}	
	uint64_t Total = runlength * ORam.GetValidBlockCount();
	for (uint64_t k = blockCount; k < Total; k++)
	{
		BlockID = rand() % blockCount;
		if (worstCase)
			BlockID = k % blockCount;
		ORam.Access(BlockID, -1, BinPathORam::read, NULL);
		ORam.InsecureBackgroundEvict();
//		ORam.BackgroundEvict();
	}
	uint64_t sum = 0;
//	for (int i = 1; i < ORam.ServerSees.size(); i++)
//		sum += ORam.GetLevelCount() - HighestBit(ORam.ServerSees[i] ^ ORam.ServerSees[i-1]);
	for (int i = 0; i < ORam.ServerSees.size(); i++)
		for (int j = i+1; j < i + 2 * workingSetInMB + 1&& i < ORam.ServerSees.size(); j++)
			sum += ORam.GetLevelCount() - HighestBit(ORam.ServerSees[i] ^ ORam.ServerSees[j]);
	cout << sum * 1.0 / ORam.ServerSees.size() / workingSetInMB / 2 << endl;
	fout << sum * 1.0 / ORam.ServerSees.size() / workingSetInMB / 2 << endl;
}

int LocalCacheHist(const char * filename)
{
	srand (time(NULL));
	BinPathORam ORam;
	if (!levels)
		ORam.Configure(workingSetInMB * (uint64_t) 1024 * 1024, oramSizeInMB * (uint64_t) 1024 * 1024, dataBlockSize, dataZ);
	else
		ORam.Configure(workingSetInMB, levels, dataZ);
	ORam.Initialize();

	int BlockID;
	int validBlockCount = ORam.GetValidBlockCount();
	assert(RAND_MAX > validBlockCount);

	int RecordInterval = validBlockCount << 3;
	for (int k = 0; k < validBlockCount; k++)
	{
		ORam.Access(k, -1, BinPathORam::write, NULL);
	}	
	ORam.EnableHistogram(histSize, 2);
	uint64_t Total = runlength * validBlockCount;
	for (uint64_t k = validBlockCount; k < Total; k++)
	{
		BlockID = rand() % validBlockCount;
		if (worstCase)
			BlockID = k % validBlockCount;
		ORam.Access(BlockID, -1, BinPathORam::read, NULL);
		if (k % RecordInterval == RecordInterval - 1)
			ORam.DumpHistogram(filename);
	}
	return ORam.GetPeakLocalCacheSize();
}

void SimDRAMLatency()
{
#ifdef USE_DRAMSIM2	
	uint64_t workingSet = workingSetInMB * ((uint64_t) 1024) * 1024;
	double ORAMUtil[10] = {workingSetInMB * 1.0 / oramSizeInMB, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
	int HORamBlockSize[10] = {dataBlockSize};
	int HORamBlocksPerBucket[10] = {dataZ};	 
	for (int i = 1; i < 10; i++)
		HORamBlockSize[i] = posBlockSize, HORamBlocksPerBucket[i] = posZ;

	HierBinPathORam HORam;
	HORam.Configure(workingSet, ORAMUtil, HORamBlockSize, HORamBlocksPerBucket, 100*1024);	//FINAL_POS_MAP_SIZE
//	cout << "Final position map size = " << HORam.GetFinalPositionMapSize() / 1024 << " KB" << endl;

	int AES_pad_in_bits = 8 * 8;
	char Dir[30] = "./";
	HORam.DRAMConfigure(AES_pad_in_bits, 1.0, Dir, (char*) "DDR3_micron_16M_8B_x8_sg15_ORAM.ini", (char*) "system_ORAM.ini");
	if (!runlength)	return;
	char strategies[10] = {7};
	for (int i = 0; i < 1; i++)
	{
		HORam.Simulator->Strategy = strategies[i];
		HORam.Simulator->Reset();
		HORam.SimulateLatency(runlength);	
	}		
	
	fout.open("HORAM_overhead_breakdown", ios_base::app);
	char tab = '\t';
	fout << workingSetInMB << tab << oramSizeInMB << tab << dataZ << tab << posZ << tab << posBlockSize << tab << integrity << tab;
	fout << HORam.GetHitDelay() << tab << HORam.GetReadyDelay();
	fout << endl;
	fout.close();
#endif	
}

void Miscellany()
{
	BinPathORam ORam;
	
	workingSetInMB = 15;
	oramSizeInMB = 21;
	dataBlockSize = 1;
	dataZ = 1;
	ORam.Configure(workingSetInMB * (int64_t) 1024 * 1024, oramSizeInMB * (int64_t) 1024 * 1024, dataBlockSize, dataZ);
	//ORam.Configure(workingSetInMB, oramSizeInMB, dataBlockSize, dataZ);

	int x, y, z;
	for (int i = 0; i < 20; i++)
	{
		x = ORam.GenerateRandLeaf();
		y = ORam.GenerateRandLeaf();
		ORam.CurPathLength = HighestBit(y + 1);
		z = ORam.PosDiff(x, y);
		cout << x << ' ' << y << ' ' << ORam.CurPathLength - z - 1 << endl;
	}
}




