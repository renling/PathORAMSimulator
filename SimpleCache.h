#include <stdint.h>
#include <cstdio>

#define NDEBUG
#undef NDEBUG	// checking asserts
#include <assert.h>

class SimpleCacheline
{
public:
	SimpleCacheline(int valid = 0, uint64_t tag = 0, uint64_t last_time = 0)
	{
		this->valid = valid;
		this->tag = tag;
		this->last_time = last_time;
	}
	
	bool valid;
	uint64_t tag;
	uint64_t last_time;
};

class SimpleCache
{
public:
	SimpleCache(int capacity, int blockSize, int ways);
	~SimpleCache();
	int Read(uint64_t addr);
	int Write(uint64_t addr);
	
	SimpleCacheline * evicted;
private:
	int Capacity;	// in Bytes, power of 2
	int BlockSize;
	int Ways;
	int Blocks;
	int Sets;	

	SimpleCacheline * cache;
	
	uint64_t numAccess;
	short debug;
};
