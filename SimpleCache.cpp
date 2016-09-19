#include "SimpleCache.h"

SimpleCache::SimpleCache(int capacity, int blockSize, int ways)
{
	Capacity = capacity;
	BlockSize =  blockSize;
	Ways = ways;
	
	Blocks = Capacity / BlockSize;
	Sets = Blocks / Ways;
	
	cache = new SimpleCacheline [Blocks];		
	evicted = new SimpleCacheline;

	numAccess = 0;	
	debug = 0;
}

SimpleCache::~SimpleCache() 
{ 
	delete []cache;
	delete evicted;
}

int SimpleCache::Read(uint64_t addr)
{
	numAccess++;
	int set = addr % Sets;
	int hit_at = -1;
	for (int j = 0; j < Ways; j++)
	{
		SimpleCacheline * current = cache + set * Ways + j;
		if (current->valid && current->tag == addr)
		{
			current->last_time = numAccess;
			hit_at = set * Ways + j;
			break;
		}
	}
	
	if (debug)
		printf("read %ld, result = %d\n", addr, hit_at);
	
	return hit_at;
}

int SimpleCache::Write(uint64_t addr)
{	
	assert(Read(addr) < 0);
	numAccess--;		
	int set = addr % Sets;
	uint64_t least_recent_time = cache[set * Ways].last_time;
	int replacement = 0;
	SimpleCacheline * current;
	for (int j = 0; j < Ways; j++)
	{
		current = cache + set * Ways + j;
		if (!current->valid)
		{
			replacement = j;
			break;
		}
		else if (current->last_time < least_recent_time)
		{
			replacement = j;
			least_recent_time = current->last_time;
		}
	}
	current = cache + set * Ways + replacement;
	*evicted = *current;
	*current = (SimpleCacheline) {1, addr, numAccess};	
		
	if (debug)
	{
		printf("write %ld, ", addr);
		if (evicted->valid)
			printf("evicting %ld\n", evicted->tag);
		else
			printf("no eviction\n");
	}	
	
	return evicted->valid;		// -1 no eviction; -2 evicted some other line
}
