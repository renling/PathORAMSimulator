CXX = g++ 
DRAMSim2 = ../../DRAMSim2
SRC = BinPathORam_test.cpp BinPathORam.cpp HierBinPathORam.cpp UnifiedBinPathORam.cpp	SimpleCache.cpp
PARAM = -O3 -I. -I$(DRAMSim2) -L$(DRAMSim2) -ldramsim -Wl,-rpath=../../DRAMSim2

all: BinPathORam_test 
BinPathORam_test: $(SRC) *.h
	$(CXX) -o $@ $(SRC) $(PARAM)
clean:
	rm BinPathORam_test
