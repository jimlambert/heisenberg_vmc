# Simple Makefile
MAIN = HEIS_VMC
AUX = AuxiliaryHamiltonians
GEO = Geoms
HIL = Hilbert
OBS = Observables
PAR = Parameters
UTL = Utils
WAV = Wavefunction
CXX = g++ 
CXX_FLAGS = -std=c++17 -O2 -Wall -fPIC -Wextra -lboost_program_options -lboost_system
CPP_FILES = $(wildcard src/*.cpp) \
					  $(wildcard src/$(AUX)/*.cpp) \
					  $(wildcard src/$(GEO)/*.cpp) \
					  $(wildcard src/$(HIL)/*.cpp) \
						$(wildcard src/$(OBS)/*.cpp) \
						$(wildcard src/$(PAR)/*.cpp) \
					  $(wildcard src/$(UTL)/*.cpp) \
					  $(wildcard src/$(WAV)/*.cpp)	
OBJ_FILES = $(addprefix obj/, $(notdir $(CPP_FILES:.cpp=.o)))
INCLUDE = -I/usr/include/eigen3/ -I./inc/$(AUX)/ -I./inc/$(GEO)/ -I./inc/$(HIL)/ -I./inc/$(OBS)/ -I./inc/$(PAR)/ -I./inc/$(UTL)/ -I./inc/$(WAV)/ 

$(MAIN): $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -o $@ $^

obj/%.o: src/%.cpp  
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(AUX)/%.cpp inc/$(AUX)/%.h
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(GEO)/%.cpp inc/$(GEO)/%.h
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(HIL)/%.cpp inc/$(HIL)/%.h
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(OBS)/%.cpp inc/$(OBS)/%.h 
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(PAR)/%.cpp inc/$(PAR)/%.h
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(UTL)/%.cpp inc/$(UTL)/%.h  
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/%.o: src/$(WAV)/%.cpp inc/$(WAV)/%.h 
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

clean:	
	$(RM) $(OBJ_FILES) $(MAIN)


