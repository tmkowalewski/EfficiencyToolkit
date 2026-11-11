CXX ?= c++
CPPFLAGS += -Iinclude $(shell root-config --cflags)
CXXFLAGS ?= -O2 -Wall -Wextra
LDLIBS += $(shell root-config --libs) -lboost_program_options

.PHONY: all clean

all: make_eff_curves get_sim_efficiencies

make_eff_curves: src/make_eff_curves.cpp src/efficiency_plots.cpp include/efficiency_plots.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/make_eff_curves.cpp src/efficiency_plots.cpp $(LDFLAGS) $(LDLIBS) -o $@

get_sim_efficiencies: src/get_sim_efficiencies.cpp src/simulation.cpp include/simulation.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/get_sim_efficiencies.cpp src/simulation.cpp $(LDFLAGS) $(LDLIBS) -o $@

compare_layout: tests/compare_layout.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) make_eff_curves get_sim_efficiencies compare_layout
