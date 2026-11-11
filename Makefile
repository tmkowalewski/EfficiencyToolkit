CXX ?= c++
CPPFLAGS += -Iinclude $(shell root-config --cflags)
CXXFLAGS ?= -O2 -Wall -Wextra
LDLIBS += $(shell root-config --libs) -lboost_program_options

.PHONY: all clean make_eff_curves get_sim_efficiencies compare_layout

all: bin/make_eff_curves bin/get_sim_efficiencies

make_eff_curves: bin/make_eff_curves

get_sim_efficiencies: bin/get_sim_efficiencies

compare_layout: bin/compare_layout

bin:
	mkdir -p $@

bin/make_eff_curves: src/make_eff_curves.cpp src/efficiency_plots.cpp include/efficiency_plots.hpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/make_eff_curves.cpp src/efficiency_plots.cpp $(LDFLAGS) $(LDLIBS) -o $@

bin/get_sim_efficiencies: src/get_sim_efficiencies.cpp src/simulation.cpp include/simulation.hpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/get_sim_efficiencies.cpp src/simulation.cpp $(LDFLAGS) $(LDLIBS) -o $@

bin/compare_layout: tests/compare_layout.cpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) bin/make_eff_curves bin/get_sim_efficiencies bin/compare_layout
