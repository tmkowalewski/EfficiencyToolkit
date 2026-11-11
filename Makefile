CXX ?= c++
CPPFLAGS += -Iinclude $(shell root-config --cflags)
CXXFLAGS ?= -O2 -Wall -Wextra
LDLIBS += $(shell root-config --libs) -lboost_program_options

.PHONY: all clean make_eff_curves get_sim_efficiencies compare_layout test

all: bin/make_eff_curves bin/get_sim_efficiencies

make_eff_curves: bin/make_eff_curves

get_sim_efficiencies: bin/get_sim_efficiencies

compare_layout: bin/compare_layout

test: bin/make_eff_curves bin/check_absolute bin/check_radware bin/check_plot_output
	bash tests/test_absolute.sh

bin:
	mkdir -p $@

bin/make_eff_curves: src/make_eff_curves.cpp src/efficiency_plots.cpp src/radware_fit.cpp include/efficiency_plots.hpp include/radware_fit.hpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/make_eff_curves.cpp src/efficiency_plots.cpp src/radware_fit.cpp $(LDFLAGS) $(LDLIBS) -o $@

bin/get_sim_efficiencies: src/get_sim_efficiencies.cpp src/simulation.cpp include/simulation.hpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/get_sim_efficiencies.cpp src/simulation.cpp $(LDFLAGS) $(LDLIBS) -o $@

bin/compare_layout: tests/compare_layout.cpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

bin/check_absolute: tests/check_absolute.cpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

bin/check_radware: tests/check_radware.cpp src/radware_fit.cpp include/radware_fit.hpp include/efficiency_plots.hpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/check_radware.cpp src/radware_fit.cpp $(LDFLAGS) $(LDLIBS) -o $@

bin/check_plot_output: tests/check_plot_output.cpp | bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) bin/make_eff_curves bin/get_sim_efficiencies bin/compare_layout bin/check_absolute bin/check_radware bin/check_plot_output
