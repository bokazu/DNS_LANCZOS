gcc_options = -std=c++17 -Wall --pedantic-errors -DMKL_ILP64  -I"${MKLROOT}/include" -g
l_b = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

program : main.o dns_eigen.o dns_lanczos.o func.o
	g++ -o $@ $^ $(l_b)

mainv2.o : main.cpp
	g++ -c $(gcc_options) $< $(l_b)

dns_eigen.o : ./dns_lanczos/dns_eigen.cpp
	g++ -c $(gcc_options) $< $(l_b)

dns_lanczos.o : ./dns_lanczos/dns_lanczos.cpp
	g++ -c $(gcc_options) $< $(l_b)

func.o : ./dns_lanczos/func.cpp
	g++ -c $(gcc_options) $< $(l_b)

run : program
	./program

clean:
	rm -f ./program

.PHONY : run clean
