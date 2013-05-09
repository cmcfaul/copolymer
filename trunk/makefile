#$Id$
.PHONY: dos 

det: deterministic.cpp
	g++ deterministic.cpp -o deterministic

markov: markov.cpp
	g++ markov.cpp -o markov

#make the files readable by KGraph
dos: 
	unix2dos *.dat

test: deterministic.cpp test_suite
	g++ deterministic.cpp -o deterministic
	./test_suite
	touch test

clean:
	rm test
	rm *dat
