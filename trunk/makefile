#$Id$
.PHONY: dos 

markov: markov.cpp
	g++ markov.cpp -o markov

#make the files readable by KGraph
dos: 
	unix2dos *.dat

test: markov.cpp test_suite
	g++ markov.cpp -o markov
	./test_suite
	touch test

clean:
	rm test
	rm *dat
