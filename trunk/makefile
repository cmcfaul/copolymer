#$Id$

markov: markov.cpp
	g++ markov.cpp -o markov

#make the files readable by KGraph
dos: markov test.dat test_stats.dat
	unix2dos *.dat
