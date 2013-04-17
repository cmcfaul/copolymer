#$Id$

markov: markov.cpp
	g++ markov.cpp -o markov

#make the files readable by KGraph
dos: markov *.dat 
	unix2dos *.dat
