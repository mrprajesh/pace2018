all:
	g++ -std=c++14 -Wall -o "st-heuristic" "st-heuristic.cpp" -Ofast
	g++ -std=c++14 -Wall -o "st-exact" "st-exact.cpp" -Ofast
