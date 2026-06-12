CXX=clang++
FLAGS= -std=c++23 -O3 -fopenmp -ffast-math -Wall -Werror -Wextra

main : count.cpp 
	$(CXX) $(FLAGS) count.cpp -o count   

dbg : count.cpp 
	$(CXX) -g -fsanitize=address -Wall -Werror -Wextra count.cpp -o count_g
	
clean : 
	rm -f count count_g
