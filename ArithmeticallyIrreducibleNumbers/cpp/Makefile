CXX=clang++ 
FLAGS=-std=c++23 -O3 -fopenmp -ffast-math -Wall -Werror -Wextra

main : as_irreducible.cpp 
	$(CXX) $(FLAGS) as_irreducible.cpp -o as_irreducible   

dbg : as_irreducible.cpp 
	$(CXX) -g -fsanitize=address -Wall -Werror -Wextra as_irreducible.cpp -o as_irreducible_g
	
clean : 
	rm -f as_irreducible as_irreducible_g
