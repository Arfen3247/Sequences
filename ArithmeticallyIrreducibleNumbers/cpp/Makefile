CXX=g++ -std=c++23
FLAGS=-O3 -Wall -Werror -Wextra

main : AS_irreducible.cpp 
	$(CXX) $(FLAGS) AS_irreducible.cpp -o AS_irreducible   

dbg : AS_irreducible.cpp 
	$(CXX) -g -fsanitize=address -Wall -Werror -Wextra AS_irreducible.cpp -o AS_irreducible_g
	
clean : 
	rm -f AS_irreducible AS_irreducible_g
