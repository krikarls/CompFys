
#include <iostream>

void GQ();		// function for both Guassian Quadrature methods
void MC();		// function for both Monte carlo method

int main()
{

	int M;

	std::cout << "For Gaussian quadrature methods enter: 0" << std::endl;
	std::cout << "For Monte carlo methods enter: 1" << std::endl;
	std::cin >> M ;

	if(M == 0) GQ();
	else if(M == 1) MC();
  	else std::cout << "No valid option choosen.";





}