
#include "Solver.h"


#include <iostream>
#include <string>


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*                                                                           */
/*		Questions, remarks													 */
/*																			 */
/*		- In which points the surface tension term is discretized? at (i,j)	 */
/*		  or at	(i+1/2,j)?													 */
/*		- In the term (U)dot(DivU) should be used upwind? (nonsense...)		 */
/*		- How the level set is initialized?									 */
/*		- Boundary conditions for Poisson solver?							 */
/*																			 */
/*																			 */
/*		- Use higher-order in time and space								 */
/*		- In Surface tension I return (0,0) when norm == 0 (singularity at	 */
/*		  the origin of the coordinates)									 */
/*																			 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

using namespace Solver;


/*****************************************************************************/
/*                                                                           */
/*		Set the domain dimensions : [ax, bx] x [ay, by]						 */
/*																			 */
/*		[ax, bx] = Interval in x-direction									 */
/*		[ay, by] = Interval in y-direction									 */
/*																			 */
/*****************************************************************************/


/************************************************/
/*                                              */
/*		Drive-lid cavity						*/
/*												*/
/************************************************/
double const ax = 0.0;
double const bx = 1.0;

double const ay = 0.0;
double const by = 1.0;

/************************************************/
/*                                              */
/*		Rising bubble, falling droplet			*/
/*												*/
/************************************************/
//double const ax = 0.0;
//double const bx = 0.12;
//
//double const ay = 0.0;
//double const by = 0.11;

//double const ax = 0.0;
//double const bx = 1.0;
//double const ay = 0.0;
//double const by = 1.0;
//int const refine = 64;
//unsigned const nx = refine * 16;
//unsigned const ny = nx;

/*****************************************************************************/
/*                                                                           */
/*		Set the number of grid points										 */
/*																			 */
/*		nx = Number of grid points in x-direction							 */
/*		ny = Number of grid points in y-direction							 */
/*																			 */
/*****************************************************************************/
unsigned const nx = 80;
unsigned const ny = nx;

int const NT0	= 0;
int const NT	= 800;


int main() {

	

	Solver::IncompressibleTwoPhaseFlow Solution(nx, ax, bx, ny, ay, by);

	//Solution.EOC();
	//return 0;


	/************************************************/
	/*                                              */
	/*		Drive-lid cavity						*/
	/*												*/
	/************************************************/
	Solution.setTimeIncrement(0.002);
	Solution.setAdaptiveTimeStep(false);

	/************************************************/
	/*                                              */
	/*		Rising bubble, falling droplet			*/
	/*												*/
	/************************************************/
	//Solution.setTimeIncrement(1e0);
	//Solution.setAdaptiveTimeStep(true);


	//Solution.LevelSetExport("C:\\Users\\pgali\\Desktop\\levelset\\levelset" + std::to_string(NT0) + ".txt");
	////Solution.pressureExport("C:\\Users\\pgali\\Desktop\\levelset\\pressure" + std::to_string(NT0) + ".txt");
	//Solution.velocityExport("C:\\Users\\pgali\\Desktop\\levelset\\velocity" + std::to_string(NT0) + ".txt");
	//Solution.vorticityExport("C:\\Users\\pgali\\Desktop\\levelset\\vorticity" + std::to_string(NT0) + ".txt");

	Solution.LevelSetExport("D:\\simulations\\levelset\\levelset" + std::to_string(NT0) + ".txt");
	//Solution.pressureExport("D:\simulations\\levelset\\pressure" + std::to_string(NT0) + ".txt");
	Solution.velocityExport("D:\\simulations\\levelset\\velocity" + std::to_string(NT0) + ".txt");
	Solution.vorticityExport("D:\\simulations\\levelset\\vorticity" + std::to_string(NT0) + ".txt");

	//Solution.meshExport("C:\\Users\\pgali\\Desktop\\levelset\\X.txt", "C:\\Users\\pgali\\Desktop\\levelset\\Y.txt");

	for (int nt = NT0 + 1; nt <= NT; nt++) {

		Solution.setTimeLevel(nt);

		Solution.solve();

		//if (nt % 150 == 0) {
			Solution.LevelSetExport("D:\\simulations\\levelset\\levelset" + std::to_string(nt) + ".txt");
			//Solution.pressureExport("D:\simulations\\levelset\\pressure" + std::to_string(nt) + ".txt");
			Solution.velocityExport("D:\\simulations\\levelset\\velocity" + std::to_string(nt) + ".txt");
			Solution.vorticityExport("D:\\simulations\\levelset\\vorticity" + std::to_string(nt) + ".txt");
		//}
	}

	/*int nt = NT0 + 1;
	while(Solution.getTime()<1.5) {

		Solution.setTimeLevel(nt);
		Solution.solve();
		Solution.LevelSetExport("C:\\Users\\pgali\\Desktop\\levelset\\levelset" + std::to_string(nt) + ".txt");

		nt++;

	}*/

	   	 

	return 0;

}



/*
void generate_points(Grid & grid) {


	
	
	//		Insert inner points													
	
	
	for (unsigned i = 1; i < nx - 1; i++) {

		double const x = ax + i * (bx - ax) / (nx - 1);

		for (unsigned j = 1; j < ny - 1; j++) {

			double const y = ay + j * (by - ay) / (ny - 1);

			grid.insert_point<PointMarker::Inner>(Point(x, y));

		}
	}

	
	//                                                                           
	//		Insert Dirichlet points												 
	//																			
	
	for (unsigned i = 0; i < nx; i++) {

		double const x = ax + i * (bx - ax) / (nx - 1);

		grid.insert_point<PointMarker::Dirichlet>(Point(x, ay));
		grid.insert_point<PointMarker::Dirichlet>(Point(x, by));

	}
	for (unsigned j = 0; j < ny; j++) {

		double const y = ay + j * (by - ay) / (ny - 1);

		grid.insert_point<PointMarker::Dirichlet>(Point(ax, y));
		grid.insert_point<PointMarker::Dirichlet>(Point(bx, y));

	}

}
*/