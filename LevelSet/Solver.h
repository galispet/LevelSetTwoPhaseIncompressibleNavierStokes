#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>


#include <Eigen/Sparse>
#include <Eigen/Dense>




namespace Solver {

	
	/************************************************/
	/*                                              */
	/*		Template parameter 						*/
	/*												*/
	/*		0 for Column major						*/
	/*		1 for Row major							*/
	/*												*/
	/************************************************/
	typedef Eigen::SparseMatrix<double, 0> SparseMatrix;



	class IncompressibleTwoPhaseFlow {


	private:


		/************************************************/
		/*                                              */
		/*		Physical quantities						*/
		/*												*/
		/************************************************/
		double const Pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
		
		/************************************************/
		/*                                              */
		/*		Driven-lid cavity						*/
		/*												*/
		/************************************************/
		/*double const gx = 0.0;
		double const gy = 0.0;

		double const rho1 = 1.0;
		double const rho2 = 1.0;

		double const mu1 = 0.00002;
		double const mu2 = 0.00002;
		
		double const sigma = 0.0;*/

		/************************************************/
		/*                                              */
		/*		Rising bubble, falling droplet			*/
		/*												*/
		/************************************************/
		double const gx = 0.0;
		double const gy = 0.0*-9.81;

		double const rho2 = 1.0;
		double const rho1 = 1.0 * 1e1;

		double const mu2 = 2.5 * 1e-4;
		double const mu1 = 5.0 * 1e-4;

		double const sigma = 5.0 * 1e-3;





		/************************************************/
		/*                                              */
		/*		Geometry and time setup					*/
		/*												*/
		/************************************************/
		double dx;
		double dy;
		double dt;
		double dtau;
		int nt;

		bool AdaptiveTimeStep = false;

		

		


		/************************************************/
		/*                                              */
		/*		Smoothing parameter						*/
		/*												*/
		/************************************************/
		double alpha;


		/************************************************/
		/*                                              */
		/*		Sparse matrix for poisson solver. 		*/
		/*												*/
		/************************************************/
		std::vector<int>					ptr;
		std::vector<int>					col;
		std::vector<double>					val;
		std::vector<Eigen::Triplet<double>> triplet;

		SparseMatrix					A;
		Eigen::VectorXd					pressureRightHandSide;
		Eigen::VectorXd					pressureSolution;
		Eigen::BiCGSTAB<SparseMatrix>	poissonSystemSolverBiCGSTAB;


		Eigen::MatrixXd p;

		Eigen::MatrixXd phi;
		Eigen::MatrixXd phiZero;
		Eigen::MatrixXd phiStar;
		Eigen::MatrixXd phiStarStar;
		Eigen::MatrixXd phiPrev;

		Eigen::Matrix<bool,Eigen::Dynamic, Eigen::Dynamic> fixed;

		Eigen::MatrixXd kappa;


		Eigen::MatrixXd prevU;
		Eigen::MatrixXd prevV;
		double dtprev;
		
		
		int nx;
		int ny;






		int NumberOfPointsX;
		int NumberOfPointsY;

		int NumberOfCellX;
		int NumberOfCellY;

		int NumberOfPhysicalCellX;
		int NumberOfPhysicalCellY;




		/************************************************/
		/*                                              */
		/*		Velocities along x(u) and y(v) axis		*/
		/*												*/
		/*		uprev, vprev : (n-1)th time level		*/
		/*		unow, vnow	 : (n)th   time level		*/
		/*		unext, vnext : (n+1)th time level		*/
		/*												*/
		/************************************************/
		Eigen::MatrixXd u;
		Eigen::MatrixXd v;

		Eigen::MatrixXd ustar;
		Eigen::MatrixXd vstar;


		/************************************************/
		/*                                              */
		/*		Grid points								*/
		/*												*/
		/************************************************/
		Eigen::VectorXd x;
		Eigen::VectorXd y;

		
		/************************************************/
		/*                                              */
		/*		Smoothed functions						*/
		/*												*/
		/************************************************/
		double Heaviside(double const Phi) const {

			if (Phi <= -alpha) return 0.0;
			if (Phi >= +alpha) return 1.0;

			return 0.5 * (1.0 + Phi / alpha + sin(Pi * Phi / alpha) / Pi);

		};
		double Dirac(double const Phi) const {

			if (abs(Phi) >= alpha)
				return 0.0;

			return (0.5 / alpha)* (1.0 + cos(Pi * Phi / alpha));
		
		};
		double Signum(double const Phi) const {

			//return Phi / sqrt(Phi*Phi + dx * dy);

			return Phi / sqrt(Phi*Phi + alpha*alpha);

		};
		double Signum(int const i, int const j, Eigen::MatrixXd const & Phi) const {


			/*double const dphidx = 0.5 * (phiStar(i + 1, j) - phiStar(i - 1, j)) / dx;
			double const dphidy = 0.5 * (phiStar(i, j + 1) - phiStar(i, j - 1)) / dy;

			double const norm2 = sqrt(dphidx * dphidx + dphidy * dphidy);

			return phiStar(i, j) / sqrt(phiStar(i, j)*phiStar(i, j) + norm2 * dx * dy);*/



			/*double const dphidx = 0.5 * (Phi(i + 1, j) - Phi(i - 1, j)) / dx;
			double const dphidy = 0.5 * (Phi(i, j + 1) - Phi(i, j - 1)) / dy;

			double const norm2 = sqrt(dphidx * dphidx + dphidy * dphidy);

			return Phi(i, j) / sqrt(Phi(i, j)*Phi(i, j) + norm2 * dx * dy);*/


			
			
			double const a = (Phi(i, j) - Phi(i - 1, j)) / dx;
			double const b = (Phi(i + 1, j) - Phi(i, j)) / dx;

			double const c = (Phi(i, j) - Phi(i, j - 1)) / dy;
			double const d = (Phi(i, j + 1) - Phi(i, j)) / dy;

			double norm2 = 0.0;

			if (Phi(i, j) > 0.0) {

				double const ap = std::max(a, 0.0);
				double const bm = std::min(b, 0.0);

				double const cp = std::max(c, 0.0);
				double const dm = std::min(d, 0.0);

				norm2 = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm));

			}
			else if (Phi(i, j) < 0.0) {

				double const am = std::min(a, 0.0);
				double const bp = std::max(b, 0.0);

				double const cm = std::min(c, 0.0);
				double const dp = std::max(d, 0.0);

				norm2 = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp));

			}

			return Phi(i, j) / sqrt(Phi(i, j)*Phi(i, j) + norm2 * dx * dy);
			

		};

		double mu(double const Phi) const {

			return mu2 + (mu1 - mu2) * Heaviside(Phi);

		};
		double rho(double const Phi) const {

			return rho2 + (rho1 - rho2) * Heaviside(Phi);

		};



		/*****************************************************************************/
		/*                                                                           */
		/*		MOMENTUM EQUATION OPERATORS											 */
		/*																			 */
		/*																			 */
		/*		UDotDivU_x()	Computes [U dot (Div U)] 1st row					 */
		/*																			 */
		/*		UDotDivU_y()	Computes [U dot (Div U)] 2nd row					 */
		/*																			 */
		/*		DivDotTau_x()	Computes [Div dot Tau] 1st row						 */
		/*																			 */
		/*		DivDotTau_y()	Computes [Div dot Tau] 2nd row						 */
		/*                                                                           */
		/*      SurfaceTensionForce()  Computes surface tension term                 */
		/*                                                                           */
		/*		GravityForce()	Computes gravity force								 */
		/*                                                                           */
		/*****************************************************************************/
		double UDotDivU_x(int const i, int const j) const {


			/************************************************/
			/*                                              */
			/*		Gradient form	1						*/
			/*												*/
			/************************************************/
			/*double const ududx1 = 0.5 * (u(i - 1, j) + u(i, j))*(u(i, j) - u(i - 1, j)) / dx;
			double const ududx2 = 0.5 * (u(i + 1, j) + u(i, j))*(u(i + 1, j) - u(i, j)) / dx;

			double const vdudy1 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1))*(u(i, j) - u(i, j - 1)) / dy;
			double const vdudy2 = 0.5 * (v(i, j) + v(i + 1, j))*(u(i, j + 1) - u(i, j)) / dy;

			return 0.5 * (ududx1 + ududx2) + 0.5*(vdudy1 + vdudy2);*/


			/************************************************/
			/*                                              */
			/*		Gradient form 2							*/
			/*												*/
			/************************************************/
			/*double const vmean = 0.25 * (v(i, j) + v(i, j - 1) + v(i + 1, j) + v(i + 1, j - 1));
			double const ududx = u(i, j) * (0.5 / dx) * (u(i + 1, j) - u(i - 1, j));
			double const vdudy = vmean * (0.5 / dy) * (u(i, j + 1) - u(i, j - 1));

			return ududx + vdudy;*/


			/************************************************/
			/*                                              */
			/*		Divergence form	1						*/
			/*												*/
			/************************************************/
			/*double const ur = 0.5 * (u(i + 1, j) + u(i, j));
			double const ul = 0.5 * (u(i, j) + u(i - 1, j));

			double const u01 = 0.5 * (u(i, j + 1) + u(i, j));
			double const v01 = 0.5 * (v(i + 1, j) + v(i, j));

			double const u11 = 0.5 * (u(i, j) + u(i, j - 1));
			double const v11 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));

			double const duudx = (ur * ur - ul * ul) / dx;
			double const duvdy = (u01 * v01 - u11 * v11) / dy;

			return duudx + duvdy;*/


			
			/************************************************/
			/*                                              */
			/*		Upwind Gradient form	1a				*/
			/*												*/
			/************************************************/
			double const vmean = 0.25 * (v(i, j) + v(i, j - 1) + v(i + 1, j) + v(i + 1, j - 1));
			double const dudxm = (u(i + 1, j) - u(i, j)) / dx;
			double const dudxp = (u(i, j) - u(i - 1, j)) / dx;
			double const dudym = (u(i, j + 1) - u(i, j)) / dy;
			double const dudyp = (u(i, j) - u(i, j - 1)) / dy;

			double const ududx = u(i, j) >= 0.0 ? u(i, j) * dudxp : u(i, j) * dudxm;
			double const vdudy = vmean   >= 0.0 ? vmean * dudyp : vmean * dudym;

			return vdudy + ududx;


			/************************************************/
			/*                                              */
			/*		Upwind Gradient form	1b				*/
			/*												*/
			/************************************************/
			/*double const vmean = 0.25 * (v(i, j) + v(i, j - 1) + v(i + 1, j) + v(i + 1, j - 1));
			double const dudxm = (u(i + 1, j) - u(i, j)) / dx;
			double const dudxp = (u(i, j) - u(i - 1, j)) / dx;

			double const ududx = u(i, j) >= 0.0 ? u(i, j) * dudxp : u(i, j) * dudxm;
			double const vdudy = vmean * 0.5 * (u(i, j + 1) - u(i, j - 1)) / dy;

			return vdudy + ududx;*/


			/************************************************/
			/*                                              */
			/*		Upwind Divergence form					*/
			/*												*/
			/************************************************/
			/*double const uc = 0.5 * (u(i, j) + u(i - 1, j));
			double const ur = 0.5 * (u(i, j) + u(i + 1, j));
			double const uu1 = 0.5 * (u(i, j) + u(i - 1, j));
			double const uu2 = 0.5 * (u(i, j) - u(i - 1, j));
			double const uu3 = 0.5 * (u(i, j) + u(i + 1, j));
			double const uu4 = 0.5 * (u(i + 1, j) - u(i, j));

			double const uuc = uc * uu1 - std::abs(uc) * uu2;
			double const uur = ur * uu3 - std::abs(ur) * uu4;
			double const duudx = (uur - uuc) / dx;

			
			double const v01 = 0.5 * (v(i, j) + v(i + 1, j));
			double const v11 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));

			double const vu1 = 0.5 * (u(i, j) + u(i, j - 1));
			double const vu2 = 0.5 * (u(i, j) - u(i, j - 1));
			double const vu3 = 0.5 * (u(i, j) + u(i, j + 1));
			double const vu4 = 0.5 * (u(i, j + 1) - u(i, j));

			double const vuup = v01 * vu3 - std::abs(v01) * vu4;
			double const vud  = v11 * vu1 - std::abs(v11) * vu2;

			double const dvudy = (vuup - vud) / dy;*/


			/*double const vc = 0.5 * (v(i, j) + v(i, j - 1));
			double const umeanu = 0.25 * (u(i, j) + u(i, j + 1) + u(i - 1, j + 1) + u(i - 1, j));
			double const umeand = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j - 1) + u(i, j - 1));
			double const vuc = vc * 0.5*(umeanu + umeand) - std::abs(vc) * 0.5 * (umeanu - umeand);

			double const vup = 0.5 * (v(i, j) + v(i, j + 1));
			double const umeanup = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j - 1) + u(i, j - 1));
			double const vuc = vc * 0.5*(umeanu + umeand) - std::abs(vc) * 0.5 * (umeanu - umeand);*/


			//double const uc = 0.5 * (u(i, j) + u(i - 1, j));
			//double const uuc = uc >= 0.0 ? uc * u(i - 1, j) : uc * u(i, j);
			//double const ur = 0.5 * (u(i, j) + u(i + 1, j));
			//double const uur = ur >= 0.0 ? ur * u(i, j) : ur * u(i + 1, j);

			//double const duudx = (uur - uuc) / dx;

			/*double const vn = 0.5 * (v(i, j) + v(i + 1, j));
			double const vs = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));
			double const vuc = vc >= 0.0 ? vc * u(i - 1, j) : uc * u(i, j);
			double const vu = 0.5 * (v(i, j) + v(i, j + 1));
			double const uur = ur >= 0.0 ? ur * u(i, j) : ur * u(i + 1, j);

			double const dvudy*/












			/*double const vmean = 0.25 * (v(i, j) + v(i, j - 1) + v(i + 1, j) + v(i + 1, j - 1));
			double const ududx = u(i, j) > 0.0 ? u(i, j) * (1.0 / dx) * (u(i, j) - u(i - 1, j)) : u(i, j) * (1.0 / dx) * (u(i + 1, j) - u(i, j));
			double const vdudy = vmean > 0.0 ? vmean * (1.0 / dy) * (u(i, j) - u(i, j - 1)) : vmean * (1.0 / dy) * (u(i, j + 1) - u(i, j));

			return ududx + vdudy;*/

	/*		double const umeanUP = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
			double const umeanDOWN = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));*/
			/*double const ududx = 0.5 * (u(i, j) + u(i - 1, j))*(u(i, j) - u(i - 1, j)) / dx;
			double const vdudy = 0.5 * (v(i, j) + v(i, j - 1))*(u(i, j) - u(i, j - 1)) / dy;

			return ududx + vdudy;*/

		};
		double UDotDivU_y(int const i, int const j) const {


			/************************************************/
			/*                                              */
			/*		Gradient form 1							*/
			/*												*/
			/************************************************/
			/*double const udvdx1 = 0.5 * (u(i - 1, j) + u(i - 1, j + 1))*(v(i, j) - v(i - 1, j)) / dx;
			double const udvdx2 = 0.5 * (u(i, j) + u(i, j + 1))*(v(i + 1, j) - v(i, j)) / dx;

			double const vdvdy1 = 0.5 * (v(i, j - 1) + v(i, j))*(v(i, j) - v(i, j - 1)) / dy;
			double const vdvdy2 = 0.5 * (v(i, j + 1) + v(i, j))*(v(i, j + 1) - v(i, j)) / dy;

			return 0.5 * (udvdx1 + udvdx2) + 0.5*(vdvdy1 + vdvdy2);*/


			/************************************************/
			/*                                              */
			/*		Gradient form 2							*/
			/*												*/
			/************************************************/
			/*double const umean = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
			double const udvdx = umean * (0.5 / dx) * (v(i + 1, j) - v(i - 1, j));
			double const vdvdy = v(i, j) * (0.5 / dy) * (v(i, j + 1) - v(i, j - 1));

			return udvdx + vdvdy;*/


			/************************************************/
			/*                                              */
			/*		Divergence form	1						*/
			/*												*/
			/************************************************/
			/*double const vu = 0.5 * (v(i, j + 1) + v(i, j));
			double const vd = 0.5 * (v(i, j) + v(i, j - 1));

			double const u00 = 0.5 * (u(i - 1, j + 1) + u(i - 1, j));
			double const v00 = 0.5 * (v(i - 1, j) + v(i, j));

			double const u01 = 0.5 * (u(i, j + 1) + u(i, j));
			double const v01 = 0.5 * (v(i + 1, j) + v(i, j));

			double const duvdx = (u01 * v01 - u00 * v00) / dx;
			double const dvvdy = (vu * vu - vd * vd) / dy;

			return duvdx + dvvdy;*/


			/************************************************/
			/*                                              */
			/*		Upwind Gradient form	1a				*/
			/*												*/
			/************************************************/
			double const umean = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
			double const dvdxm = (v(i + 1, j) - v(i, j)) / dx;
			double const dvdxp = (v(i, j) - v(i - 1, j)) / dx;
			double const dvdym = (v(i, j + 1) - v(i, j)) / dy;
			double const dvdyp = (v(i, j) - v(i, j - 1)) / dy;

			double const udvdx = umean   >= 0.0 ? umean * dvdxp : umean * dvdxm;
			double const vdvdy = v(i, j) >= 0.0 ? v(i, j) * dvdyp : v(i, j) * dvdym;

			return udvdx + vdvdy;


			/************************************************/
			/*                                              */
			/*		Upwind Gradient form	1b				*/
			/*												*/
			/************************************************/
			/*double const umean = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
			double const dvdym = (v(i, j + 1) - v(i, j)) / dy;
			double const dvdyp = (v(i, j) - v(i, j - 1)) / dy;

			double const udvdx = umean * 0.5 * (v(i + 1, j) - v(i - 1, j)) / dx;
			double const vdvdy = v(i, j) >= 0.0 ? v(i, j) * dvdyp : v(i, j) * dvdym;

			return udvdx + vdvdy;*/
			

			/************************************************/
			/*                                              */
			/*		Upwind Divergence form					*/
			/*												*/
			/************************************************/
			/*double const vc = 0.5 * (u(i, j) + u(i - 1, j));
			double const vn = 0.5 * (u(i, j) + u(i + 1, j));
			double const uu1 = 0.5 * (u(i, j) + u(i - 1, j));
			double const uu2 = 0.5 * (u(i, j) - u(i - 1, j));
			double const uu3 = 0.5 * (u(i, j) + u(i + 1, j));
			double const uu4 = 0.5 * (u(i + 1, j) - u(i, j));

			double const uuc = uc * uu1 - std::abs(uc) * uu2;
			double const uur = ur * uu3 - std::abs(ur) * uu4;
			double const duudx = (uur - uuc) / dx;


			double const v01 = 0.5 * (v(i, j) + v(i + 1, j));
			double const v11 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));

			double const vu1 = 0.5 * (u(i, j) + u(i, j - 1));
			double const vu2 = 0.5 * (u(i, j) - u(i, j - 1));
			double const vu3 = 0.5 * (u(i, j) + u(i, j + 1));
			double const vu4 = 0.5 * (u(i, j + 1) - u(i, j));

			double const vuup = v01 * vu3 - std::abs(v01) * vu4;
			double const vud = v11 * vu1 - std::abs(v11) * vu2;

			double const dvudy = (vuup - vud) / dy;*/






			/*double const vdvdy = v(i, j) > 0.0 ? v(i, j) * (1.0 / dy) * (v(i, j) - v(i, j - 1)) : v(i, j) * (1.0 / dy) * (v(i, j + 1) - v(i, j));
			double const umean = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
			double const udvdx = umean > 0.0 ? umean * (1.0 / dx) * (v(i, j) - v(i - 1, j)) : umean *(1.0 / dx) * (v(i + 1, j) - v(i, j));

			return udvdx + vdvdy;*/

			/*double const udvdx = 0.5 * (u(i, j) + u(i - 1, j))*(v(i, j) - v(i - 1, j)) / dx;
			double const vdvdy = 0.5 * (v(i, j) + v(i, j - 1))*(v(i, j) - v(i, j - 1)) / dy;

			return udvdx + vdvdy;*/

		};

		double DivDotTau_x(int const i, int const j) const {

			/*double const phi_01 = 0.25 * (phi(i, j) + phi(i + 1, j) + phi(i + 1, j + 1) + phi(i, j + 1));
			double const phi_11 = 0.25 * (phi(i, j) + phi(i, j - 1) + phi(i + 1, j - 1) + phi(i + 1, j));

			double const mu_r = mu(phi(i + 1, j));
			double const mu_c = mu(phi(i, j));

			double const mu_01 = mu(phi_01);
			double const mu_11 = mu(phi_11);

			double const dudx_r = (u(i + 1, j) - u(i, j)) / dx;
			double const dudx_c = (u(i, j) - u(i - 1, j)) / dx;

			double const dudy_01 = (u(i, j + 1) - u(i, j)) / dy;
			double const dvdx_01 = (v(i + 1, j) - v(i, j)) / dx;

			double const dudy_11 = (u(i, j) - u(i, j - 1)) / dy;
			double const dvdx_11 = (v(i + 1, j - 1) - v(i, j - 1)) / dx;


			double const divTau_00 = 2.0 * (mu_r * dudx_r - mu_c * dudx_c) / dx;
			double const divTau_01 = (mu_01 * (dudy_01 + dvdx_01) - mu_11 * (dudy_11 + dvdx_11)) / dy;


			return divTau_00 + divTau_01;*/




			double const phi_01 = 0.25 * (phi(i, j) + phi(i + 1, j) + phi(i + 1, j + 1) + phi(i, j + 1));
			double const phi_11 = 0.25 * (phi(i, j) + phi(i, j - 1) + phi(i + 1, j - 1) + phi(i + 1, j));

			double const mu_01 = mu(phi_01);
			double const mu_11 = mu(phi_11);

			/************************************************/
			/*                                              */
			/*		[ 2*mu*[u_x] ]_x						*/
			/*												*/
			/************************************************/
			double const mu_r = mu(phi(i + 1, j));
			double const mu_c = mu(phi(i, j));

			double const dudx_r = (u(i + 1, j) - u(i, j)) / dx;
			double const dudx_c = (u(i, j) - u(i - 1, j)) / dx;

			double const divTau_00 = 2.0 * (mu_r * dudx_r - mu_c * dudx_c) / dx;

			/************************************************/
			/*                                              */
			/*		[ mu*[u_y] ]_y							*/
			/*												*/
			/************************************************/
			double const dudy_01_a = (u(i, j + 1) - u(i, j)) / dy;
			double const dudy_11_a = (u(i, j) - u(i, j - 1)) / dy;

			double const divTau_01_a = ( mu_01 * dudy_01_a  - mu_11 * dudy_11_a) / dy;
			
			/************************************************/
			/*                                              */
			/*		[ mu*[v_x] ]_y							*/
			/*												*/
			/************************************************/
			double const dvdx_01_b = (v(i + 1, j) - v(i, j)) / dx;
			double const dvdx_11_b = (v(i + 1, j - 1) - v(i, j - 1)) / dx;

			double const divTau_01_b = (mu_01 * dvdx_01_b - mu_11 * dvdx_11_b) / dy;



			return divTau_00 + (divTau_01_a + divTau_01_b);

		};
		double DivDotTau_y(int const i, int const j) const {

			/*double const phi_00 = 0.25 * (phi(i, j) + phi(i, j + 1) + phi(i - 1, j + 1) + phi(i - 1, j));
			double const phi_01 = 0.25 * (phi(i, j) + phi(i + 1, j) + phi(i + 1, j + 1) + phi(i, j + 1));

			double const mu_c = mu(phi(i, j));
			double const mu_u = mu(phi(i, j + 1));

			double const mu_00 = mu(phi_00);
			double const mu_01 = mu(phi_01);

			double const dvdy_u = (v(i, j + 1) - v(i, j)) / dy;
			double const dvdy_c = (v(i, j) - v(i, j - 1)) / dy;

			double const dudy_00 = (u(i - 1, j + 1) - u(i - 1, j)) / dy;
			double const dvdx_00 = (v(i, j) - v(i - 1, j)) / dx;

			double const dudy_01 = (u(i, j + 1) - u(i, j)) / dy;
			double const dvdx_01 = (v(i + 1, j) - v(i, j)) / dx;


			double const divTau_10 = (mu_01 * (dudy_01 + dvdx_01) - mu_00 * (dudy_00 + dvdx_00)) / dx;
			double const divTau_11 = 2.0 * (mu_u * dvdy_u - mu_c * dvdy_c) / dy;


			return divTau_10 + divTau_11;*/


			double const phi_00 = 0.25 * (phi(i, j) + phi(i, j + 1) + phi(i - 1, j + 1) + phi(i - 1, j));
			double const phi_01 = 0.25 * (phi(i, j) + phi(i + 1, j) + phi(i + 1, j + 1) + phi(i, j + 1));

			double const mu_00 = mu(phi_00);
			double const mu_01 = mu(phi_01);

			/************************************************/
			/*                                              */
			/*		[ 2*mu*[v_y] ]_y						*/
			/*												*/
			/************************************************/
			double const mu_c = mu(phi(i, j));
			double const mu_u = mu(phi(i, j + 1));

			double const dvdy_u = (v(i, j + 1) - v(i, j)) / dy;
			double const dvdy_c = (v(i, j) - v(i, j - 1)) / dy;

			double const divTau_11 = 2.0 * (mu_u * dvdy_u - mu_c * dvdy_c) / dy;

			/************************************************/
			/*                                              */
			/*		[ mu*[u_y] ]_x							*/
			/*												*/
			/************************************************/
			double const dudy_01_a = (u(i, j + 1) - u(i, j)) / dy;
			double const dudy_00_a = (u(i - 1, j + 1) - u(i - 1, j)) / dy;

			double const divTau_10_a = (mu_01 * dudy_01_a - mu_00 * dudy_00_a) / dx;

			/************************************************/
			/*                                              */
			/*		[ mu*[v_x] ]_x							*/
			/*												*/
			/************************************************/
			double const dvdx_01_b = (v(i + 1, j) - v(i, j)) / dx;
			double const dvdx_00_b = (v(i, j) - v(i - 1, j)) / dx;

			double const divTau_10_b = (mu_01 * dvdx_01_b - mu_00 * dvdx_00_b) / dx;



			return (divTau_10_a + divTau_10_b) + divTau_11;


		};
		
		double SurfaceTensionForce_x(int const i, int const j) {


			double const phi_r = 0.5 * (phi(i, j) + phi(i + 1, j));

			// Because Dirac delta function is zero
			if (fabs(phi_r) > alpha)
				return 0.0;

			double const norm1 = NormGradPhi(i, j);
			double const norm2 = NormGradPhi(i + 1, j);

			//double const curvature = 0.5 * (DivGradPhi(i, j) + DivGradPhi(i + 1, j)) * (phi(i + 1, j) - phi(i, j)) / dx;
			double const curvature = 0.5 * (DivGradPhi(i, j) * DPhiDx(i, j) / norm1 + DivGradPhi(i + 1, j) * DPhiDx(i + 1, j) / norm2);

			return sigma * Dirac(phi_r) * curvature;


		};
		double SurfaceTensionForce_y(int const i, int const j) const {


			double const phi_u = 0.5 * (phi(i, j) + phi(i, j + 1));

			// Because Dirac delta function is zero
			if (fabs(phi_u) > alpha)
				return 0.0;

			double const norm1 = NormGradPhi(i, j);
			double const norm2 = NormGradPhi(i, j + 1);

			//double const curvature = 0.5 * (DivGradPhi(i, j) + DivGradPhi(i, j + 1))  * (phi(i, j + 1) - phi(i, j)) / dy;
			double const curvature = 0.5 * (DivGradPhi(i, j) * DPhiDy(i, j) / norm1 + DivGradPhi(i, j + 1) * DPhiDy(i, j + 1) / norm2);

			return sigma * Dirac(phi_u) * curvature;

		};


		double DPhiDx(int const i, int const j) const {

			if (i == 0)
				return (phi(i + 1, j) - phi(i, j)) / dx;
			if (i == nx + 1)
				return (phi(i, j) - phi(i - 1, j)) / dx;

			return 0.5 * (phi(i + 1, j) - phi(i - 1, j)) / dx;

		};
		double DPhiDy(int const i, int const j) const {

			if (j == 0)
				return (phi(i, j + 1) - phi(i, j)) / dy;
			if (j == ny + 1)
				return (phi(i, j) - phi(i, j - 1)) / dy;

			return 0.5 * (phi(i, j + 1) - phi(i, j - 1)) / dy;

		};
		double NormGradPhi(int const i, int const j) const {

			double dphidx = DPhiDx(i, j);
			double dphidy = DPhiDy(i, j);

			return sqrt(dphidx*dphidx + dphidy*dphidy);

		};

		double DivGradPhi(int const i, int const j) const {


			double val_x;
			double val_y;

			/************************************************/
			/*                                              */
			/*		x - part of the divergence				*/
			/*												*/
			/************************************************/
			if (i == 0) {

				double const dphidx_r = DPhiDx(i + 1, j);
				double const dphidy_r = DPhiDy(i + 1, j);
				double const norm_r = sqrt(dphidx_r*dphidx_r + dphidy_r * dphidy_r);

				double const dphidx_l = DPhiDx(i, j);
				double const dphidy_l = DPhiDy(i, j);
				double const norm_l = sqrt(dphidx_l*dphidx_l + dphidy_l * dphidy_l);

				val_x = (1.0 / dx) * (dphidx_r / norm_r - dphidx_l / norm_l);

			}
			else if (i == nx + 1) {

				double const dphidx_r = DPhiDx(i, j);
				double const dphidy_r = DPhiDy(i, j);
				double const norm_r = sqrt(dphidx_r*dphidx_r + dphidy_r * dphidy_r);

				double const dphidx_l = DPhiDx(i - 1, j);
				double const dphidy_l = DPhiDy(i - 1, j);
				double const norm_l = sqrt(dphidx_l*dphidx_l + dphidy_l * dphidy_l);

				val_x = (1.0 / dx) * (dphidx_r / norm_r - dphidx_l / norm_l);

			}
			else {

				double const dphidx_r = DPhiDx(i + 1, j);
				double const dphidy_r = DPhiDy(i + 1, j);
				double const norm_r   = sqrt(dphidx_r*dphidx_r + dphidy_r * dphidy_r);

				double const dphidx_l = DPhiDx(i - 1, j);
				double const dphidy_l = DPhiDy(i - 1, j);
				double const norm_l   = sqrt(dphidx_l*dphidx_l + dphidy_l * dphidy_l);

				val_x = (0.5 / dx) * (dphidx_r / norm_r - dphidx_l / norm_l);

			}

			/************************************************/
			/*                                              */
			/*		y - part of the divergence				*/
			/*												*/
			/************************************************/
			if (j == 0) {

				double const dphidx_u = DPhiDx(i, j + 1);
				double const dphidy_u = DPhiDy(i, j + 1);
				double const norm_u = sqrt(dphidx_u*dphidx_u + dphidy_u * dphidy_u);

				double const dphidx_d = DPhiDx(i, j);
				double const dphidy_d = DPhiDy(i, j);
				double const norm_d = sqrt(dphidx_d*dphidx_d + dphidy_d * dphidy_d);

				val_y = (1.0 / dy) * (dphidy_u / norm_u - dphidy_d / norm_d);

			}
			else if (j == ny + 1) {

				double const dphidx_u = DPhiDx(i, j);
				double const dphidy_u = DPhiDy(i, j);
				double const norm_u = sqrt(dphidx_u*dphidx_u + dphidy_u * dphidy_u);

				double const dphidx_d = DPhiDx(i, j - 1);
				double const dphidy_d = DPhiDy(i, j - 1);
				double const norm_d = sqrt(dphidx_d*dphidx_d + dphidy_d * dphidy_d);

				val_y = (1.0 / dy) * (dphidy_u / norm_u - dphidy_d / norm_d);

			}
			else {

				double const dphidx_u = DPhiDx(i, j + 1);
				double const dphidy_u = DPhiDy(i, j + 1);
				double const norm_u   = sqrt(dphidx_u*dphidx_u + dphidy_u * dphidy_u);

				double const dphidx_d = DPhiDx(i, j - 1);
				double const dphidy_d = DPhiDy(i, j - 1);
				double const norm_d   = sqrt(dphidx_d*dphidx_d + dphidy_d * dphidy_d);

				val_y = (0.5 / dy) * (dphidy_u / norm_u - dphidy_d / norm_d);

			}

			return val_x + val_y;






			/*double const dphidx			= DPhiDx(i, j);
			double const dphidy			= DPhiDy(i, j);

			double const ddphidxx		= (phi(i + 1, j) - 2.0 * phi(i, j) + phi(i - 1, j)) / (dx*dx);
			double const ddphidyy		= (phi(i, j + 1) - 2.0 * phi(i, j) + phi(i, j - 1)) / (dy*dy);

			double const ddphidxy		= 0.25 * ((phi(i + 1, j + 1) - phi(i - 1, j + 1)) - (phi(i + 1, j - 1) - phi(i - 1, j - 1))) / (dx*dy);

			double const normGradPhi	= NormGradPhi(i, j);

			return (dphidx*dphidx * ddphidyy - 2.0*dphidx*dphidy*ddphidxy + dphidy * dphidy*ddphidxx) / (normGradPhi*normGradPhi*normGradPhi);*/


		};


		/*****************************************************************************/
		/*                                                                           */
		/*		LEVEL SET OPERATORS													 */
		/*																			 */
		/*																			 */
		/*		UDotDivPhi()	Computes [U dot (Div Phi)]							 */
		/*																			 */
		/*		UpwindPhi_x()	Upwind [Div Phi] 1st row							 */
		/*																			 */
		/*		UpwindPhi_y()	Upwind [Div Phi] 2nd row							 */
		/*                                                                           */
		/*                                                                           */
		/*****************************************************************************/
		/*double UDotDivPhi(int const i, int const j) {



			double uij;
			double vij;

			//uij = 0.5 * (u(i, j) + u(i - 1, j));
			//vij = 0.5 * (v(i, j) + v(i, j - 1));


			//uij = 0.5 * (ustar(i, j) + ustar(i - 1, j));
			//vij = 0.5 * (vstar(i, j) + vstar(i, j - 1));



			if (i == 1)
				uij = (5.0 * u(i - 1, j) + 15.0* u(i, j) - 5.0*u(i + 1, j) + u(i + 2, j)) / 16.0;
			else if (i == nx)
				uij = (u(i - 3, j) - 5.0* u(i - 2, j) + 15.0*u(i - 1, j) + 5.0*u(i, j)) / 16.0;
			else
				uij = (-u(i - 2, j) + 9.0* u(i - 1, j) + 9.0*u(i, j) - u(i + 1, j)) / 16.0;


			if (j == 1)
				vij = (5.0 * v(i, j - 1) + 15.0* v(i, j) - 5.0*v(i, j + 1) + v(i, j + 2)) / 16.0;
			else if (j == ny)
				vij = (v(i, j - 3) - 5.0* v(i, j - 2) + 15.0*v(i, j - 1) + 5.0*v(i, j)) / 16.0;
			else
				vij = (-v(i, j - 2) + 9.0* v(i, j - 1) + 9.0*v(i, j) - v(i, j + 1)) / 16.0;

			//if (i == 1)
			//	uij = (5.0 * ustar(i - 1, j) + 15.0* ustar(i, j) - 5.0*ustar(i + 1, j) + ustar(i + 2, j)) / 16.0;
			//else if (i == nx)
			//	uij = (ustar(i - 3, j) - 5.0* ustar(i - 2, j) + 15.0*ustar(i - 1, j) + 5.0*ustar(i, j)) / 16.0;
			//else
			//	uij = (-ustar(i - 2, j) + 9.0* ustar(i - 1, j) + 9.0*ustar(i, j) - ustar(i + 1, j)) / 16.0;
			//
			//
			//if (j == 1)
			//	vij = (5.0 * vstar(i, j - 1) + 15.0* vstar(i, j) - 5.0*vstar(i, j + 1) + vstar(i, j + 2)) / 16.0;
			//else if (j == ny)
			//	vij = (vstar(i, j - 3) - 5.0* vstar(i, j - 2) + 15.0*vstar(i, j - 1) + 5.0*vstar(i, j)) / 16.0;
			//else
			//	vij = (-vstar(i, j - 2) + 9.0* vstar(i, j - 1) + 9.0*vstar(i, j) - vstar(i, j + 1)) / 16.0;


			double const dphidx = UpwindPhi_x(uij, i, j);
			double const dphidy = UpwindPhi_y(vij, i, j);

			return uij * dphidx + vij * dphidy;

		};

		double UpwindPhi_x(double const uij, int const i, int const j) {

			if (uij > 0.0)
				return (phi(i, j) - phi(i - 1, j)) / dx;
			if (uij < 0.0)
				return (phi(i + 1, j) - phi(i, j)) / dx;

			return 0.0;
			
		};
		double UpwindPhi_y(double const vij, int const i, int const j) {

			if (vij > 0.0)
				return (phi(i, j) - phi(i, j - 1)) / dy;
			if (vij < 0.0)
				return (phi(i, j + 1) - phi(i, j)) / dy;

			return 0.0;

		};
		*/

		double superbeeLimiter(double const x, double const y) {


			double const absx = std::abs(x);
			double const absy = std::abs(y);

			double const signx = x > 0.0 ? 1.0 : -1.0;

			if (0.5*absx <= absy && absy <= 2.0*absx && x*y > 0.0)
				return signx * std::max(absx, absy);
			if ((0.5*absx >= absy || absy >= 2.0 * absx) && x*y > 0.0)
				return 2.0*signx*std::min(absx, absy);

			return 0.0;

		};

		double F(int const i, int const j, Eigen::MatrixXd const & Phi) {


			double const dphidx_forward_c  = (Phi(i + 1, j) - Phi(i, j)) / dx;
			double const dphidx_backward_c = (Phi(i, j) - Phi(i - 1, j)) / dx;

			double const dphidx_forward_r  = (Phi(i + 2, j) - Phi(i + 1, j)) / dx;
			double const dphidx_backward_r = (Phi(i + 1, j) - Phi(i, j)) / dx;


			double const sx_c = superbeeLimiter(dphidx_forward_c, dphidx_backward_c);
			double const sx_r = superbeeLimiter(dphidx_forward_r, dphidx_backward_r);


			double const phiMinus_x = Phi(i, j) + 0.5 * dx * sx_c;
			double const phiPlus_x  = Phi(i + 1, j) - 0.5 * dx * sx_r;


			return std::max(u(i, j), 0.0) * phiMinus_x + std::min(u(i, j), 0.0) * phiPlus_x;

		};
		double G(int const i, int const j, Eigen::MatrixXd const & Phi) {


			double const dphidy_forward_c  = (Phi(i, j + 1) - Phi(i, j)) / dy;
			double const dphidy_backward_c = (Phi(i, j) - Phi(i, j - 1)) / dy;

			double const dphidy_forward_u  = (Phi(i, j + 2) - Phi(i, j + 1)) / dy;
			double const dphidy_backward_u = (Phi(i, j + 1) - Phi(i, j)) / dy;


			double const sy_c = superbeeLimiter(dphidy_forward_c, dphidy_backward_c);
			double const sy_u = superbeeLimiter(dphidy_forward_u, dphidy_backward_u);


			double const phiMinus_y = Phi(i, j) + 0.5*dy *sy_c;
			double const phiPlus_y  = Phi(i, j + 1) - 0.5*dy *sy_u;


			return std::max(v(i, j), 0.0) * phiMinus_y + std::min(v(i, j), 0.0) * phiPlus_y;

		};

		Eigen::MatrixXd OperatorF(Eigen::MatrixXd const & Phi) {


			Eigen::MatrixXd result(Phi.rows(), Phi.cols());

			result.setZero();

			for (int i = 2; i < nx; i++) {
				for (int j = 2; j < ny; j++) {

					double const Value1 = -(F(i, j, Phi) - F(i - 1, j, Phi)) / dx;
					double const Value2 = -(G(i, j, Phi) - G(i, j - 1, Phi)) / dy;

					result.coeffRef(i, j) = Value1 + Value2;

				}
			}

			return result;

		};


		double UDotGradPhi(int const i, int const j) {



			double uij;
			double vij;

			//uij = 0.5 * (u(i, j) + u(i - 1, j));
			//vij = 0.5 * (v(i, j) + v(i, j - 1));

			if (i == 1)
				uij = (5.0 * u(i - 1, j) + 15.0* u(i, j) - 5.0*u(i + 1, j) + u(i + 2, j)) / 16.0;
			else if (i == nx)
				uij = (u(i - 3, j) - 5.0* u(i - 2, j) + 15.0*u(i - 1, j) + 5.0*u(i, j)) / 16.0;
			else
				uij = (-u(i - 2, j) + 9.0* u(i - 1, j) + 9.0*u(i, j) - u(i + 1, j)) / 16.0;


			if (j == 1)
				vij = (5.0 * v(i, j - 1) + 15.0* v(i, j) - 5.0*v(i, j + 1) + v(i, j + 2)) / 16.0;
			else if (j == ny)
				vij = (v(i, j - 3) - 5.0* v(i, j - 2) + 15.0*v(i, j - 1) + 5.0*v(i, j)) / 16.0;
			else
				vij = (-v(i, j - 2) + 9.0* v(i, j - 1) + 9.0*v(i, j) - v(i, j + 1)) / 16.0;



			double const dphidx = UpwindPhi_x(uij, i, j);
			double const dphidy = UpwindPhi_y(vij, i, j);

			return uij * dphidx + vij * dphidy;

		};

		double UpwindPhi_x(double const uij, int const i, int const j) {


			/*if (uij > 0.0)
				return (phi(i, j) - phi(i - 1, j)) / dx;
			if (uij < 0.0)
				return (phi(i + 1, j) - phi(i, j)) / dx;
			
			return 0.0;*/

			/**/
			if (uij > 0.0) {

				if (i == 1 || i == 2 || i == nx || i == nx - 1)
					return (phi(i, j) - phi(i - 1, j)) / dx;

				double const Value1 = omega1minus_x(i, j)*(q1minus_x(i, j) / 3.0 - 7.0 * q2minus_x(i, j) / 6.0 + 11.0 * q3minus_x(i, j) / 6.0);
				double const Value2 = omega2minus_x(i, j)*(-q2minus_x(i, j) / 6.0 + 5.0 * q3minus_x(i, j) / 6.0 + q4minus_x(i, j) / 3.0);
				double const Value3 = omega3minus_x(i, j)*(q3minus_x(i, j) / 3.0 + 5.0 * q4minus_x(i, j) / 6.0 - q5minus_x(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}
			else if (uij < 0.0) {

				if (i == 1 || i == 2 || i == nx - 1 || i == nx)
					return (phi(i + 1, j) - phi(i, j)) / dx;

				double const Value1 = omega1plus_x(i, j)*(q1plus_x(i, j) / 3.0 - 7.0 * q2plus_x(i, j) / 6.0 + 11.0 * q3plus_x(i, j) / 6.0);
				double const Value2 = omega2plus_x(i, j)*(-q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0);
				double const Value3 = omega3plus_x(i, j)*(q3plus_x(i, j) / 3.0 + 5.0 * q4plus_x(i, j) / 6.0 - q5plus_x(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}

			return 0.0;
			/**/


			if (uij > 0.0) {

				if (i == 1)
					return q3minus_x(i, j) / 3.0 + 5.0 * q4minus_x(i, j) / 6.0 - q5minus_x(i, j) / 6.0;
				else if (i == 2 || i == nx - 1)
					return -q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0;
				else if (i == nx)
					return q1minus_x(i, j) / 3.0 - 7.0 * q2minus_x(i, j) / 6.0 + 11.0 * q3minus_x(i, j) / 6.0;


				double const Value1 = omega1minus_x(i, j)*(q1minus_x(i, j) / 3.0 - 7.0 * q2minus_x(i, j) / 6.0 + 11.0 * q3minus_x(i, j) / 6.0);
				double const Value2 = omega2minus_x(i, j)*(-q2minus_x(i, j) / 6.0 + 5.0 * q3minus_x(i, j) / 6.0 + q4minus_x(i, j) / 3.0);
				double const Value3 = omega3minus_x(i, j)*(q3minus_x(i, j) / 3.0 + 5.0 * q4minus_x(i, j) / 6.0 - q5minus_x(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}
			else if (uij < 0.0) {

				if (i == 1)
					return q1plus_x(i, j) / 3.0 - 7.0 * q2plus_x(i, j) / 6.0 + 11.0 * q3plus_x(i, j) / 6.0;
				else if (i == 2 || i == nx - 1)
					return -q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0;
				else if (i == nx)
					return q3plus_x(i, j) / 3.0 + 5.0 * q4plus_x(i, j) / 6.0 - q5plus_x(i, j) / 6.0;

				double const Value1 = omega1plus_x(i, j)*(q1plus_x(i, j) / 3.0 - 7.0 * q2plus_x(i, j) / 6.0 + 11.0 * q3plus_x(i, j) / 6.0);
				double const Value2 = omega2plus_x(i, j)*(-q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0);
				double const Value3 = omega3plus_x(i, j)*(q3plus_x(i, j) / 3.0 + 5.0 * q4plus_x(i, j) / 6.0 - q5plus_x(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}

			return 0.0;

		};
		double UpwindPhi_y(double const vij, int const i, int const j) {

			/*if (vij > 0.0)
				return (phi(i, j) - phi(i, j - 1)) / dy;
			if (vij < 0.0)
				return (phi(i, j + 1) - phi(i, j)) / dy;
			
			return 0.0;*/

			/**/
			if (vij > 0.0) {

				if (j == 1 || j == 2 || j == ny - 1 || j == ny)
					return (phi(i, j) - phi(i, j - 1)) / dy;

				double const Value1 = omega1minus_y(i, j)*(q1minus_y(i, j) / 3.0 - 7.0 * q2minus_y(i, j) / 6.0 + 11.0 * q3minus_y(i, j) / 6.0);
				double const Value2 = omega2minus_y(i, j)*(-q2minus_y(i, j) / 6.0 + 5.0 * q3minus_y(i, j) / 6.0 + q4minus_y(i, j) / 3.0);
				double const Value3 = omega3minus_y(i, j)*(q3minus_y(i, j) / 3.0 + 5.0 * q4minus_y(i, j) / 6.0 - q5minus_y(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}
			else if (vij < 0.0) {

				if (j == 1 || j == 2 || j == ny - 1 || j == ny)
					return (phi(i, j + 1) - phi(i, j)) / dy;

				double const Value1 = omega1plus_y(i, j)*(q1plus_y(i, j) / 3.0 - 7.0 * q2plus_y(i, j) / 6.0 + 11.0 * q3plus_y(i, j) / 6.0);
				double const Value2 = omega2plus_y(i, j)*(-q2plus_y(i, j) / 6.0 + 5.0 * q3plus_y(i, j) / 6.0 + q4plus_y(i, j) / 3.0);
				double const Value3 = omega3plus_y(i, j)*(q3plus_y(i, j) / 3.0 + 5.0 * q4plus_y(i, j) / 6.0 - q5plus_y(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}

			return 0.0;
			/**/


			if (vij > 0.0) {

				if (j == 1)
					return q3minus_y(i, j) / 3.0 + 5.0 * q4minus_y(i, j) / 6.0 - q5minus_y(i, j) / 6.0;
				else if (j == 2 || j == ny - 1)
					return -q2minus_y(i, j) / 6.0 + 5.0 * q3minus_y(i, j) / 6.0 + q4minus_y(i, j) / 3.0;
				else if (j == ny)
					return q1minus_y(i, j) / 3.0 - 7.0 * q2minus_y(i, j) / 6.0 + 11.0 * q3minus_y(i, j) / 6.0;

				double const Value1 = omega1minus_y(i, j)*(q1minus_y(i, j) / 3.0 - 7.0 * q2minus_y(i, j) / 6.0 + 11.0 * q3minus_y(i, j) / 6.0);
				double const Value2 = omega2minus_y(i, j)*(-q2minus_y(i, j) / 6.0 + 5.0 * q3minus_y(i, j) / 6.0 + q4minus_y(i, j) / 3.0);
				double const Value3 = omega3minus_y(i, j)*(q3minus_y(i, j) / 3.0 + 5.0 * q4minus_y(i, j) / 6.0 - q5minus_y(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}
			else if (vij < 0.0) {

				if (j == 1)
					return q1plus_y(i, j) / 3.0 - 7.0 * q2plus_y(i, j) / 6.0 + 11.0 * q3plus_y(i, j) / 6.0;
				else if (j == 2 || j == ny - 1)
					return -q2plus_y(i, j) / 6.0 + 5.0 * q3plus_y(i, j) / 6.0 + q4plus_y(i, j) / 3.0;
				else if (j == ny)
					return q3plus_y(i, j) / 3.0 + 5.0 * q4plus_y(i, j) / 6.0 - q5plus_y(i, j) / 6.0;

				double const Value1 = omega1plus_y(i, j)*(q1plus_y(i, j) / 3.0 - 7.0 * q2plus_y(i, j) / 6.0 + 11.0 * q3plus_y(i, j) / 6.0);
				double const Value2 = omega2plus_y(i, j)*(-q2plus_y(i, j) / 6.0 + 5.0 * q3plus_y(i, j) / 6.0 + q4plus_y(i, j) / 3.0);
				double const Value3 = omega3plus_y(i, j)*(q3plus_y(i, j) / 3.0 + 5.0 * q4plus_y(i, j) / 6.0 - q5plus_y(i, j) / 6.0);

				return Value1 + Value2 + Value3;

			}

			return 0.0;

		};

		double q1plus_x(int const i, int const j) {

			return (phi(i + 3, j) - phi(i + 2, j)) / dx;

		};
		double q2plus_x(int const i, int const j) {

			return (phi(i + 2, j) - phi(i + 1, j)) / dx;

		};
		double q3plus_x(int const i, int const j) {

			return (phi(i + 1, j) - phi(i, j)) / dx;

		};
		double q4plus_x(int const i, int const j) {

			return (phi(i, j) - phi(i - 1, j)) / dx;

		};
		double q5plus_x(int const i, int const j) {

			return (phi(i - 1, j) - phi(i - 2, j)) / dx;

		};

		double q1minus_x(int const i, int const j) {

			return (phi(i - 2, j) - phi(i - 3, j)) / dx;

		};
		double q2minus_x(int const i, int const j) {

			return (phi(i - 1, j) - phi(i - 2, j)) / dx;

		};
		double q3minus_x(int const i, int const j) {

			return (phi(i, j) - phi(i - 1, j)) / dx;

		};
		double q4minus_x(int const i, int const j) {

			return (phi(i + 1, j) - phi(i, j)) / dx;

		};
		double q5minus_x(int const i, int const j) {

			return (phi(i + 2, j) - phi(i + 1, j)) / dx;

		};


		double q1plus_y(int const i, int const j) {

			return (phi(i, j + 3) - phi(i, j + 2)) / dy;

		};
		double q2plus_y(int const i, int const j) {

			return (phi(i, j + 2) - phi(i, j + 1)) / dy;

		};
		double q3plus_y(int const i, int const j) {

			return (phi(i, j + 1) - phi(i, j)) / dy;

		};
		double q4plus_y(int const i, int const j) {

			return (phi(i, j) - phi(i, j - 1)) / dy;

		};
		double q5plus_y(int const i, int const j) {

			return (phi(i, j - 1) - phi(i, j - 2)) / dy;

		};

		double q1minus_y(int const i, int const j) {

			return (phi(i, j - 2) - phi(i, j - 3)) / dy;

		};
		double q2minus_y(int const i, int const j) {

			return (phi(i, j - 1) - phi(i, j - 2)) / dy;

		};
		double q3minus_y(int const i, int const j) {

			return (phi(i, j) - phi(i, j - 1)) / dy;

		};
		double q4minus_y(int const i, int const j) {

			return (phi(i, j + 1) - phi(i, j)) / dy;

		};
		double q5minus_y(int const i, int const j) {

			return (phi(i, j + 2) - phi(i, j + 1)) / dy;

		};



		double IS1minus_x(int const i, int const j) {

			double const Value1 = q1minus_x(i, j) - 2.0*q2minus_x(i, j) + q3minus_x(i, j);
			double const Value2 = q1minus_x(i, j) - 4.0*q2minus_x(i, j) + 3.0*q3minus_x(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS2minus_x(int const i, int const j) {

			double const Value1 = q2minus_x(i, j) - 2.0*q3minus_x(i, j) + q4minus_x(i, j);
			double const Value2 = q2minus_x(i, j) - q4minus_x(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS3minus_x(int const i, int const j) {

			double const Value1 = q3minus_x(i, j) - 2.0*q4minus_x(i, j) + q5minus_x(i, j);
			double const Value2 = 3.0*q3minus_x(i, j) - 4.0*q4minus_x(i, j) + q5minus_x(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};

		double IS1plus_x(int const i, int const j) {

			double const Value1 = q1plus_x(i, j) - 2.0*q2plus_x(i, j) + q3plus_x(i, j);
			double const Value2 = q1plus_x(i, j) - 4.0*q2plus_x(i, j) + 3.0*q3plus_x(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS2plus_x(int const i, int const j) {

			double const Value1 = q2plus_x(i, j) - 2.0*q3plus_x(i, j) + q4plus_x(i, j);
			double const Value2 = q2plus_x(i, j) - q4plus_x(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS3plus_x(int const i, int const j) {

			double const Value1 = q3plus_x(i, j) - 2.0*q4plus_x(i, j) + q5plus_x(i, j);
			double const Value2 = 3.0*q3plus_x(i, j) - 4.0*q4plus_x(i, j) + q5plus_x(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};


		double IS1minus_y(int const i, int const j) {

			double const Value1 = q1minus_y(i, j) - 2.0*q2minus_y(i, j) + q3minus_y(i, j);
			double const Value2 = q1minus_y(i, j) - 4.0*q2minus_y(i, j) + 3.0*q3minus_y(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS2minus_y(int const i, int const j) {

			double const Value1 = q2minus_y(i, j) - 2.0*q3minus_y(i, j) + q4minus_y(i, j);
			double const Value2 = q2minus_y(i, j) - q4minus_y(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS3minus_y(int const i, int const j) {

			double const Value1 = q3minus_y(i, j) - 2.0*q4minus_y(i, j) + q5minus_y(i, j);
			double const Value2 = 3.0*q3minus_y(i, j) - 4.0*q4minus_y(i, j) + q5minus_y(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};

		double IS1plus_y(int const i, int const j) {

			double const Value1 = q1plus_y(i, j) - 2.0*q2plus_y(i, j) + q3plus_y(i, j);
			double const Value2 = q1plus_y(i, j) - 4.0*q2plus_y(i, j) + 3.0*q3plus_y(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS2plus_y(int const i, int const j) {

			double const Value1 = q2plus_y(i, j) - 2.0*q3plus_y(i, j) + q4plus_y(i, j);
			double const Value2 = q2plus_y(i, j) - q4plus_y(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};
		double IS3plus_y(int const i, int const j) {

			double const Value1 = q3plus_y(i, j) - 2.0*q4plus_y(i, j) + q5plus_y(i, j);
			double const Value2 = 3.0*q3plus_y(i, j) - 4.0*q4plus_y(i, j) + q5plus_y(i, j);

			return (13.0 / 12.0) * Value1* Value1 + 0.25* Value2*Value2;

		};



		double alpha1plus_x(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS1plus_x(i, j);

			return (1.0 / 10.0) / (Value*Value);

		};
		double alpha2plus_x(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS2plus_x(i, j);

			return (6.0 / 10.0) / (Value*Value);

		};
		double alpha3plus_x(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS3plus_x(i, j);

			return (3.0 / 10.0) / (Value*Value);

		};

		double alpha1minus_x(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS1minus_x(i, j);

			return (1.0 / 10.0) / (Value*Value);

		};
		double alpha2minus_x(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS2minus_x(i, j);

			return (6.0 / 10.0) / (Value*Value);

		};
		double alpha3minus_x(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS3minus_x(i, j);

			return (3.0 / 10.0) / (Value*Value);

		};


		double alpha1plus_y(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS1plus_y(i, j);

			return (1.0 / 10.0) / (Value*Value);

		};
		double alpha2plus_y(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS2plus_y(i, j);

			return (6.0 / 10.0) / (Value*Value);

		};
		double alpha3plus_y(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS3plus_y(i, j);

			return (3.0 / 10.0) / (Value*Value);

		};

		double alpha1minus_y(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS1minus_y(i, j);

			return (1.0 / 10.0) / (Value*Value);

		};
		double alpha2minus_y(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS2minus_y(i, j);

			return (6.0 / 10.0) / (Value*Value);

		};
		double alpha3minus_y(int const i, int const j) {

			double const eps = 1e-6;

			double const Value = eps + IS3minus_y(i, j);

			return (3.0 / 10.0) / (Value*Value);

		};



		double omega1plus_x(int const i, int const j) {

			double const Value = alpha1plus_x(i, j) + alpha2plus_x(i, j) + alpha3plus_x(i, j);

			return alpha1plus_x(i, j) / Value;

		};
		double omega2plus_x(int const i, int const j) {

			double const Value = alpha1plus_x(i, j) + alpha2plus_x(i, j) + alpha3plus_x(i, j);

			return alpha2plus_x(i, j) / Value;

		};
		double omega3plus_x(int const i, int const j) {

			double const Value = alpha1plus_x(i, j) + alpha2plus_x(i, j) + alpha3plus_x(i, j);

			return alpha3plus_x(i, j) / Value;

		};

		double omega1minus_x(int const i, int const j) {

			double const Value = alpha1minus_x(i, j) + alpha2minus_x(i, j) + alpha3minus_x(i, j);

			return alpha1minus_x(i, j) / Value;

		};
		double omega2minus_x(int const i, int const j) {

			double const Value = alpha1minus_x(i, j) + alpha2minus_x(i, j) + alpha3minus_x(i, j);

			return alpha2minus_x(i, j) / Value;

		};
		double omega3minus_x(int const i, int const j) {

			double const Value = alpha1minus_x(i, j) + alpha2minus_x(i, j) + alpha3minus_x(i, j);

			return alpha3minus_x(i, j) / Value;

		};


		double omega1plus_y(int const i, int const j) {

			double const Value = alpha1plus_y(i, j) + alpha2plus_y(i, j) + alpha3plus_y(i, j);

			return alpha1plus_y(i, j) / Value;

		};
		double omega2plus_y(int const i, int const j) {

			double const Value = alpha1plus_y(i, j) + alpha2plus_y(i, j) + alpha3plus_y(i, j);

			return alpha2plus_y(i, j) / Value;

		};
		double omega3plus_y(int const i, int const j) {

			double const Value = alpha1plus_y(i, j) + alpha2plus_y(i, j) + alpha3plus_y(i, j);

			return alpha3plus_y(i, j) / Value;

		};

		double omega1minus_y(int const i, int const j) {

			double const Value = alpha1minus_y(i, j) + alpha2minus_y(i, j) + alpha3minus_y(i, j);

			return alpha1minus_y(i, j) / Value;

		};
		double omega2minus_y(int const i, int const j) {

			double const Value = alpha1minus_y(i, j) + alpha2minus_y(i, j) + alpha3minus_y(i, j);

			return alpha2minus_y(i, j) / Value;

		};
		double omega3minus_y(int const i, int const j) {

			double const Value = alpha1minus_y(i, j) + alpha2minus_y(i, j) + alpha3minus_y(i, j);

			return alpha3minus_y(i, j) / Value;

		};

		/*****************************************************************************/
		/*                                                                           */
		/*		LEVEL SET FUNCTION													 */
		/*																			 */
		/*																			 */
		/*		LevelSetInit()		Initial condition for the level-set function phi */
		/*																			 */
		/*		LevelSetReInit()	Transform level-set function to signed distance	 */
		/*							function ( |grad Phi| = 1 )						 */
		/*																			 */
		/*		LevelSetUpdate()	Propagate level-set function in velocity field	 */
		/*																			 */
		/*		LevelSetBoundaryCondition() Update Dirichlet boundary condition		 */
		/*									of the level-set function				 */
		/*                                                                           */
		/*****************************************************************************/
		Eigen::MatrixXd OperatorL0(Eigen::MatrixXd const & Phi) {


			Eigen::MatrixXd result = Phi;


			#pragma omp parallel for schedule(static)
			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny + 1; j++) {


					double const a = (Phi(i, j) - Phi(i - 1, j)) / dx;
					double const b = (Phi(i + 1, j) - Phi(i, j)) / dx;

					double const c = (Phi(i, j) - Phi(i, j - 1)) / dy;
					double const d = (Phi(i, j + 1) - Phi(i, j)) / dy;

					//double const sign = Signum(phiZero(i, j));
					//double const sign = Signum(i, j, phi); 
					//double const sign = Signum(i, j, phiTemp);
					double const sign = Signum(i, j, phiZero);


					double G = 0.0;

					if (phiZero(i, j) > 0.0) {

						double const ap = std::max(a, 0.0);
						double const bm = std::max(-b, 0.0);

						double const cp = std::max(c, 0.0);
						double const dm = std::max(-d, 0.0);

						G = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1.0;

					}
					else if (phiZero(i, j) < 0.0) {

						double const am = std::max(-a, 0.0);
						double const bp = std::max(b, 0.0);

						double const cm = std::max(-c, 0.0);
						double const dp = std::max(d, 0.0);

						G = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1.0;

					}

					result(i, j) = -sign * G;


				}
			}

			return result;

		};
		Eigen::MatrixXd OperatorL1(Eigen::MatrixXd const & Phi) {


			Eigen::MatrixXd result = Phi;


			#pragma omp parallel for schedule(static)
			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny + 1; j++) {


					double const a = (Phi(i, j) - Phi(i - 1, j)) / dx;
					double const b = (Phi(i + 1, j) - Phi(i, j)) / dx;

					double const c = (Phi(i, j) - Phi(i, j - 1)) / dy;
					double const d = (Phi(i, j + 1) - Phi(i, j)) / dy;

					//double const sign = Signum(phiZero(i, j));
					//double const sign = Signum(i, j, phi); 
					//double const sign = Signum(i, j, phiTemp);
					double const sign = Signum(i, j, phiZero);


					double G = 0.0;

					if (phiZero(i, j) > alpha) {

						double const ap = std::max(a, 0.0);
						double const bm = std::max(-b, 0.0);

						double const cp = std::max(c, 0.0);
						double const dm = std::max(-d, 0.0);

						G = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1.0;

					}
					else if (phiZero(i, j) < -alpha) {

						double const am = std::max(-a, 0.0);
						double const bp = std::max(b, 0.0);

						double const cm = std::max(-c, 0.0);
						double const dp = std::max(d, 0.0);

						G = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1.0;

					}

					result(i, j) = -sign * G;


				}
			}

			return result;

		};

		double ArgAbsMin(double const x, double const y) {

			return (std::abs(x) < std::abs(y)) ? (x) : (y);

		};
		double SignFunction(double const x) {

			return (x < 0.0) ? (-1.0) : (+1.0);

		};

		void InitializeValues() {

			//std::cout << std::setprecision(3) << phiZero << std::endl << std::endl << std::endl;
			//std::cout << std::setprecision(3) << phi << std::endl << std::endl << std::endl;

			fixed.setConstant(false);

			/************************************************/
			/*                                              */
			/*		Boundary values are fixed				*/
			/*												*/
			/************************************************/
			/*for (int j = 0; j < ny + 2; j++) {

				fixed(0, j)		 = true;
				fixed(nx + 1, j) = true;

			}
			for (int i = 0; i < nx + 2; i++) {

				fixed(i, 0)		 = true;
				fixed(i, ny + 1) = true;

			}*/
			
			//They are not fixed, do the recomputation for the boundary value
			//change boundarz value to -1 / 1

			//std::cout << std::setprecision(3) << fixed << std::endl << std::endl << std::endl;


			/************************************************/
			/*                                              */
			/*		Initialize values						*/
			/*												*/
			/************************************************/
			for (int i = 0; i < nx + 2; i++)
				for (int j = 0; j < ny + 2; j++)
					phi(i, j) = phiZero(i, j) >= 0.0 ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();


			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {


					double const phi0 = phiZero(i, j);
					double const sign0 = SignFunction(phi0);

					if (i != 0) {
						if (phi0 * phiZero(i - 1, j) <= 0.0) {

							double const t = sign0 * phi0 * dx / (phi0 - phiZero(i - 1, j));

							if (std::abs(phi(i, j)) > std::abs(t))
								phi(i, j) = t;

							fixed(i, j) = true;

						}
					}
					if (i != nx + 1) {
						if (phi0 * phiZero(i + 1, j) <= 0.0) {

							double const t = sign0 * phi0 * dx / (phi0 - phiZero(i + 1, j));

							if (std::abs(phi(i, j)) > std::abs(t))
								phi(i, j) = t;

							fixed(i, j) = true;

						}
					}
					if (j != 0) {
						if (phi0 * phiZero(i, j - 1) <= 0.0) {

							double const t = sign0 * phi0 * dy / (phi0 - phiZero(i, j - 1));

							if (std::abs(phi(i, j)) > std::abs(t))
								phi(i, j) = t;

							fixed(i, j) = true;

						}
					}
					if (j != ny + 1) {
						if (phi0 * phiZero(i, j + 1) <= 0.0) {

							double const t = sign0 * phi0 * dy / (phi0 - phiZero(i, j + 1));

							if (std::abs(phi(i, j)) > std::abs(t))
								phi(i, j) = t;

							fixed(i, j) = true;

						}
					}

				}
			}

			//for (int i = 0; i < nx + 1; i++) {
			//	for (int j = 0; j < ny + 1; j++) {

			//		double const phi0 = phiZero(i, j);
			//		double const phi0x = phiZero(i + 1, j);
			//		double const phi0y = phiZero(i, j + 1);

			//		double const sign0 = SignFunction(phi0);
			//		//double const sign0 = Signum(phi0);

			//		bool const iSign = phi0 * phi0x <= 0.0;
			//		bool const jSign = phi0 * phi0y <= 0.0;


			//		if (iSign) {

			//			double const t = sign0 * phi0 * dx / (phi0 - phi0x);

			//			if (std::abs(phi(i, j)) > std::abs(t))
			//				phi(i, j) = t;
			//			if (std::abs(phi(i + 1, j)) > std::abs(t - sign0 * dx))
			//				phi(i + 1, j) = t - sign0 * dx;

			//			fixed(i, j)		= true;
			//			fixed(i + 1, j) = true;

			//		}

			//		if (jSign) {

			//			double const t = sign0 * phi0 * dy / (phi0 - phi0y);

			//			if (std::abs(phi(i, j)) > std::abs(t))
			//				phi(i, j) = t;
			//			if (std::abs(phi(i, j + 1)) > std::abs(t - sign0 * dy))
			//				phi(i, j + 1) = t - sign0 * dy;

			//			fixed(i, j)		= true;
			//			fixed(i, j + 1) = true;

			//		}

			//	}
			//}

			/*int const I = nx + 1;
			int const J = ny + 1;

			for (int i = 0; i < nx + 1; i++) {

				double const phi0 = phiZero(i, J);
				double const phix = phiZero(i + 1, J);

				double const sign0 = SignFunction(phi0);

				if (phi0 * phix < 0.0) {

					double const t = sign0 * phi0 * dx / (phi0 - phix);

					if (std::abs(phi(i, J)) > std::abs(t))
						phi(i, J) = t;
					if (std::abs(phi(i + 1, J)) > std::abs(t - sign0 * dx))
						phi(i + 1, J) = t - sign0 * dx;

					fixed(i, J) = true;
					fixed(i + 1, J) = true;

				}
			}

			for (int j = 0; j < ny + 1; j++) {

				double const phi0 = phiZero(I, j);
				double const phiy = phiZero(I, j + 1);

				double const sign0 = SignFunction(phi0);

				if (phi0 * phiy < 0.0) {

					double const t = sign0 * phi0 * dy / (phi0 - phiy);

					if (std::abs(phi(I, j)) > std::abs(t))
						phi(I, j) = t;
					if (std::abs(phi(I, j + 1)) > std::abs(t - sign0 * dy))
						phi(I, j + 1) = t - sign0 * dy;

					fixed(I, j) = true;
					fixed(I, j + 1) = true;

				}
			}*/


			//LevelSetBoundaryCondition(phi);

			//std::cout << std::setprecision(3) << phiZero << std::endl << std::endl << std::endl;
			//std::cout << std::setprecision(3) << phi << std::endl << std::endl << std::endl;
			//std::cout << fixed << std::endl;

			//phiStar		= phi;
			//phiStarStar   = phi;
			//phiZero		= phi;
			//phiPrev		= phi;

		};
		void FastSweepMethod() {


			bool changed = true;

			while (changed) {

				changed = false;

				for (int i = 0; i < nx + 2; i++)
					for (int j = 0; j < ny + 2; j++)
						if (!fixed(i, j))
							changed = UpdateCell(i, j) || changed;

				for (int i = 0; i < nx + 2; i++)
					for (int j = ny + 1; j > -1; j--)
						if (!fixed(i, j))
							changed = UpdateCell(i, j) || changed;

				for (int i = nx + 1; i > -1; i--)
					for (int j = ny + 1; j > -1; j--)
						if (!fixed(i, j))
							changed = UpdateCell(i, j) || changed;

				for (int i = nx + 1; i > -1; i--)
					for (int j = 0; j < ny + 2; j++)
						if (!fixed(i, j))
							changed = UpdateCell(i, j) || changed;
				
			}

			//std::cout << fixed << std::endl << std::endl << std::endl;

			//std::cout << std::setprecision(3) << phi << std::endl << std::endl << std::endl;
			//std::cout << phiZero << std::endl << std::endl << std::endl;

			//std::cout << phi << std::endl << std::endl;
			//std::cout << phiZero << std::endl << std::endl << std::endl;

			LevelSetBoundaryCondition(phi);


			phiStar		= phi;
			phiStarStar = phi;
			phiZero		= phi;
			phiPrev		= phi;

		}
		bool UpdateCell(int const i, int const j) {


			double a;
			double b;

			double const originalValue = phi(i, j);

			if (i == 0)				a = phi(i + 1, j);
			else if (i == nx + 1)	a = phi(i - 1, j);
			else					a = ArgAbsMin(phi(i - 1, j), phi(i + 1, j));

			if (j == 0)				b = phi(i, j + 1);
			else if (j == ny + 1)	b = phi(i, j - 1);
			else					b = ArgAbsMin(phi(i, j - 1), phi(i, j + 1));			  

			if (fabs(a) == std::numeric_limits<double>::max() && fabs(b) == std::numeric_limits<double>::max())
				return false;

			double valuesAndSteps[4] = { a, b, dx, dy};
			double newValue = GetNewValue(valuesAndSteps, originalValue);


			phi(i, j) = newValue;

			if (fabs(originalValue - phi(i, j)) > 0.001*dx)
				return true;
			else
				return false;

		};
		double GetNewValue(double ValuesAndSteps[], double const OriginalValue) {


			sortMinims(ValuesAndSteps);

			double const u1 = ValuesAndSteps[0];
			double const u2 = ValuesAndSteps[1];
			double const h1 = ValuesAndSteps[2];
			double const h2 = ValuesAndSteps[3];

			double const newValue1 = u1 + SignFunction(OriginalValue) * h1;
			//double const newValue1 = u1 + Signum(OriginalValue) * h1;

			if (fabs(newValue1) < fabs(u2))
				return ArgAbsMin(OriginalValue, newValue1);
			

			double const norm2 = h1 * h1 + h2 * h2;
			double const amb2  = (u2 - u1)*(u2 - u1);

			double const newValue2 = (u2*h1*h1 + u1 * h2*h2 + SignFunction(OriginalValue) * h1*h2 * sqrt(norm2 - amb2)) / norm2;
			//double const newValue2 = (u2*h1*h1 + u1 * h2*h2 + Signum(OriginalValue) * h1*h2 * sqrt(norm2 - amb2)) / norm2;


			return ArgAbsMin(OriginalValue, newValue2);

		};
		void sortMinims(double pom[]) {

			double temp[4] = { 0.0,0.0,0.0,0.0 };

			if (fabs(pom[0]) <= fabs(pom[1])) {

				temp[0] = pom[0];
				temp[1] = pom[1];

				temp[2] = pom[2];
				temp[3] = pom[3];

			}
			else {

				temp[0] = pom[1];
				temp[1] = pom[0];

				temp[2] = pom[3];
				temp[3] = pom[2];

			}

			for (unsigned int i = 0; i < 4; i++)
				pom[i] = temp[i];
			
			return;

			/*double tmp[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };

			if (fabs(pom[0]) <= fabs(pom[1]) && fabs(pom[1]) <= fabs(pom[2])) {
				tmp[0] = pom[0]; tmp[1] = pom[1]; tmp[2] = pom[2];
				tmp[3] = pom[3]; tmp[4] = pom[4]; tmp[5] = pom[5];

			}
			else if (fabs(pom[0]) <= fabs(pom[2]) && fabs(pom[2]) <= fabs(pom[1])) {
				tmp[0] = pom[0]; tmp[1] = pom[2]; tmp[2] = pom[1];
				tmp[3] = pom[3]; tmp[4] = pom[5]; tmp[5] = pom[4];
			}
			else if (fabs(pom[1]) <= fabs(pom[0]) && fabs(pom[0]) <= fabs(pom[2])) {
				tmp[0] = pom[1]; tmp[1] = pom[0]; tmp[2] = pom[2];
				tmp[3] = pom[4]; tmp[4] = pom[3]; tmp[5] = pom[5];
			}
			else if (fabs(pom[1]) <= fabs(pom[2]) && fabs(pom[2]) <= fabs(pom[0])) {
				tmp[0] = pom[1]; tmp[1] = pom[2]; tmp[2] = pom[0];
				tmp[3] = pom[4]; tmp[4] = pom[5]; tmp[5] = pom[3];
			}
			else if (fabs(pom[2]) <= fabs(pom[0]) && fabs(pom[0]) <= fabs(pom[1])) {
				tmp[0] = pom[2]; tmp[1] = pom[0]; tmp[2] = pom[1];
				tmp[3] = pom[5]; tmp[4] = pom[3]; tmp[5] = pom[4];
			}
			else if (fabs(pom[2]) <= fabs(pom[1]) && fabs(pom[1]) <= fabs(pom[0])) {
				tmp[0] = pom[2]; tmp[1] = pom[1]; tmp[2] = pom[0];
				tmp[3] = pom[5]; tmp[4] = pom[4]; tmp[5] = pom[3];
			}

			for (unsigned int i = 0; i < 6; i++)
			{
				pom[i] = tmp[i];
			}*/
			
		};


		void LevelSetInit() {


			/*for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					double const X = x(i) - 0.06;
					double const Y = y(j) - 0.045;

					double const r = sqrt(X * X + Y * Y);

					phi(i, j) = r <= 0.025 ? +1.0 : -1.0;

				}
			}*/


			/************************************************/
			/*                                              */
			/*		Rising bubble							*/
			/*												*/
			/************************************************/
			/*for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					double const X = x(i) - 0.06;
					double const Y = y(j) - 0.045;

					phi(i, j) = 0.025 - sqrt(X * X + Y * Y);

				}
			}*/


			/************************************************/
			/*                                              */
			/*		Rising bubble to the surface			*/
			/*												*/
			/************************************************/
			phi.setConstant(-1.0);

			//std::cout << 0.5*(x[0] + x[nx+1]) << std::endl;

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					/*double const X = x(i) - 0.06;
					double const Y = y(j) - 0.06;

					phi(i, j) = 0.008 - sqrt(X * X + Y * Y);*/
					//phi(i, j) = 0.02 * 0.02 - (X * X + Y * Y);


					double const X = x(i) - 0.5;
					double const Y = y(j) - 0.75;

					phi(i, j) = 0.15 - sqrt(X * X + Y * Y);

					/*if (sqrt(X * X + Y * Y) <= 0.02)
						phi(i, j) = +1.0;*/

				}
			}
			/*for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					if (y(j) <= 0.027)
						phi(i, j) = +1.0;

				}
			}*/


			/************************************************/
			/*                                              */
			/*		Falling droplet							*/
			/*												*/
			/************************************************/
			/*for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					double const X = x(i) - 0.06;
					double const Y = y(j) - 0.060;

					phi(i, j) = 0.025 - sqrt(X * X + Y * Y);

				}
			}
			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					if (y(j) <= 0.02)
						phi(i, j) = +1.0;

				}
			}*/

			/************************************************/
			/*                                              */
			/*		Lid-driven cavity						*/
			/*												*/
			/************************************************/
			//phi.setConstant(1.0);
			//


			LevelSetBoundaryCondition(phi);

			phiStar		= phi;
			phiStarStar = phi;
			phiZero		= phi;
			phiPrev		= phi;
			
			//LevelSetReInit_0();
			//LevelSetReInit_alpha();
			//LevelSetReInit_OperatorL0();
			//LevelSetReInit_OperatorL1();

			//LevelSetReInit_Zhao();

			//LevelSetExport("D:\\simulations\\levelset\\levelset" + std::to_string(nt) + ".txt");
			


		};
		void LevelSetUpdate() {

			phiPrev = phi;
			phiStar = phi;

			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					phiStar(i, j) = phi(i, j) - dt * UDotGradPhi(i, j);

			phi = phiStar;			


			/*phiStar = phi + dt * OperatorF(phi);
			phiStar = phiStar + dt * OperatorF(phiStar);

			phi = 0.5 * (phi + phiStar);*/

			//phiStar = phi + dt * OperatorF(phi);
			//phiStarStar = phiStar + dt * OperatorF(phiStar);
			//
			//phi = 0.5 * (phi + phiStarStar);

		};
		void LevelSetBoundaryCondition(Eigen::MatrixXd & Phi) {


			/************************************************/
			/*                                              */
			/*		Neumann homogennous GradPhi x N = 0		*/
			/*												*/
			/************************************************/
			for (int j = 0; j < ny + 2; j++) {

				Phi(0, j)		= Phi(1, j);
				Phi(nx + 1, j)  = Phi(nx, j);

			}
			for (int i = 0; i < nx + 2; i++) {

				Phi(i, 0)		= Phi(i, 1);
				Phi(i, ny + 1)  = Phi(i, ny);

			}

			/************************************************/
			/*                                              */
			/*		Neumann homogennous GradPhi x N = +-1	*/
			/*												*/
			/************************************************/
			/*for (int j = 0; j < ny + 2; j++) {

				Phi(0, j)		= Phi(1, j) + SignFunction(phi(1, j))*dx;
				Phi(nx + 1, j)	= Phi(nx, j) + SignFunction(phi(nx, j))*dx;

			}
			for (int i = 0; i < nx + 2; i++) {

				Phi(i, 0)		= Phi(i, 1) + SignFunction(phi(i, 1))*dy;
				Phi(i, ny + 1)	= Phi(i, ny) + SignFunction(phi(i, ny))*dy;

			}*/
			/*for (int j = 0; j < ny + 2; j++) {

				Phi(0, j) = Phi(1, j) + dx;
				Phi(nx + 1, j) = Phi(nx, j) + dx;

			}
			for (int i = 0; i < nx + 2; i++) {

				Phi(i, 0) = Phi(i, 1) + dy;
				Phi(i, ny + 1) = Phi(i, ny) + dy;

			}*/


			/************************************************/
			/*                                              */
			/*		Dirichlet homogenous					*/
			/*												*/
			/************************************************/
			/*for (int j = 0; j < ny + 2; j++) {

				Phi(0, j)	   = -1.0;
				Phi(nx + 1, j) = -1.0;

			}
			for (int i = 0; i < nx + 2; i++) {

				Phi(i, 0)	   = -1.0;
				Phi(i, ny + 1) = -1.0;

			}*/

		};

		void LevelSetReInit_0() {

			phiZero				    = phi;
			Eigen::MatrixXd phiTemp = phi;

			unsigned M	 = 0;
			double Error = 0.0;

			
			do {

				//#pragma omp parallel for schedule(static)
				for (int i = 1; i < nx + 1; i++) {
					for (int j = 1; j < ny + 1; j++) {


						/*double const a = i != 0 ? (phiTemp(i, j) - phiTemp(i - 1, j)) / dx : (phiTemp(i + 1, j) - phiTemp(i, j)) / dx;
						double const b = i != nx + 1 ? (phiTemp(i + 1, j) - phiTemp(i, j)) / dx : (phiTemp(i, j) - phiTemp(i - 1, j)) / dx;

						double const c = j != 0 ? (phiTemp(i, j) - phiTemp(i, j - 1)) / dy : (phiTemp(i, j + 1) - phiTemp(i, j)) / dy;
						double const d = j != ny + 1 ? (phiTemp(i, j + 1) - phiTemp(i, j)) / dy : (phiTemp(i, j) - phiTemp(i, j - 1)) / dy;*/

						double const a = (phiTemp(i, j) - phiTemp(i - 1, j)) / dx;
						double const b = (phiTemp(i + 1, j) - phiTemp(i, j)) / dx;

						double const c = (phiTemp(i, j) - phiTemp(i, j - 1)) / dy;
						double const d = (phiTemp(i, j + 1) - phiTemp(i, j)) / dy;

						double const sign = SignFunction(phiZero(i, j));


						double G = 0.0;

						if (phiZero(i, j) > 0.0) {

							double const ap = std::max(a, 0.0);
							double const bm = std::min(b, 0.0);

							double const cp = std::max(c, 0.0);
							double const dm = std::min(d, 0.0);

							G = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1.0;

						}
						else if (phiZero(i, j) < 0.0) {

							double const am = std::min(a, 0.0);
							double const bp = std::max(b, 0.0);

							double const cm = std::min(c, 0.0);
							double const dp = std::max(d, 0.0);

							G = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1.0;

						}

						//std::cout << G << std::endl;

						phi(i, j) = phiTemp(i, j) - sign * dtau * G;

					}
				}

				//LevelSetExport("D:\\simulations\\levelset\\levelset" + std::to_string(nt) + ".txt");

				LevelSetBoundaryCondition(phi);

				Error = 0;
				M	  = 1;

				for (int i = 1; i < nx + 1; i++) {
					for (int j = 1; j < ny + 1; j++) {

						//if (abs(phiZero(i, j)) < alpha) {
						if (fabs(phiTemp(i, j)) < alpha) {

							M++;
							Error += fabs(phi(i, j) - phiTemp(i, j));

						}
					}
				}

				if (M == 0) 
					Error = INFINITY;
				else		
					Error /= M;


				phiTemp = phi;


			} while (Error > dtau * dx * dy);

		};
		void LevelSetReInit_alpha() {

			phiZero					= phi;
			Eigen::MatrixXd phiTemp = phi;


			unsigned M	 = 0;
			double Error = 0.0;


			do {

				#pragma omp parallel for schedule(static)
				for (int i = 1; i < nx + 1; i++) {
					for (int j = 1; j < ny + 1; j++) {


						double const a = (phiTemp(i, j) - phiTemp(i - 1, j)) / dx;
						double const b = (phiTemp(i + 1, j) - phiTemp(i, j)) / dx;

						double const c = (phiTemp(i, j) - phiTemp(i, j - 1)) / dy;
						double const d = (phiTemp(i, j + 1) - phiTemp(i, j)) / dy;

						/*
						double a, b, c, d;

						if (i == 1 || i == 2 || i == nx || i == nx - 1)
							a = (phi(i, j) - phi(i - 1, j)) / dx;
						else {

							double const Value1 = omega1minus_x(i, j)*(q1minus_x(i, j) / 3.0 - 7.0 * q2minus_x(i, j) / 6.0 + 11.0 * q3minus_x(i, j) / 6.0);
							double const Value2 = omega2minus_x(i, j)*(-q2minus_x(i, j) / 6.0 + 5.0 * q3minus_x(i, j) / 6.0 + q4minus_x(i, j) / 3.0);
							double const Value3 = omega3minus_x(i, j)*(q3minus_x(i, j) / 3.0 + 5.0 * q4minus_x(i, j) / 6.0 - q5minus_x(i, j) / 6.0);

							a = Value1 + Value2 + Value3;

						}
						if (i == 1 || i == 2 || i == nx - 1 || i == nx)
							b = (phi(i + 1, j) - phi(i, j)) / dx;
						else {


							double const Value1 = omega1plus_x(i, j)*(q1plus_x(i, j) / 3.0 - 7.0 * q2plus_x(i, j) / 6.0 + 11.0 * q3plus_x(i, j) / 6.0);
							double const Value2 = omega2plus_x(i, j)*(-q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0);
							double const Value3 = omega3plus_x(i, j)*(q3plus_x(i, j) / 3.0 + 5.0 * q4plus_x(i, j) / 6.0 - q5plus_x(i, j) / 6.0);

							b = Value1 + Value2 + Value3;

						}
						if (j == 1 || j == 2 || j == ny - 1 || j == ny)
							c = (phi(i, j) - phi(i, j - 1)) / dy;
						else {

							double const Value1 = omega1minus_y(i, j)*(q1minus_y(i, j) / 3.0 - 7.0 * q2minus_y(i, j) / 6.0 + 11.0 * q3minus_y(i, j) / 6.0);
							double const Value2 = omega2minus_y(i, j)*(-q2minus_y(i, j) / 6.0 + 5.0 * q3minus_y(i, j) / 6.0 + q4minus_y(i, j) / 3.0);
							double const Value3 = omega3minus_y(i, j)*(q3minus_y(i, j) / 3.0 + 5.0 * q4minus_y(i, j) / 6.0 - q5minus_y(i, j) / 6.0);

							c = Value1 + Value2 + Value3;
						}
						if (j == 1 || j == 2 || j == ny - 1 || j == ny)
							d = (phi(i, j + 1) - phi(i, j)) / dy;
						else {

							double const Value1 = omega1plus_y(i, j)*(q1plus_y(i, j) / 3.0 - 7.0 * q2plus_y(i, j) / 6.0 + 11.0 * q3plus_y(i, j) / 6.0);
							double const Value2 = omega2plus_y(i, j)*(-q2plus_y(i, j) / 6.0 + 5.0 * q3plus_y(i, j) / 6.0 + q4plus_y(i, j) / 3.0);
							double const Value3 = omega3plus_y(i, j)*(q3plus_y(i, j) / 3.0 + 5.0 * q4plus_y(i, j) / 6.0 - q5plus_y(i, j) / 6.0);

							d = Value1 + Value2 + Value3;

						}
						*/


						double const sign = Signum(phiZero(i, j));
						//double const sign = SignFunction(phiZero(i, j));

						double G = 0.0;

						if (phiZero(i, j) > alpha) {

							double const ap = std::max(a, 0.0);
							double const bm = std::min(b, 0.0);

							double const cp = std::max(c, 0.0);
							double const dm = std::min(d, 0.0);

							G = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1.0;

						}
						else if (phiZero(i, j) < -alpha) {

							double const am = std::min(a, 0.0);
							double const bp = std::max(b, 0.0);

							double const cm = std::min(c, 0.0);
							double const dp = std::max(d, 0.0);

							G = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1.0;

						}

						phi(i, j) = phiTemp(i, j) - dtau * sign * G;

					}
				}

				Error = 0;
				M = 0;

				for (int i = 1; i < nx + 1; i++) {
					for (int j = 1; j < ny + 1; j++) {

						if (abs(phiTemp(i, j)) < alpha) {

							M++;
							Error += abs(phi(i, j) - phiTemp(i, j));

						}
					}
				}

				if (M == 0) Error = INFINITY;
				else		Error /= M;

				phiTemp = phi;


			} while (Error > dtau * dx * dy);

		};
		void LevelSetReInit_OperatorL0() {

			phiZero					= phi;
			Eigen::MatrixXd phiTemp = phi;

			double Error = 0.0;
			double tau   = 0.0;

			do {

				Error = 0;

				phiStar		= phi + dtau * (OperatorL0(phi));
				phiStarStar = phi + dtau * (OperatorL0(phi) + OperatorL0(phiStar)) / 4.0;
				phi			= phi + dtau * (OperatorL0(phi) + 4.0 * OperatorL0(phiStarStar) + OperatorL0(phiStar)) / 6.0;

				for (int i = 0; i < nx + 2; i++)
					for (int j = 0; j < ny + 2; j++)
						Error += abs(phiTemp(i, j) - phi(i, j)) * dx * dy;

				tau += dtau;

			} while (tau < 5.0*alpha);
			//while (Error > 1e-5 * dtau);; (tt < 2.0*alpha)

			std::cout << "Reinitialization error: " << Error << std::endl << std::endl;

		};
		void LevelSetReInit_OperatorL1() {

			phiZero				    = phi;
			Eigen::MatrixXd phiTemp = phi;

			double Error = 0.0;
			double tau   = 0.0;

			do {

				Error = 0;

				phiStar		= phi + dtau * (OperatorL1(phi));
				phiStarStar = phi + dtau * (OperatorL1(phi) + OperatorL1(phiStar)) / 4.0;
				phi			= phi + dtau * (OperatorL1(phi) + 4.0 * OperatorL1(phiStarStar) + OperatorL1(phiStar)) / 6.0;

				for (int i = 0; i < nx + 2; i++)
					for (int j = 0; j < ny + 2; j++)
						Error += abs(phiTemp(i, j) - phi(i, j)) * dx * dy;

				tau += dtau;

			} while (tau < 5.0*alpha);
			//while (Error > 1e-5 * dtau);; (tt < 2.0*alpha)

			std::cout << "Reinitialization error: " << Error << std::endl << std::endl;

		};
		void LevelSetReInit_Zhao() {

			phiZero = phi;

			InitializeValues();
			FastSweepMethod();

		};
				
		public:

		void EOC() {


			Eigen::MatrixXd phiAnalytical((nx + 2), (ny + 2));

			phi.setConstant(-1.0);

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					//double const X = x(i) - 0.06;
					//double const Y = y(j) - 0.0375;

					//phi(i, j) = 0.02 * 0.02 - (X * X + Y * Y);
					//phi(i, j) = 0.02 - sqrt(X * X + Y * Y);

					double const X = x(i) - 0.0;
					double const Y = y(j) - 0.0;
					phi(i, j) = -0.75*0.75 + (X * X + Y * Y);

					/*double const X = x(i) - 0.0;
					double const Y = y(j) - 0.0;
					if (sqrt(X*X + Y * Y) < 0.8)
						phi(i, j) = sin(2.0*Pi *sqrt(X*X + Y * Y) / 0.4);
					else
						phi(i, j) = sqrt(X*X + Y * Y) - 0.8;*/

				}
			}


			LevelSetReInit_Zhao();

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					/*
					double const X = x(i) - 0.06;
					double const Y = y(j) - 0.0375;

					phiAnalytical(i, j) = 0.02 - sqrt(X * X + Y * Y);*/

					double const X = x(i) - 0.0;
					double const Y = y(j) - 0.0;
					phiAnalytical(i, j) = -0.75 + sqrt(X * X + Y * Y);

					/*double const X = x(i) - 0.0;
					double const Y = y(j) - 0.0;
					if (sqrt(X*X + Y * Y) < 0.8)
						phiAnalytical(i, j) = sin(2.0*Pi *sqrt(X*X + Y * Y) / 0.4);
					else
						phiAnalytical(i, j) = (X*X + Y * Y) - 0.8*0.8;*/
				}
			}

			//std::ofstream txtfile2;
			//txtfile2.open("C:\\Users\\pgali\\Desktop\\levelset\\levelset_" + std::to_string(nx) + ".txt");
			//for (int i = 0; i < nx + 2; i++)
			//	for (int j = 0; j < ny + 2; j++)
			//		txtfile2 << x(i) << " " << y(j) << " " << phi(i, j) << std::endl;
			//txtfile2.close();


			double normMax = 0.0;
			double normL1 = 0.0;
			double normL2 = 0.0;

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {

					double const absDifference = fabs(phi(i, j) - phiAnalytical(i, j));

					normL1 += absDifference * dx * dy;
					normL2 += absDifference * absDifference * dx * dy;

					if (absDifference > normMax)
						normMax = absDifference;

				}
			}

			std::ofstream txtfile;
			txtfile.open("C:\\Users\\pgali\\Desktop\\levelset\\eocphi" + std::to_string(nx) +  ".txt");

			txtfile << "#L1 L2 MAX" << std::endl;

			txtfile << std::setprecision(20) << normL1 << std::endl;
			txtfile << std::setprecision(20) << sqrt(normL2) << std::endl;
			txtfile << std::setprecision(20) << normMax << std::endl;

			txtfile.close();

		};

		/*****************************************************************************/
		/*                                                                           */
		/*		POISSON EQUATION FOR PRESSURES										 */
		/*																			 */
		/*																			 */
		/*		PoissonSystemInit()	Allocate memory for the Poisson system solver.	 */
		/*							Initialize position	of elements in the matrix	 */
		/*							and extracts sparsity pattern of the matrix 	 */
		/*							used in Eigen solver BiCGStab					 */
		/*																			 */
		/*		PoissonSystemAssembly()		Update elements in the matrix			 */
		/*																			 */
		/*		PoissonSystemSolve()		Solve the Poisson system for pressure	 */
		/*									correction with iterative solver 		 */
		/*									BiCGStab using diagonal preconditioner.	 */
		/*                                                                           */
		/*****************************************************************************/
		void PoissonSystemInit() {


			double const dummyValue = 0.0;

			/************************************************/
			/*                                              */
			/*		Using CSR sparse matrix format			*/
			/*												*/
			/************************************************/
			ptr.clear();
			col.clear();
			val.clear();


			/************************************************/
			/*                                              */
			/*		5-point stencil							*/
			/*                                              */
			/*      At most (nx+2)*(ny+2) * 5 nonzero       */
			/*		elements.								*/
			/*                                              */
			/************************************************/

			ptr.reserve((nx + 2) * (ny + 2) + 1);
			col.reserve((nx + 2) * (ny + 2) * 5);
			val.reserve((nx + 2) * (ny + 2) * 5);

			ptr.push_back(0);


			/************************************************/
			/*                                              */
			/*		Allocate memory for the matrix			*/
			/*												*/
			/************************************************/
			A.resize((nx + 2) * (ny + 2), (nx + 2) * (ny + 2));
			A.reserve((nx + 2) * (ny + 2) * 5);


			/************************************************/
			/*                                              */
			/*		Allocate memory for the right-hand		*/
			/*		side of the system and the solution		*/
			/*												*/
			/************************************************/
			pressureRightHandSide.resize((nx + 2) * (ny + 2));
			pressureSolution	 .resize((nx + 2) * (ny + 2));


			int k = 0;

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {


					if (i == 1 && j == 1) {

						col.push_back(k);
						val.push_back(dummyValue);
						ptr.push_back(col.size());

						k++;

						continue;

					}

					/************************************************/
					/*                                              */
					/*		Neumann boundary conditions				*/
					/*												*/
					/************************************************/
					if (i == 0) {

						col.push_back(k); 
						col.push_back(k + (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (i == nx + 1) {

						col.push_back(k);
						col.push_back(k - (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == 0) {

						col.push_back(k);
						col.push_back(k + 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == ny + 1) {

						col.push_back(k);
						col.push_back(k - 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					/************************************************/
					/*                                              */
					/*		Inner cells								*/
					/*												*/
					/************************************************/
					else {

						col.push_back(k - (ny + 2));
						col.push_back(k - 1);
						col.push_back(k);
						col.push_back(k + 1);
						col.push_back(k + (ny + 2));
						
						val.push_back(dummyValue);
						val.push_back(dummyValue);
						val.push_back(dummyValue);
						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}

					k++;

					ptr.push_back(col.size());

				}
			}


			/************************************************/
			/*                                              */
			/*		Transform CSR format into Coordinate	*/
			/*		format which is used by Eigen 			*/
			/*												*/
			/*		The transformation uses 'slice'			*/
			/*												*/
			/************************************************/
			for (int i = 0; i < (nx + 2) * (ny + 2); i++) {

				int const row_start = ptr[i];
				int const row_end	= ptr[i + 1];

				for (int col_index = row_start; col_index < row_end; col_index++) {

					int const j		   = col[col_index];
					double const value = val[col_index];

					triplet.push_back(Eigen::Triplet<double>(i, j, value));

				}
			}


			/************************************************/
			/*                                              */
			/*		Assemble matrix from triplets			*/
			/*												*/
			/************************************************/
			A.setFromTriplets(triplet.begin(), triplet.end());

			triplet.clear();


			/************************************************/
			/*                                              */
			/*		Analyze sparsity pattern of the system	*/
			/*		matrix. It is used by the solver.		*/
			/*												*/
			/************************************************/
			poissonSystemSolverBiCGSTAB.analyzePattern(A);

		};
		void PoissonSystemAssembly() {


			double const idx2 = 1.0 / (dx * dx);
			double const idy2 = 1.0 / (dy * dy);

			/************************************************/
			/*                                              */
			/*		Clear old values						*/
			/*												*/
			/************************************************/
			val		.clear();
			triplet	.clear();


			/************************************************/
			/*                                              */
			/*		Assemble the matrix in CSR format		*/
			/*												*/
			/************************************************/
			int k = 0;

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {


					/************************************************/
					/*                                              */
					/*		Because there are only Neumann boundary */
					/*		conditions prescribed, there is no 		*/
					/*      unique solution. Therefore the matrix   */
					/*      is singular and thus we need to assign  */
					/*      arbitrary value to one unknown. For     */
					/*      instance we set p11 = 0                 */
					/*                                              */
					/************************************************/
					if (i == 1 && j == 1) {

						val.push_back(1.0);
						pressureRightHandSide[k] = 0.0;

						k++;

						continue;

					}						
					
					/************************************************/
					/*                                              */
					/*		Neumann boundary conditions				*/
					/*												*/
					/************************************************/
					if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {

						val.push_back(1.0);
						val.push_back(-1.0);
						pressureRightHandSide[k] = 0.0;

					}
					/************************************************/
					/*                                              */
					/*		Inner cells								*/
					/*												*/
					/************************************************/
					else {

						/*double const phi_l = 0.5*(phiPrev(i, j) + phiPrev(i - 1, j));
						double const phi_d = 0.5*(phiPrev(i, j) + phiPrev(i, j - 1));
						double const phi_u = 0.5*(phiPrev(i, j) + phiPrev(i, j + 1));
						double const phi_r = 0.5*(phiPrev(i, j) + phiPrev(i + 1, j));*/

						double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
						double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
						double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
						double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

						double const irho_l = 1.0 / rho(phi_l);
						double const irho_d = 1.0 / rho(phi_d);
						double const irho_u = 1.0 / rho(phi_u);
						double const irho_r = 1.0 / rho(phi_r);

						double const val_c = idx2 * irho_l + idx2 * irho_r + idy2 * irho_d + idy2 * irho_u;


						val.push_back(idx2 * irho_l);
						val.push_back(idy2 * irho_d);
						val.push_back(-val_c);
						val.push_back(idy2 * irho_u);
						val.push_back(idx2 * irho_r);


						double const dudx = (ustar(i, j) - ustar(i - 1, j)) / dx;
						double const dvdy = (vstar(i, j) - vstar(i, j - 1)) / dy;

						pressureRightHandSide[k] = (dudx + dvdy) / dt;
						//pressureRightHandSide[k] = (PoissonQ(i, j) - PoissonQ(i - 1, j)) / dx + (PoissonR(i, j) - PoissonR(i, j - 1)) / dy;
						//pressureRightHandSide[k] = (Q(i, j) - Q(i - 1, j)) / dx + (R(i, j) - R(i, j - 1)) / dy;

					}

					k++;

				}
			}


			/************************************************/
			/*                                              */
			/*		Transform CSR format into coordinate	*/
			/*		format which is used by Eigen			*/
			/*												*/
			/*		The transformation uses 'slice'			*/
			/*												*/
			/************************************************/
			for (int i = 0; i < (nx + 2)*(ny + 2); i++) {

				int const row_start = ptr[i];
				int const row_end   = ptr[i + 1];

				for (int col_index = row_start; col_index < row_end; col_index++) {

					int const j		   = col[col_index];
					double const value = val[col_index];

					triplet.push_back(Eigen::Triplet<double>(i, j, value));

				}
			}


			/************************************************/
			/*                                              */
			/*		Update values in the system matrix		*/
			/*												*/
			/************************************************/
			A.setFromTriplets(triplet.begin(), triplet.end());


		};

		void PoissonSystemSolve_iteration() {


			double const idx2 = 1.0 / (dx * dx);
			double const idy2 = 1.0 / (dy * dy);

			Eigen::MatrixXd pt = p;

			for (int it = 1; it < 30000; it++) {


				/************************************************/
				/*                                              */
				/*		Neumann boundary conditions				*/
				/*												*/
				/************************************************/
				for (int j = 0; j < ny + 2; j++) {

					pt(0, j)	  = pt(1, j);
					pt(nx + 1, j) = pt(nx, j);

				}
				for (int i = 0; i < nx + 2; i++) {

					pt(i, 0)	  = pt(i, 1);
					pt(i, ny + 1) = pt(i, ny);

				}

				//#pragma omp parallel for schedule(static)
				for (int i = 1; i < nx + 1; i++) {
					for (int j = 1; j < ny + 1; j++) {


						double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
						double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
						double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
						double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

						double const irho_l = 1.0 / rho(phi_l);
						double const irho_d = 1.0 / rho(phi_d);
						double const irho_u = 1.0 / rho(phi_u);
						double const irho_r = 1.0 / rho(phi_r);


						double const A = (idx2 * irho_l + idx2 * irho_r + idy2 * irho_d + idy2 * irho_u);

						double const Q = (1.0 / dt) * (ustar(i, j) - ustar(i - 1, j)) / dx;
						double const R = (1.0 / dt) * (vstar(i, j) - vstar(i, j - 1)) / dy;

						pt(i, j) = -(Q + R) / A + ( idx2 * irho_l * pt(i - 1, j) + idx2 * irho_r * pt(i + 1, j) +
													idy2 * irho_d * pt(i, j - 1) + idy2 * irho_u * pt(i, j + 1) ) / A;

					}
				}

				double Error = 0.0;

				Eigen::MatrixXd const temp = pt - p;

				for (int i = 1; i < nx + 1; i++)
					for (int j = 1; j < ny + 1; j++)
						Error += std::abs(temp(i, j));

				if (Error < DBL_EPSILON) {

					std::cout << "Pressure system error: "<< Error << std::endl;
					break;

				}

				p = pt;

			}

		};
		void PoissonSystemSolve_matrixassemble() {


			/************************************************/
			/*                                              */
			/*		Set initial guess for the pressures		*/
			/*		on the previous time level				*/
			/*												*/
			/************************************************/
			for (int i = 0; i < nx + 2; i++)
				for (int j = 0; j < ny + 2; j++)
					pressureSolution[i*(ny + 2) + j] = p(i, j);


			/************************************************/
			/*                                              */
			/*		Solve the Poisson equation		 		*/
			/*												*/
			/************************************************/
			PoissonSystemAssembly();

			//std::cout << pressureRightHandSide << std::endl;

			poissonSystemSolverBiCGSTAB.factorize(A);

			pressureSolution = poissonSystemSolverBiCGSTAB.solveWithGuess(pressureRightHandSide, pressureSolution);



			/************************************************/
			/*                                              */
			/*		Update Pressures				 		*/
			/*												*/
			/************************************************/
			for (int i = 0; i < nx + 2; i++)
				for (int j = 0; j < ny + 2; j++)
					p(i, j) = pressureSolution[i*(ny + 2) + j];

			std::cout << "Pressure solution error: " << poissonSystemSolverBiCGSTAB.error() << std::endl;

		};


		double PoissonQ(int const i, int const j) {


			double const UDivU_x = UDotDivU_x(i, j);
			double const DivTau_x = DivDotTau_x(i, j);
			double const surfaceTension_x = SurfaceTensionForce_x(i, j);

			double const phi_r = 0.5 *(phi(i + 1, j) + phi(i, j));

			double const val_x = dt * UDivU_x + (dt / rho(phi_r)) * (DivTau_x + surfaceTension_x) + dt * gx;

			return val_x;

		};
		double PoissonR(int const i, int const j) {


			double const UDivU_y = UDotDivU_y(i, j);
			double const DivTau_y = DivDotTau_y(i, j);
			double const surfaceTension_y = SurfaceTensionForce_y(i, j);

			double const phi_u = 0.5 *(phi(i, j + 1) + phi(i, j));

			double const val_y = dt * UDivU_y + (dt / rho(phi_u)) * (DivTau_y + surfaceTension_y) + dt * gy;

			return val_y;
	
		}

		/*****************************************************************************/
		/*                                                                           */
		/*		VELOCITIES															 */
		/*																			 */
		/*																			 */
		/*		VelocityUpdate()	Get intermediate velocities ustar, vstar		 */
		/*																			 */
		/*		VelocityCorrect()	Correct velocities with the pressure gradient  	 */
		/*							so they are divergence-free (Div dot U = 0)		 */
		/*																			 */
		/*		VelocityBoundaryCondition()	Update Dirichlet boundary conditions for */
		/*									the velocities u,v						 */
		/*																			 */
		/*****************************************************************************/
		void VelocityUpdate() {

			if (nt < 400) {

				for (int i = 0; i < nx + 1; i++)
					for (int j = 0; j < ny + 2; j++)
						u(i, j) = sin(Pi * x(i))*sin(Pi * x(i)) * sin(2.0 * Pi * y(j));

				for (int i = 0; i < nx + 2; i++)
					for (int j = 0; j < ny + 1; j++)
						v(i, j) = -sin(Pi * y(j))*sin(Pi * y(j)) * sin(2.0 * Pi * x(i));

			}
			else {

				for (int i = 0; i < nx + 1; i++)
					for (int j = 0; j < ny + 2; j++)
						u(i, j) = -sin(Pi * x(i))*sin(Pi * x(i)) * sin(2.0 * Pi * y(j));

				for (int i = 0; i < nx + 2; i++)
					for (int j = 0; j < ny + 1; j++)
						v(i, j) = sin(Pi * y(j))*sin(Pi * y(j)) * sin(2.0 * Pi * x(i));

			}

			
			Eigen::MatrixXd LnU(nx + 1, ny + 2);
			Eigen::MatrixXd LnV(nx + 2, ny + 1);

			LnU.setZero();
			LnV.setZero();

			
			/*
			for (int i = 1; i < nx; i++) {
				for (int j = 1; j < ny + 1; j++) {

					double const UDivU_x		  = UDotDivU_x(i, j);
					double const DivTau_x		  = DivDotTau_x(i, j);
					double const surfaceTension_x = SurfaceTensionForce_x(i, j);

					double const phi_r = 0.5 *(phi(i + 1, j) + phi(i, j));

					//double const rhoo = 0.5 *(rho(phi(i + 1, j)) + rho(phi(i, j)));

					ustar(i, j) = u(i, j) - dt * UDivU_x + (dt / rho(phi_r)) * (DivTau_x + surfaceTension_x) + dt * gx;

				}
			}

			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny; j++) {

					double const UDivU_y = UDotDivU_y(i, j);
					double const DivTau_y = DivDotTau_y(i, j);
					double const surfaceTension_y = SurfaceTensionForce_y(i, j);

					double const phi_u = 0.5 *(phi(i, j + 1) + phi(i, j));

					//double const rhoo = 0.5 *(rho(phi(i, j + 1)) + rho(phi(i, j)));

					vstar(i, j) = v(i, j) - dt * UDivU_y + (dt / rho(phi_u)) * (DivTau_y + surfaceTension_y) + dt * gy;

				}
			}
			*/
			

			
			//for (int i = 1; i < nx; i++) {
			//	for (int j = 1; j < ny + 1; j++) {

			//		double const UDivU_x = UDotDivU_x(i, j);
			//		double const DivTau_x = DivDotTau_x(i, j);
			//		double const surfaceTension_x = SurfaceTensionForce_x(i, j);

			//		double const phi_r = 0.5 *(phi(i + 1, j) + phi(i, j));

			//		//double const rhoo = 0.5 *(rho(phi(i + 1, j)) + rho(phi(i, j)));

			//		LnU(i, j) = -UDivU_x + (1.0 / rho(phi_r)) * (DivTau_x + surfaceTension_x) + gx;

			//		ustar(i, j) = u(i,j) + dt * 0.5 * ((dt + 2.0 * dtprev) / dtprev * LnU(i, j) - dt / dtprev * prevU(i, j));

			//	}
			//}

			//for (int i = 1; i < nx + 1; i++) {
			//	for (int j = 1; j < ny; j++) {

			//		double const UDivU_y		  = UDotDivU_y(i, j);
			//		double const DivTau_y		  = DivDotTau_y(i, j);
			//		double const surfaceTension_y = SurfaceTensionForce_y(i, j);

			//		double const phi_u = 0.5 *(phi(i, j + 1) + phi(i, j));

			//		//double const rhoo = 0.5 *(rho(phi(i, j + 1)) + rho(phi(i, j)));

			//		LnV(i, j) = -UDivU_y + (1.0 / rho(phi_u)) * (DivTau_y + surfaceTension_y) + gy;

			//		vstar(i, j) = v(i,j) + dt * 0.5 * ((dt + 2.0 * dtprev) / dtprev * LnV(i, j) - dt / dtprev * prevV(i, j));

			//	}
			//}

			prevU = LnU;
			prevV = LnV;


		};
		void VelocityCorrect() {


			for (int i = 0; i < nx + 1; i++)
				for (int j = 0; j < ny + 2; j++) {

					double const phi_r = 0.5 *(phi(i + 1, j) + phi(i, j));

					u(i, j) = ustar(i, j) - (dt / rho(phi_r)) * (p(i + 1, j) - p(i, j)) / dx;

				}
					

			for (int i = 0; i < nx + 2; i++)
				for (int j = 0; j < ny + 1; j++) {

					double const phi_u = 0.5 *(phi(i, j + 1) + phi(i, j));
					
					v(i, j) = vstar(i, j) - (dt / rho(phi_u)) * (p(i, j + 1) - p(i, j)) / dy;

				}
					

		};

		void VelocityBoundaryCondition() {


			/************************************************/
			/*                                              */
			/*		Time on the new time level n+1			*/
			/*												*/
			/************************************************/
			double const t = (nt + 1) * dt;


			/************************************************/
			/*                                              */
			/*		Tangential velocity						*/
			/*												*/
			/************************************************/
			/************************************************/
			/*                                              */
			/*		Uses 'ghost' values outside the domain	*/
			/*		computed from averaging using boundary	*/
			/*		value									*/
			/*												*/
			/*		u_bc    = (u_last + u_ghost) / 2		*/
			/*		u_ghost = 2 * u_bc - u_last				*/
			/*												*/
			/************************************************/
			for (int i = 0; i < nx + 1; i++) {

				u(i, 0)		 = 2.0 * VelocitySouth(x(i), t) - u(i, 1);
				u(i, ny + 1) = 2.0 * VelocityNorth(x(i), t) - u(i, ny);

			}
			for (int j = 0; j < ny + 1; j++) {

				v(0, j)		 = 2.0 * VelocityWest(y(j), t) - v(1, j);
				v(nx + 1, j) = 2.0 * VelocityEast(y(j), t) - v(nx, j);

			}


			/************************************************/
			/*                                              */
			/*		Normal velocity							*/
			/*												*/
			/************************************************/
			for (int j = 0; j < ny + 2; j++) {

				//u(0, j) = -1.0 * SignFunction();
				//u(nx, j) = 0.0;

				u(0, j)  = u(1, j);
				u(nx, j) = u(nx - 1, j);

			}
			for (int i = 0; i < nx + 2; i++) {

				//v(i, 0) = 0.0;
				//v(i, ny) = 0.0;

				//v(i, 0) = -dt * gy;
				//v(i, ny) = -dt * gy;

				v(i, 0)  = v(i, 1);
				v(i, ny) = v(i, ny - 1);

			}


		};


		/*double Q(int const i, int const j) {

			Eigen::Vector2d const Force = GravityForce(i, j) - SurfaceTensionForce(i, j) / rho(phi(i, j));

			double const UDivU_x = UDotDivU_x(i, j);
			double const DivTau_x = DivDotTau_x(i, j);

			double const phi_r = 0.5 *(phi(i + 1, j) + phi(i, j));

			

			return dt * (DivTau_x / rho(phi_r) - UDivU_x + Force(0));

		};
		double R(int const i, int const j) {

			Eigen::Vector2d const Force = GravityForce(i, j) - SurfaceTensionForce(i, j) / rho(phi(i, j));

			double const UDivU_y = UDotDivU_y(i, j);
			double const DivTau_y = DivDotTau_y(i, j);

			double const phi_u = 0.5 *(phi(i, j + 1) + phi(i, j));



			return dt * (DivTau_y / rho(phi_u) - UDivU_y + Force(1));

		};*/

		/************************************************/
		/*                                              */
		/*		Tangential Boundary conditions			*/
		/*												*/
		/************************************************/
		double VelocityWest(double const y, double const t) {

			return 0.0;
			//return dt * gy;

		};
		double VelocityEast(double const y, double const t) {

			return 0.0;
			//return dt * gy;

		};
		double VelocityNorth(double const x, double const t) {

			return 0.0;
			//return 1.0;

		};
		double VelocitySouth(double const x, double const t) {

			return 0.0;

		};


		/************************************************/
		/*                                              */
		/*		Normal Boundary conditions				*/
		/*												*/
		/************************************************/

		


	public:

		
		void setTimeLevel(int const n) {

			nt = n;

		};
		void setTimeIncrement(double const dtime) {

			dt = dtime;
			dtprev = dt;

		};
		double getTime() const {

			return nt * dt;

		};

		void setAdaptiveTimeStep(bool const b) {

			AdaptiveTimeStep = b;

		};

		void checkTimeStep() {


			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					kappa(i, j) = DivGradPhi(i, j);


			double const maxu1 = std::abs(u.maxCoeff());
			double const maxu2 = std::abs(u.minCoeff());
			double const maxv1 = std::abs(v.maxCoeff());
			double const maxv2 = std::abs(v.minCoeff());

			double maxU = std::max(maxu1, maxu2);
			double maxV = std::max(maxv1, maxv2);

			double const maxkappa1 = std::abs(kappa.maxCoeff());
			double const maxkappa2 = std::abs(kappa.minCoeff());

			double const maxKappa = std::max(maxkappa1, maxkappa2);


			/************************************************/
			/*                                              */
			/*		Convective restriction					*/
			/*												*/
			/************************************************/
			double const convectiveU = maxU / dx;
			double const convectiveV = maxV / dy;

			/************************************************/
			/*                                              */
			/*		Convective restriction gravity			*/
			/*												*/
			/************************************************/
			double const gravityU = 4.0 * std::abs(gx) / dx;
			double const gravityV = 4.0 * std::abs(gy) / dy;

			/************************************************/
			/*                                              */
			/*		Viscous restriction						*/
			/*												*/
			/************************************************/
			double const viscousBound = std::max(mu1 / rho1, mu2 / rho2) * (2.0 / (dx*dx) + 2.0 / (dy*dy));

			/************************************************/
			/*                                              */
			/*		Surface tension restriction				*/
			/*												*/
			/************************************************/
			double const tensionU = 8.0 * maxKappa * sigma / (alpha * (rho2 + rho1) * dx * dx);
			double const tensionV = 8.0 * maxKappa * sigma / (alpha * (rho2 + rho1) * dy * dy);


			double const valU = convectiveU + viscousBound;
			double const valV = convectiveV + viscousBound;

			double const boundU = 1.0 / (valU + sqrt(valU*valU + gravityU + tensionU));
			double const boundV = 1.0 / (valV + sqrt(valV*valV + gravityV + tensionV));


			double const Xi = 0.5;

			dt = Xi * 2.0 * std::min(boundU, boundV);



		};

		void solve() {

			

			VelocityBoundaryCondition();
			
			//std::cout << std::setprecision(3) << phi << std::endl << std::endl;

			/************************************************/
			/*                                              */
			/*		Get intermediate non-divergence free	*/
			/*		velocity U^star							*/
			/*												*/
			/************************************************/
			VelocityUpdate();
			/*VelocityBoundaryCondition();*/


			/************************************************/
			/*                                              */
			/*		Level-set update and reinitialization	*/
			/*												*/
			/************************************************/
			LevelSetUpdate();

			LevelSetBoundaryCondition(phi);

			//LevelSetReInit_0();
			//LevelSetReInit_alpha();
			//LevelSetReInit_OperatorL0();
			//LevelSetReInit_OperatorL1();
			LevelSetReInit_Zhao();

			LevelSetBoundaryCondition(phi);
			
			//std::cout << std::setprecision(3) << phi << std::endl << std::endl;
			//std::cout << phiZero << std::endl << std::endl << std::endl;
			

			/************************************************/
			/*                                              */
			/*		Get pressure corrections solving		*/
			/*		Poisson equation						*/
			/*												*/
			/************************************************/
			//PoissonSystemSolve_iteration();
			PoissonSystemSolve_matrixassemble();


			/************************************************/
			/*                                              */
			/*		Correct the velocities so they are		*/
			/*		diveregence-free						*/
			/*												*/
			/************************************************/
			VelocityCorrect();
			//VelocityBoundaryCondition();



			if (AdaptiveTimeStep) {

				dtprev = dt;
				checkTimeStep();

			}
			


		};

		void clear() {

			NumberOfPointsX = 0;
			NumberOfPointsY = 0;

			//velocityU.resize(0);
			//velocityV.resize(0);
			//pressureP.resize(0);

			x.resize(0);
			y.resize(0);

			dx = 0;
			dy = 0;
			dt = 0;
		};

		void LevelSetExport(std::string const fileName) {


			std::ofstream txtfile;
			txtfile.open(fileName);


			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					txtfile << x(i) << " " << y(j) << " " << phi(i, j) << std::endl;

			//for (int i = 0; i < nx + 1; i++)
			//	for (int j = 0; j < ny + 1; j++)
			//		txtfile << x(i) << " " << y(j) << " " << 0.25 *( phi(i, j) + phi(i + 1, j) + phi(i, j + 1) + phi(i+1, j+1)) << std::endl;

			//Eigen::MatrixXd const temp = 0.5 *(phi.block(0, 0, nx + 1, ny + 1) + phi.block(1, 1, nx + 1, ny + 1));


			//txtfile << temp;


			//for (int i = 1; i < nx + 2; i++) {
			//	for (int j = 1; j < ny + 2; j++) {
			//
			//		txtfile << x(i) - 0.5*dx << " " << y(j) - 0.5*dy << " " << 0.25 * (phi(i, j) + phi(i - 1, j) + phi(i - 1, j - 1) + phi(i, j - 1)) << std::endl;

			//	}
			//}

			txtfile.close();


		};
		void pressureExport(std::string const fileName) {


			std::ofstream txtfile;
			txtfile.open(fileName);

			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					txtfile << x(i) << " " << y(j) << " " << p(i, j) << std::endl;

			//Eigen::MatrixXd const temp = 0.5 *(p.block(0, 0, nx + 1, ny + 1) + p.block(1, 1, nx + 1, ny + 1));

			//txtfile << temp;


			txtfile.close();


		};
		void velocityExport(std::string const fileName) {


			std::ofstream txtfile;
			txtfile.open(fileName);


			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					txtfile << x(i) << " " << y(j) << " " << 0.5 * (u(i, j) + u(i - 1, j)) << " " << 0.5 * (v(i, j) + v(i, j - 1)) << std::endl;

	/*		Eigen::MatrixXd const tempU = 0.5 *(u.block(0, 0, nx + 1, ny + 1) + u.block(0, 1, nx + 1, ny + 1));
			Eigen::MatrixXd const tempV = 0.5 *(v.block(0, 0, nx + 1, ny + 1) + v.block(0, 1, nx + 1, ny + 1));

			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					txtfile << x(i) << " " << y(j) << " " << 0.5 * (u(i, j) + u(i - 1, j)) << " " << 0.5 * (v(i, j) + v(i, j - 1)) << std::endl;

	
			txtfile.close();*/


		};
		void vorticityExport(std::string const fileName) {


			std::ofstream txtfile;
			txtfile.open(fileName);

			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny + 1; j++) {

					double const vr = 0.25 * (v(i, j) + v(i, j - 1) + v(i + 1, j - 1) + v(i + 1, j));
					double const vl = 0.25 * (v(i, j) + v(i - 1, j) + v(i - 1, j - 1) + v(i, j - 1));

					double const uu = 0.25 * (u(i, j) + u(i, j + 1) + u(i - 1, j + 1) + u(i - 1, j));
					double const ud = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j - 1) + u(i, j - 1));

					double const dvdx = (vr - vl) / dx;
					double const dudy = (uu - ud) / dy;

					txtfile << x(i) << " " << y(j) << " " << dvdx - dudy << std::endl;
				}
			}
			txtfile.close();


		};

		void meshExport(std::string const fileNameX, std::string const fileNameY) {


			std::ofstream txtfile;
			txtfile.open(fileNameX);

			for (int j = 0; j < ny + 1; j++) {

				for (int i = 0; i < nx + 1; i++) {

					txtfile << x(i) << " ";

				}

				txtfile << std::endl;
			}

			txtfile.close();

			txtfile.open(fileNameY);

			for (int j = 0; j < ny + 1; j++) {

				for (int i = 0; i < nx + 1; i++) {

					txtfile << y(j) << " ";

				}

				txtfile << std::endl;
			}

			txtfile.close();


		};

		void velocityExport_u(std::string const fileName) {


			std::ofstream txtfile;
			txtfile.open(fileName);

			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					txtfile << x(i) << " " << y(j) << " " << 0.5 * (u(i, j) + u(i - 1, j)) << std::endl;

			txtfile.close();


		};
		void velocityExport_v(std::string const fileName) {


			std::ofstream txtfile;
			txtfile.open(fileName);

			for (int i = 1; i < nx + 1; i++)
				for (int j = 1; j < ny + 1; j++)
					txtfile << x(i) << " " << y(j) << " " << 0.5 * (u(i, j) + u(i - 1, j)) << std::endl;

			txtfile.close();


		};


		IncompressibleTwoPhaseFlow(int const numCellX, double const ax, double const bx, int const numCellY, double const ay, double const by) {


			nx = numCellX;
			ny = numCellY;

			NumberOfPhysicalCellX = numCellX;
			NumberOfPhysicalCellY = numCellY;

			NumberOfPointsX = numCellX + 1;
			NumberOfPointsY = numCellY + 1;

			NumberOfCellX = numCellX + 2;
			NumberOfCellY = numCellY + 2;



			u.resize((numCellX + 1), (numCellY + 2));
			v.resize((numCellX + 2), (numCellY + 1));

			u.setZero();
			v.setZero();

			ustar.resize((numCellX + 1), (numCellY + 2));
			vstar.resize((numCellX + 2), (numCellY + 1));

			ustar.setZero();
			vstar.setZero();


			p.resize((numCellX + 2), (numCellY + 2));

			phi.resize((numCellX + 2), (numCellY + 2));
			phiZero.resize((numCellX + 2), (numCellY + 2));
			phiStar.resize((numCellX + 2), (numCellY + 2));
			phiStarStar.resize((numCellX + 2), (numCellY + 2));
			phiPrev.resize((numCellX + 2), (numCellY + 2));

			p.setZero();

			phi.setConstant(-1.0);
			phiZero.setConstant(-1.0);
			phiStar.setConstant(-1.0);
			phiStarStar.setConstant(-1.0);
			phiPrev.setConstant(-1.0);

			fixed.resize((numCellX + 2), (numCellY + 2));

			kappa.resize((numCellX + 2), (numCellY + 2));
			kappa.setZero();


			prevU.resize((numCellX + 1), (numCellY + 2));
			prevV.resize((numCellX + 2), (numCellY + 1));

			prevU.setZero();
			prevV.setZero();


			x.resize(numCellX + 2);
			y.resize(numCellY + 2);

			//x.resize(numCellX + 1);
			//y.resize(numCellY + 1);

			x.setZero();
			y.setZero();


			dx = (bx - ax) / (numCellX);
			dy = (by - ay) / (numCellY);

			alpha = 1.5 * std::max(dx, dy);

			//dtau = (dx*dx) / 2.0;
			//dtau = pow(dx,1.2) / 2.0;
			dtau = std::min(dx, dy) * 1e-1;


			for (int i = 0; i < numCellX + 2; i++)
				x(i) = ax + (i - 0.5)* dx;
			for (int j = 0; j < numCellY + 2; j++)
				y(j) = ay + (j - 0.5)* dy;

			//for (int i = 0; i < numCellX + 1; i++)
			//	x(i) = ax + i * dx;
			//for (int j = 0; j < numCellY + 1; j++)
			//	y(j) = ay + j * dy;


			LevelSetInit();
			PoissonSystemInit();

		};
		~IncompressibleTwoPhaseFlow() {

			clear();

		};

	};






}

