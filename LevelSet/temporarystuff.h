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

	ptr.reserve((nx + 0) * (ny + 0) + 1);
	col.reserve((nx + 0) * (ny + 0) * 5);
	val.reserve((nx + 0) * (ny + 0) * 5);

	ptr.push_back(0);


	/************************************************/
	/*                                              */
	/*		Allocate memory for the matrix			*/
	/*												*/
	/************************************************/
	A.resize((nx + 0) * (ny + 0), (nx + 0) * (ny + 0));
	A.reserve((nx + 0) * (ny + 0) * 5);

	/************************************************/
	/*                                              */
	/*		Allocate memory for the right-hand		*/
	/*		side of the system and the solution		*/
	/*												*/
	/************************************************/
	pressureRightHandSide.resize((nx + 0) * (ny + 0));
	pressureSolution.resize((nx + 0) * (ny + 0));


	int diagonalIndex = 0;

	for (int i = 0; i < nx + 0; i++) {
		for (int j = 0; j < ny + 0; j++) {


			/************************************************/
			/*                                              */
			/*		Neumann boundary condition 				*/
			/*			: DivP dot n = 0 					*/
			/*		requires boundary 'ghost' points to		*/
			/*		be zeros								*/
			/*												*/
			/************************************************/
			if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

			}


			/*if (i == 0 && j== 0) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

			}
			else if (i == nx + 1 && j == 0) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

			}
			else if (j == ny + 1 && i == nx + 1) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

			}
			else if (j == ny + 1 && i == 0) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

			}*/


			/*else if (i == 0) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex + (ny + 2));

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (i == nx + 1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex - (ny + 2));

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (j == 0) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex + 1);

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (j == ny + 1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex - 1);

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}*/


			/*else if (i == 1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex + (NumberOfPhysicalCellY + 2));

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (i == (NumberOfPhysicalCellX + 2) - 1-1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex - (NumberOfPhysicalCellY + 2));

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (j == 1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex + 1);

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (j == (NumberOfPhysicalCellY + 2) - 1-1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex - 1);

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}*/

			else {


				/************************************************/
				/*                                              */
				/*		First element in the row				*/
				/*												*/
				/************************************************/
				col.push_back(diagonalIndex - (ny + 0));
				val.push_back(dummyValue);

				/************************************************/
				/*                                              */
				/*		Second element in the row				*/
				/*												*/
				/************************************************/
				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

				/************************************************/
				/*                                              */
				/*		Diagonal element in the row				*/
				/*												*/
				/************************************************/
				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				/************************************************/
				/*                                              */
				/*		Fourth element in the row				*/
				/*												*/
				/************************************************/
				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

				/************************************************/
				/*                                              */
				/*		Fifth element in the row				*/
				/*												*/
				/************************************************/
				col.push_back(diagonalIndex + (ny + 0));
				val.push_back(dummyValue);


			}

			diagonalIndex++;

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
	for (int i = 0; i < (nx + 0) * (ny + 0); i++) {

		int const row_start = ptr[i];
		int const row_end = ptr[i + 1];

		for (int col_index = row_start; col_index < row_end; col_index++) {

			int const j = col[col_index];
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
	val.clear();
	triplet.clear();


	/************************************************/
	/*                                              */
	/*		Assemble matrix in CSR format			*/
	/*												*/
	/************************************************/
	int k = 0;

	for (int i = 1; i <= nx + 0; i++) {
		for (int j = 1; j <= ny + 0; j++) {


			/*if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {

				val.push_back(1.0);

				pressureRightHandSide[k] = 0.0;

			}*/


			if (i == 1) {

				val.push_back(1.0);

				pressureRightHandSide[k] = dt * p(i - 1, j);

			}
			else if (i == nx) {

				val.push_back(1.0);

				pressureRightHandSide[k] = dt * p(i + 1, j);

			}
			else if (j == 1) {

				val.push_back(1.0);

				pressureRightHandSide[k] = dt * p(i, j - 1);

			}
			else if (j == ny) {

				val.push_back(1.0);

				pressureRightHandSide[k] = dt * p(i, j + 1);

			}

			/*if (i == 0 && j == 0) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}
			else if (i == nx + 1 && j == 0) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}
			else if (j == ny + 1 && i == nx + 1) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}
			else if (j == ny + 1 && i == 0) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}*/


			/*else if (i == 0) {

				val.push_back(1.0 / dx);
				val.push_back(-1.0 / dx);

				pressureRightHandSide[k] = 0.0;

			}
			else if (i == nx + 1) {

				val.push_back(1.0 / dx);
				val.push_back(-1.0 / dx);

				pressureRightHandSide[k] = 0.0;

			}
			else if (j == 0) {

				val.push_back(1.0 / dy);
				val.push_back(-1.0 / dy);

				pressureRightHandSide[k] = 0.0;

			}
			else if (j == ny + 1) {

				val.push_back(1.0 / dx);
				val.push_back(-1.0 / dx);

				pressureRightHandSide[k] = 0.0;

			}*/


			/*else if (i == 1) {

				val.push_back(+1.0 / dx);
				val.push_back(-1.0 / dx);

				Rhs[k] = 0.0;

			}
			else if (i == (NumberOfPhysicalCellX + 2) - 1-1) {

				val.push_back(+1.0 / dx);
				val.push_back(-1.0 / dx);

				Rhs[k] = 0.0;

			}
			else if (j == 1) {

				val.push_back(+1.0 / dy);
				val.push_back(-1.0 / dy);

				Rhs[k] = 0.0;

			}
			else if (j == (NumberOfPhysicalCellY + 2) - 1-1) {

				val.push_back(+1.0 / dy);
				val.push_back(-1.0 / dy);

				Rhs[k] = 0.0;

			}*/


			else {

				double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
				double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
				double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
				double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

				double const irho_l = 1.0 / rho(phi_l);
				double const irho_d = 1.0 / rho(phi_d);
				double const irho_u = 1.0 / rho(phi_u);
				double const irho_r = 1.0 / rho(phi_r);


				/************************************************/
				/*                                              */
				/*		First element in the row				*/
				/*												*/
				/************************************************/
				val.push_back(idx2 * irho_l);


				/************************************************/
				/*                                              */
				/*		Second element in the row				*/
				/*												*/
				/************************************************/
				val.push_back(idy2 * irho_d);


				/************************************************/
				/*                                              */
				/*		Diagonal element in the row				*/
				/*												*/
				/************************************************/
				double const val_c = idx2 * irho_l + idx2 * irho_r + idy2 * irho_u + idy2 * irho_d;

				val.push_back(-val_c);

				/************************************************/
				/*                                              */
				/*		Fourth element in the row				*/
				/*												*/
				/************************************************/
				val.push_back(idy2 * irho_u);


				/************************************************/
				/*                                              */
				/*		Fifth element in the row				*/
				/*												*/
				/************************************************/
				val.push_back(idx2 * irho_r);


				double const dudx = (ustar(i, j) - ustar(i - 1, j)) / dx;
				double const dvdy = (vstar(i, j) - vstar(i, j - 1)) / dy;

				pressureRightHandSide[k] = dudx + dvdy;

			}

			k++;

		}
	}


	/************************************************/
	/*                                              */
	/*		Transform CSR format into coordinates	*/
	/*		format which Eigen uses					*/
	/*												*/
	/*		The transformation uses 'slice'			*/
	/*												*/
	/************************************************/
	for (int i = 0; i < (nx + 0)*(ny + 0); i++) {

		int const row_start = ptr[i];
		int const row_end = ptr[i + 1];

		for (int col_index = row_start; col_index < row_end; col_index++) {

			int const j = col[col_index];
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

	//std::cout << A.toDense() << std::endl;


};

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
		
			ptr.reserve((nx + 0) * (ny + 0) + 1);
			col.reserve((nx + 0) * (ny + 0) * 5);
			val.reserve((nx + 0) * (ny + 0) * 5);
		
			ptr.push_back(0);
		
		
			/************************************************/
			/*                                              */
			/*		Allocate memory for the matrix			*/
			/*												*/
			/************************************************/
			A.resize((nx + 0) * (ny + 0), (nx + 0) * (ny + 0));
			A.reserve((nx + 0) * (ny + 0) * 5);
		
			/************************************************/
			/*                                              */
			/*		Allocate memory for the right-hand		*/
			/*		side of the system and the solution		*/
			/*												*/
			/************************************************/
			pressureRightHandSide.resize((nx + 0) * (ny + 0));
			pressureSolution.resize((nx + 0) * (ny + 0));
		
		
			int diagonalIndex = 0;
		
			for (int i = 1; i <= nx + 0; i++) {
				for (int j = 1; j <= ny + 0; j++) {
		
		
					/************************************************/
					/*                                              */
					/*		Neumann boundary condition 				*/
					/*			: DivP dot n = 0 					*/
					/*		requires boundary 'ghost' points to		*/
					/*		be zeros								*/
					/*												*/
					/************************************************/
		
					/*if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {
		
						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + 1);
						col.push_back(diagonalIndex + ny);
		
						val.push_back(dummyValue);
						val.push_back(dummyValue);
						val.push_back(dummyValue);
		
					}*/
					/*else if (i == 0) {
		
						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + (ny + 2));
		
						val.push_back(dummyValue);
						val.push_back(dummyValue);
		
					}
					else if (i == nx + 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == 0) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == ny + 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}*/
					
					/*if (i == 0 && j== 0) {

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);

						col.push_back(diagonalIndex + 1);
						val.push_back(dummyValue);

					}
					else if (i == nx + 1 && j == 0) {

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);

						col.push_back(diagonalIndex + 1);
						val.push_back(dummyValue);

					}
					else if (j == ny + 1 && i == nx + 1) {

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);

						col.push_back(diagonalIndex - 1);
						val.push_back(dummyValue);

					}
					else if (j == ny + 1 && i == 0) {

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);

						col.push_back(diagonalIndex - 1);
						val.push_back(dummyValue);

					}*/
					
					/*else if (i == 0) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (i == nx + 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == 0) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == ny + 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}*/

					/*else if (i == 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + (NumberOfPhysicalCellY + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (i == (NumberOfPhysicalCellX + 2) - 1-1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - (NumberOfPhysicalCellY + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == (NumberOfPhysicalCellY + 2) - 1-1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - 1);
	
						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}*/


					if (i == 1 && j == 1) {

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);

					}
					else {


						if (i != 1) {

							col.push_back(diagonalIndex - (ny + 0));
							val.push_back(dummyValue);

						}
						if (i != nx) {

							col.push_back(diagonalIndex + (ny + 0));
							val.push_back(dummyValue);

						}
						if (j != 1) {

							col.push_back(diagonalIndex - 1);
							val.push_back(dummyValue);

						}
						if (j != ny) {

							col.push_back(diagonalIndex + 1);
							val.push_back(dummyValue);

						}

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);
					
					}


					diagonalIndex++;

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
			for (int i = 0; i < (nx + 0) * (ny + 0); i++) {

				int const row_start = ptr[i];
				int const row_end = ptr[i + 1];

				for (int col_index = row_start; col_index < row_end; col_index++) {

					int const j = col[col_index];
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
			val.clear();
			triplet.clear();


			/************************************************/
			/*                                              */
			/*		Assemble matrix in CSR format			*/
			/*												*/
			/************************************************/
			int k = 0;

			for (int i = 1; i <= nx + 0; i++) {
				for (int j = 1; j <= ny + 0; j++) {


					double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
					double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
					double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
					double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

					double const irho_l = 1.0 / rho(phi_l);
					double const irho_d = 1.0 / rho(phi_d);
					double const irho_u = 1.0 / rho(phi_u);
					double const irho_r = 1.0 / rho(phi_r);


					double value = 0;

					if (i == 1 && j == 1) {

						/************************************************/
						/*                                              */
						/*		Neumann boundary conditions				*/
						/*												*/
						/************************************************/
						val.push_back(1.0);
						pressureRightHandSide[k] = 0.0;

					}
					else {

						if (i != 1) {

							val.push_back(idx2 * irho_l);

							value += idx2 * irho_l;

						}
						if (i != nx) {

							val.push_back(idx2 * irho_r);

							value += idx2 * irho_r;

						}
						if (j != 1) {

							val.push_back(idy2 * irho_d);

							value += idy2 * irho_d;

						}
						if (j != ny) {

							val.push_back(idy2 * irho_u);

							value += idy2 * irho_u;

						}


						//double const val_c = idx2 * irho_l + idx2 * irho_r + idy2 * irho_d + idy2 * irho_u;
						val.push_back(-value);


						double const dudx = (ustar(i, j) - ustar(i - 1, j)) / dx;
						double const dvdy = (vstar(i, j) - vstar(i, j - 1)) / dy;

						pressureRightHandSide[k] = (dudx + dvdy) / dt;

					}

					k++;

				}
			}


			/************************************************/
			/*                                              */
			/*		Transform CSR format into coordinates	*/
			/*		format which Eigen uses					*/
			/*												*/
			/*		The transformation uses 'slice'			*/
			/*												*/
			/************************************************/
			for (int i = 0; i < (nx + 0)*(ny + 0); i++) {

				int const row_start = ptr[i];
				int const row_end = ptr[i + 1];

				for (int col_index = row_start; col_index < row_end; col_index++) {

					int const j = col[col_index];
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

					pt(0, j) = pt(1, j);
					pt(nx + 1, j) = pt(nx, j);

				}
				for (int i = 0; i < nx + 2; i++) {

					pt(i, 0) = pt(i, 1);
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

						pt(i, j) = -(Q + R) / A + (idx2 * irho_l * pt(i - 1, j) + idx2 * irho_r * pt(i + 1, j) +
							idy2 * irho_d * pt(i, j - 1) + idy2 * irho_u * pt(i, j + 1)) / A;

					}
				}

				double Error = 0.0;

				Eigen::MatrixXd const temp = pt - p;

				for (int i = 1; i < nx + 1; i++)
					for (int j = 1; j < ny + 1; j++)
						Error += std::abs(temp(i, j));

				if (Error < DBL_EPSILON) {

					std::cout << "Pressure system error: " << Error << std::endl;
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

			//for (int i = 1; i <= nx + 0; i++)
			//	for (int j = 1; j <= ny + 0; j++)
			//		pressureSolution[(i - 1)*(ny + 0) + j - 1] = p(i, j);


			/************************************************/
			/*                                              */
			/*		Solve the Poisson equation		 		*/
			/*												*/
			/************************************************/
			PoissonSystemAssembly();

			poissonSystemSolverBiCGSTAB.factorize(A);

			pressureSolution = poissonSystemSolverBiCGSTAB.solveWithGuess(pressureRightHandSide, pressureSolution);

			std::cout << "Pressure solution error: " << poissonSystemSolverBiCGSTAB.error() << std::endl;


			/************************************************/
			/*                                              */
			/*		Update Pressures				 		*/
			/*												*/
			/************************************************/
			for (int i = 0; i < nx + 2; i++)
				for (int j = 0; j < ny + 2; j++)
					p(i, j) = pressureSolution[i*(ny + 2) + j];

			//for (int i = 1; i <= nx + 0; i++)
			//	for (int j = 1; j <= ny + 0; j++)
			//		p(i, j) = pressureSolution[(i - 1)*(ny + 0) + j - 1];


		};



Eigen::Vector2d SurfaceTensionForce(int const i, int const j) const {


	// Because Dirac delta function is zero
	//if (abs(phi(i, j)) > alpha)
	//	return Eigen::Vector2d(0.0, 0.0);


	/*double const dphidx = 0.5 * (phi(i + 1, j) - phi(i - 1, j)) / dx;
	double const dphidy = 0.5 * (phi(i, j + 1) - phi(i, j - 1)) / dy;

	double const norm = sqrt(dphidx * dphidx + dphidy * dphidy);

	if (abs(norm) < 1e-10) {

		std::cout << "INVALID NORMAL" << std::endl;

		return Eigen::Vector2d(0.0, 0.0);

	}

	double const normalx = dphidx / norm;
	double const normaly = dphidy / norm;*/


	/************************************************/
	/*                                              */
	/*		1. discretization of curvature Kappa	*/
	/*												*/
	/************************************************/
	/*double const ddphidxx = (phi(i + 1, j) - 2.0 * phi(i, j) + phi(i - 1, j)) / (dx * dx);
	double const ddphidyy = (phi(i, j + 1) - 2.0 * phi(i, j) + phi(i, j - 1)) / (dy * dy);
	double const ddphidyx = 0.25 * ((phi(i + 1, j + 1) - phi(i - 1, j + 1)) - (phi(i + 1, j - 1) - phi(i - 1, j - 1))) / (dx * dy);

	double const curvatureKappa = ((dphidx*dphidx)*ddphidyy - 2.0 * dphidx*dphidy * ddphidyx + (dphidy*dphidy)*ddphidxx) / (norm*norm*norm);*/


	/************************************************/
	/*                                              */
	/*		2. discretization of curvature Kappa	*/
	/*												*/
	/************************************************/
	/*double const dphidx_r = 0.5 * (phi(i + 2, j) - phi(i, j)) / dx;
	double const dphidy_r = 0.5 * (phi(i + 1, j + 1) - phi(i + 1, j - 1)) / dy;

	double const norm_r = sqrt(dphidx_r * dphidx_r + dphidy_r * dphidy_r);


	double const dphidx_l = 0.5 * (phi(i, j) - phi(i - 2, j)) / dx;
	double const dphidy_l = 0.5 * (phi(i - 1, j + 1) - phi(i - 1, j - 1)) / dy;

	double const norm_l = sqrt(dphidx_l * dphidx_l + dphidy_l * dphidy_l);


	double const dphidx_u = 0.5 * (phi(i + 1, j + 1) - phi(i - 1, j + 1)) / dx;
	double const dphidy_u = 0.5 * (phi(i, j + 2) - phi(i, j)) / dy;

	double const norm_u = sqrt(dphidx_u * dphidx_u + dphidy_u * dphidy_u);


	double const dphidx_d = 0.5 * (phi(i + 1, j - 1) - phi(i - 1, j - 1)) / dx;
	double const dphidy_d = 0.5 * (phi(i, j) - phi(i, j - 2)) / dy;

	double const norm_d = sqrt(dphidx_d * dphidx_d + dphidy_d * dphidy_d);


	double const curvatureKappa = 0.5*(dphidx_r / norm_r - dphidx_l / norm_l) / dx + 0.5*(dphidy_u / norm_u - dphidy_d / norm_d);*/

	//return sigma * curvatureKappa * Dirac(phi(i, j)) * Eigen::Vector2d(normalx, normaly);



	double phi_r = 0.5 * (phi(i, j) + phi(i + 1, j));
	double phi_u = 0.5 * (phi(i, j) + phi(i, j + 1));

	double const Fx = 0.5*(DivGradPhi(i, j) * DPhiDx(i, j) + DivGradPhi(i + 1, j) * DPhiDx(i + 1, j));
	double const Fy = 0.5*(DivGradPhi(i, j) * DPhiDy(i, j) + DivGradPhi(i, j + 1) * DPhiDy(i, j + 1));

	double const val_x = std::abs(phi_r) > alpha ? 0.0 : Fx * Dirac(phi_r);
	double const val_y = std::abs(phi_u) > alpha ? 0.0 : Fy * Dirac(phi_u);


	return sigma * Eigen::Vector2d(val_x, val_y);


};
Eigen::Vector2d GravityForce(int const i, int const j) const {

	return Eigen::Vector2d(gx, gy);

};



double UDotDivU_x(int const i, int const j) const {

	/*double const ur = 0.5 * (u(i + 1, j) + u(i, j));
	double const ul = 0.5 * (u(i, j) + u(i - 1, j));

	double const u01 = 0.5 * (u(i, j + 1) + u(i, j));
	double const v01 = 0.5 * (v(i + 1, j) + v(i, j));

	double const u11 = 0.5 * (u(i, j) + u(i, j - 1));
	double const v11 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));


	double const duudx = (ur * ur - ul * ul) / dx;
	double const duvdy = (u01 * v01 - u11 * v11) / dy;


	return duudx + duvdy;*/

	/*double const ur = 0.5 * (u(i + 1, j) * u(i + 1, j) + u(i, j)*u(i, j));
	double const ul = 0.5 * (u(i, j)*u(i, j) + u(i - 1, j)*u(i - 1, j));

	double const u01 = 0.5 * (u(i, j + 1) + u(i, j));
	double const v01 = 0.5 * (v(i + 1, j) + v(i, j));

	double const u11 = 0.5 * (u(i, j) + u(i, j - 1));
	double const v11 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));


	double const duudx = (ur - ul) / dx;
	double const duvdy = (u01 * v01 - u11 * v11) / dy;


	return duudx + duvdy;*/


	/*double const ududx = u(i, j) * (0.5 / dx) * (u(i + 1, j) - u(i - 1, j));
	double const vdudy = 0.25 * ( v(i, j) + v(i, j - 1) + v(i + 1, j) + v(i + 1, j - 1) ) * (0.5 / dy) * (u(i, j + 1) - u(i, j - 1));

	return ududx + vdudy;*/


	double const ududx = u(i, j) > 0.0 ? u(i, j) * (1.0 / dx) * (u(i, j) - u(i - 1, j)) : u(i, j) * (1.0 / dx) * (u(i + 1, j) - u(i, j));
	double const vmean = 0.25 * (v(i, j) + v(i, j - 1) + v(i + 1, j) + v(i + 1, j - 1));
	double const vdudy = vmean > 0.0 ? vmean * (1.0 / dy) * (u(i, j) - u(i, j - 1)) : vmean * (1.0 / dy) * (u(i, j + 1) - u(i, j));

	return ududx + vdudy;

	/*		double const umeanUP = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
			double const umeanDOWN = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));*/
			/*double const ududx = 0.5 * (u(i, j) + u(i - 1, j))*(u(i, j) - u(i - 1, j)) / dx;
			double const vdudy = 0.5 * (v(i, j) + v(i, j - 1))*(u(i, j) - u(i, j - 1)) / dy;

			return ududx + vdudy;*/

};
double UDotDivU_y(int const i, int const j) const {


	/*double const vu = 0.5 * (v(i, j + 1) + v(i, j));
	double const vd = 0.5 * (v(i, j) + v(i, j - 1));

	double const u00 = 0.5 * (u(i - 1, j + 1) + u(i - 1, j));
	double const v00 = 0.5 * (v(i - 1, j) + v(i, j));

	double const u01 = 0.5 * (u(i, j + 1) + u(i, j));
	double const v01 = 0.5 * (v(i + 1, j) + v(i, j));

	double const duvdx = (u01 * v01 - u00 * v00) / dx;
	double const dvvdy = (vu * vu - vd * vd) / dy;


	return duvdx + dvvdy;*/

	/*double const vu = 0.5 * (v(i, j + 1)*v(i, j + 1) + v(i, j)*v(i, j));
	double const vd = 0.5 * (v(i, j)*v(i, j) + v(i, j - 1)*v(i, j - 1));

	double const u00 = 0.5 * (u(i - 1, j + 1) + u(i - 1, j));
	double const v00 = 0.5 * (v(i - 1, j) + v(i, j));

	double const u01 = 0.5 * (u(i, j + 1) + u(i, j));
	double const v01 = 0.5 * (v(i + 1, j) + v(i, j));

	double const duvdx = (u01 * v01 - u00 * v00) / dx;
	double const dvvdy = (vu - vd) / dy;


	return duvdx + dvvdy;*/


	/*double const udvdx = 0.25 * ( u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1) ) * (0.5 / dx) * (v(i + 1, j) - v(i - 1, j));
	double const vdvdy = v(i, j) * (0.5 / dy) * (v(i, j + 1) - v(i, j - 1));

	return udvdx + vdvdy;*/


	double const vdvdy = v(i, j) > 0.0 ? v(i, j) * (1.0 / dy) * (v(i, j) - v(i, j - 1)) : v(i, j) * (1.0 / dy) * (v(i, j + 1) - v(i, j));
	double const umean = 0.25 * (u(i, j) + u(i - 1, j) + u(i - 1, j + 1) + u(i, j + 1));
	double const udvdx = umean > 0.0 ? umean * (1.0 / dx) * (v(i, j) - v(i - 1, j)) : umean * (1.0 / dx) * (v(i + 1, j) - v(i, j));

	return udvdx + vdvdy;

	/*double const udvdx = 0.5 * (u(i, j) + u(i - 1, j))*(v(i, j) - v(i - 1, j)) / dx;
	double const vdvdy = 0.5 * (v(i, j) + v(i, j - 1))*(v(i, j) - v(i, j - 1)) / dy;

	return udvdx + vdvdy;*/

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
void LevelSetInit() {


	/*for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {

			double const X = x(i) - 0.075;
			double const Y = y(j) - 0.075;

			double const r = sqrt(X * X + Y * Y);

			phi(i, j) = r <= 0.025 ? +1.0 : -1.0;

		}
	}*/

	// Square
	/*phi.setConstant(-1.0);
	for (int i = 7; i < 34; i++) {
		for (int j = 7; j < 34; j++) {

			phi(i, j) = 1.0;

		}
	}*/
	// Circle
	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {

			double const X = x(i) - 0.075;
			double const Y = y(j) - 0.05;

			phi(i, j) = 0.025 - sqrt(X * X + Y * Y);

		}
	}

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {

			if (y(j) >= 0.12)
				phi(i, j) = +1.0;

		}
	}

	//for (int i = 0; i < nx + 2; i++) {
	//	for (int j = 0; j < ny + 2; j++) {

	//		if (y(j) <= 0.02)
	//			phi(i, j) = +1.0;

	//	}
	//}

	//LevelSetReInit();

	phiStar = phi;
	phiStarStar = phi;

	/*for (int i = 0; i < (NumberOfPhysicalCellX + 2); i++) {
		for (int j = 0; j < (NumberOfPhysicalCellY + 2); j++) {

			double const X = x(i) - 0.5;
			double const Y = y(j) - 0.5;

			phi(i, j) = 0.2 - sqrt(X*X + Y * Y);

		}
	}*/



};
void LevelSetReInit() {

	phiZero = phi;

	unsigned M = 0;
	double Error = 0.0;

	double tt = 0.0;

	do {

		Error = 0;

		phiStar = phi + dtau * OperatorL(phi);
		phiStarStar = phi + 0.25 * dtau *(OperatorL(phi) + OperatorL(phiStar));

		phi = phi + dtau * (OperatorL(phi) + 4.0 * OperatorL(phiStarStar) + OperatorL(phiStar)) / 6.0;

		for (int i = 0; i < (NumberOfPhysicalCellX + 2); i++)
			for (int j = 0; j < (NumberOfPhysicalCellY + 2); j++)
				Error += abs(phiStar(i, j) - phi(i, j));

		phiStar = phi;

		tt += dtau;

	} while (tt < 2.0*alpha);

	std::cout << Error << std::endl << std::endl;

	return;


	phiStar = phi;
	phiZero = phi;

	unsigned M = 0;
	double Error = 0.0;


	do {

		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny + 1; j++) {


				double const a = (phiStar(i, j) - phiStar(i - 1, j)) / dx;
				double const b = (phiStar(i + 1, j) - phiStar(i, j)) / dx;

				double const c = (phiStar(i, j) - phiStar(i, j - 1)) / dy;
				double const d = (phiStar(i, j + 1) - phiStar(i, j)) / dy;

				//double const sign = Signum(phiStar(i, j));
				double const sign = Signum(i, j);


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

				phi(i, j) = phiStar(i, j) - dtau * sign * G;

			}
		}

		Error = 0;
		M = 0;

		for (int i = 0; i < (NumberOfPhysicalCellX + 2); i++) {
			for (int j = 0; j < (NumberOfPhysicalCellY + 2); j++) {

				if (abs(phi(i, j)) < alpha) {

					M++;
					Error += abs(phiStar(i, j) - phi(i, j));

				}
			}
		}

		if (M == 0)
			Error = INFINITY;
		else
			Error /= M;

		phiStar = phi;

	} while (Error > dtau * dx * dy);


	/*
	do {

		Error = 0;
		//M = 0;

		for (int i = 1; i < (NumberOfPhysicalCellX + 2) - 1; i++) {
			for (int j = 1; j < (NumberOfPhysicalCellY + 2) - 1; j++) {


				double const a = (phiStar(i, j) - phiStar(i - 1, j)) / dx;
				double const b = (phiStar(i + 1, j) - phiStar(i, j)) / dx;

				double const c = (phiStar(i, j) - phiStar(i, j - 1)) / dy;
				double const d = (phiStar(i, j + 1) - phiStar(i, j)) / dy;

				double const sign = Signum(phiStar(i, j));
				double const phi0 = phiZero(i, j);


				double G = 0.0;

				if (phi0 > 0.0) {

					double const ap = std::max(a, 0.0);
					double const bm = std::min(b, 0.0);

					double const cp = std::max(c, 0.0);
					double const dm = std::min(d, 0.0);

					G = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1.0;

				}
				else if (phi0 < 0.0) {

					double const am = std::min(a, 0.0);
					double const bp = std::max(b, 0.0);

					double const cm = std::min(c, 0.0);
					double const dp = std::max(d, 0.0);

					G = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1.0;

				}

				phi(i, j) = phiStar(i, j) - dtau * sign * G;


				Error += std::abs(G);

				if (std::abs(phi(i, j)) < alpha)
					M++;

			}
		}

		phiStar = phi;

	} while (Error > M * dx * dy);
	*/


};
void LevelSetUpdate() {


	for (int i = 1; i < nx + 1; i++)
		for (int j = 1; j < ny + 1; j++)
			phiStar(i, j) = phi(i, j) - dt * UDotDivPhi(i, j);

	phi = phiStar;



	/*phiStar = phi + dt * OperatorF(phi);
	phiStar = phiStar + dt * OperatorF(phiStar);

	phi = 0.5 * (phi + phiStar);*/


	//phiStar = phi + dt * OperatorF(phi);
	//phiStarStar = phiStar + dt * OperatorF(phiStar);

	//phi = 0.5 * (phi + phiStarStar);

};

void LevelSetReInit() {

	phiZero = phi;
	Eigen::MatrixXd phiTemp = phi;


	unsigned M = 0;
	double Error = 0.0;

	double tt = 0.0;


	/*
	do {

		Error = 0;

		//Eigen::MatrixXd const phiTemp1 = OperatorL(phi);
		//phiStar		= phi + dtau * (phiTemp1);
		//Eigen::MatrixXd const phiTemp2 = OperatorL(phiStar);
		//phiStarStar = phi + dtau * (phiTemp1 + phiTemp2) / 4.0;
		//phi			= phi + dtau * (phiTemp1 + 4.0 * OperatorL(phiStarStar) + phiTemp2) / 6.0;



		phiStar		= phi + dtau * (OperatorL(phi));
		phiStarStar = phi + dtau * (OperatorL(phi) + OperatorL(phiStar)) / 4.0;
		phi			= phi + dtau * (OperatorL(phi) + 4.0 * OperatorL(phiStarStar) + OperatorL(phiStar)) / 6.0;

		//phi = phi + dtau * (OperatorL(phi));

		for (int i = 0; i < nx + 2; i++)
			for (int j = 0; j < ny + 2; j++)
				Error += abs(phiTemp(i, j) - phi(i, j)) * dx * dy;

		phiTemp = phi;

		tt += dtau;

	} while (tt < 2.0*alpha);
	//while (Error > 1e-5 * dtau);; (tt < 2.0*alpha)

	std::cout << "Reinitialization error: " << Error << std::endl << std::endl;

	return;*/


	/*phiStar = phi;
	do {

		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny + 1; j++) {


				double const a = (phiTemp(i, j) - phiTemp(i - 1, j)) / dx;
				double const b = (phiTemp(i + 1, j) - phiTemp(i, j)) / dx;

				double const c = (phiTemp(i, j) - phiTemp(i, j - 1)) / dy;
				double const d = (phiTemp(i, j + 1) - phiTemp(i, j)) / dy;

				//double const a = (phiZero(i, j) - phiZero(i - 1, j)) / dx;
				//double const b = (phiZero(i + 1, j) - phiZero(i, j)) / dx;

				//double const c = (phiZero(i, j) - phiZero(i, j - 1)) / dy;
				//double const d = (phiZero(i, j + 1) - phiZero(i, j)) / dy;

				//double const sign = Signum(phiStar(i, j));
				//double const sign = Signum(i, j);
				//double const sign = Signum(i, j, phi);
				//double const sign = Signum(i, j, phiTemp);
				double const sign = Signum(phiZero(i, j));



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

				phiStar(i, j) = phiTemp(i, j) - dtau * sign * G;

			}
		}

		Error = 0;
		M = 0;

		for (int i = 0; i < nx + 2; i++) {
			for (int j = 0; j < ny + 2; j++) {

				if (abs(phiTemp(i, j)) < alpha) {

					M++;
					Error += abs(phiStar(i, j) - phiTemp(i, j));

				}
			}
		}

		//if (M == 0)
		//	Error = INFINITY;
		//else
		//	Error /= M;

		Error /= M;

		phiTemp = phiStar;

	} while (Error > dtau);

	phi = phiStar;
	*/



	Error = 0;
	M = 0;

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {

			if (abs(phiTemp(i, j)) < alpha) {

				M++;
				Error += abs(phi(i, j) - phiTemp(i, j));

			}
		}
	}

	do {

#pragma omp parallel for schedule(static)
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny + 1; j++) {


				double const a = (phiTemp(i, j) - phiTemp(i - 1, j)) / dx;
				double const b = (phiTemp(i + 1, j) - phiTemp(i, j)) / dx;

				double const c = (phiTemp(i, j) - phiTemp(i, j - 1)) / dy;
				double const d = (phiTemp(i, j + 1) - phiTemp(i, j)) / dy;

				double const sign = Signum(phiZero(i, j));
				//double const sign = Signum(i, j, phi); 
				//double const sign = Signum(i, j, phiTemp);


				double G = 0.0;

				if (phiZero(i, j) > 0.0) {

					/*double const ap = std::max(a, 0.0);
					double const bm = -std::min(b, 0.0);

					double const cp = std::max(c, 0.0);
					double const dm = -std::min(d, 0.0);*/

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

				phi(i, j) = phiTemp(i, j) - dtau * sign * G;

			}
		}

		Error = 0;
		M = 0;

		for (int i = 0; i < nx + 2; i++) {
			for (int j = 0; j < ny + 2; j++) {

				if (abs(phiTemp(i, j)) < alpha) {

					M++;
					Error += abs(phi(i, j) - phiTemp(i, j));

				}
			}
		}


		double errorBound = dtau * dx * dy;

		Error /= M;

		phiTemp = phi;



	} while (Error > dtau * dx * dy);



};

void updateBoundaryConditionLevelSet() {


	/*for (int j = 0; j < ny + 2; j++) {

		phi(0, j)	   = 0.0;
		phi(nx + 1, j) = 0.0;

	}
	for (int i = 0; i < nx + 2; i++) {

		phi(i, 0)	   = 0.0;
		phi(i, ny + 1) = 0.0;

	}*/

	for (int j = 0; j < ny + 2; j++) {

		phi(0, j) = phi(1, j);
		phi(nx + 1, j) = phi(nx, j);

	}
	for (int i = 0; i < nx + 2; i++) {

		phi(i, 0) = phi(i, 1);
		phi(i, ny + 1) = phi(i, ny);

	}

};

Eigen::MatrixXd OperatorL(Eigen::MatrixXd const & Phi) {


	Eigen::MatrixXd result(Phi.rows(), Phi.cols());

	result.setZero();


	for (int i = 1; i < nx + 1; i++) {
		for (int j = 1; j < ny + 1; j++) {


			double const a = (phiStar(i, j) - phiStar(i - 1, j)) / dx;
			double const b = (phiStar(i + 1, j) - phiStar(i, j)) / dx;

			double const c = (phiStar(i, j) - phiStar(i, j - 1)) / dy;
			double const d = (phiStar(i, j + 1) - phiStar(i, j)) / dy;

			//double const sign = Signum(phiStar(i, j));
			double const sign = Signum(phiZero(i, j));
			//double const sign = Signum(i, j);


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

			result(i, j) = -sign * G;

		}
	}

	return result;

};


void LevelSetReInit0() {

	phiZero = phi;
	Eigen::MatrixXd phiTemp = phi;


	//initialize();
	//reinitialize();
	//return;



	unsigned M = 0;
	double Error = 0.0;

	double const safetyConst = 1.0;

	/*
	double tt = 0.0;

	do {

		Error = 0;

		//Eigen::MatrixXd const phiTemp1 = OperatorL(phi);
		//phiStar		= phi + dtau * (phiTemp1);
		//Eigen::MatrixXd const phiTemp2 = OperatorL(phiStar);
		//phiStarStar = phi + dtau * (phiTemp1 + phiTemp2) / 4.0;
		//phi			= phi + dtau * (phiTemp1 + 4.0 * OperatorL(phiStarStar) + phiTemp2) / 6.0;



		phiStar = phi + dtau * (OperatorL(phi));
		phiStarStar = phi + dtau * (OperatorL(phi) + OperatorL(phiStar)) / 4.0;
		phi = phi + dtau * (OperatorL(phi) + 4.0 * OperatorL(phiStarStar) + OperatorL(phiStar)) / 6.0;

		//phi = phi + dtau * (OperatorL(phi));

		for (int i = 0; i < nx + 2; i++)
			for (int j = 0; j < ny + 2; j++)
				Error += abs(phiTemp(i, j) - phi(i, j)) * dx * dy;

		phiTemp = phi;

		tt += dtau;

	} while (tt < 10.0*alpha);
	//while (Error > 1e-5 * dtau);; (tt < 2.0*alpha)

	std::cout << "Reinitialization error: " << Error << std::endl << std::endl;

	return;
	*/

	do {

#pragma omp parallel for schedule(static)
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny + 1; j++) {


				double const a = (phiTemp(i, j) - phiTemp(i - 1, j)) / dx;
				double const b = (phiTemp(i + 1, j) - phiTemp(i, j)) / dx;

				double const c = (phiTemp(i, j) - phiTemp(i, j - 1)) / dy;
				double const d = (phiTemp(i, j + 1) - phiTemp(i, j)) / dy;

				double const sign = Signum(phiZero(i, j));
				//double const sign = Signum(i, j, phi); 
				//double const sign = Signum(i, j, phiTemp);
				//double const sign = Signum(i, j, phiZero);


				/*double Gp = 0.0;
				double Gm = 0.0;

				if (sign > 0.0) {

					double const ap = std::max(a, 0.0);
					double const bm = std::min(b, 0.0);

					double const cp = std::max(c, 0.0);
					double const dm = std::min(d, 0.0);

					Gp = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm)) - 1.0;

				}
				else if (sign < 0.0) {

					double const am = std::min(a, 0.0);
					double const bp = std::max(b, 0.0);

					double const cm = std::min(c, 0.0);
					double const dp = std::max(d, 0.0);

					Gm = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp)) - 1.0;

				}

				double const sp = std::max(sign, 0.0);
				double const sm = std::min(sign, 0.0);

				phi(i, j) = phiTemp(i, j) - dtau *( sp * Gp + sm * Gm);*/

				/*double G = 0.0;

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

				phi(i, j) = phiTemp(i, j) - dtau * sign * G;*/

				double G = 0.0;

				/*if (phiZero(i, j) > alpha) {

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

				}*/

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

		if (M == 0)
			Error = INFINITY;
		else
			Error /= M;

		phiTemp = phi;


	} while (Error > safetyConst * dtau * dx * dy);


};

double UpwindPhi_x(double const uij, int const i, int const j) {


	/*if (uij > 0.0)
		return (phi(i, j) - phi(i - 1, j)) / dx;
	if (uij < 0.0)
		return (phi(i + 1, j) - phi(i, j)) / dx;

	return 0.0;*/




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
/*****************************************************************************/
/*                                                                           */
/*		OLD ASSEMBLING OF THE POISSON MATRIX	ALMOST THERE?				 */
/*																			 */
/*****************************************************************************/
/*
void PoissonSystemInit() {


	double const dummyValue = 0.0;


	ptr.clear();
	col.clear();
	val.clear();


	ptr.reserve((nx + 2) * (ny + 2) + 1);
	col.reserve((nx + 2) * (ny + 2) * 5);
	val.reserve((nx + 2) * (ny + 2) * 5);

	ptr.push_back(0);


	A.resize((nx + 2) * (ny + 2), (nx + 2) * (ny + 2));
	A.reserve((nx + 2) * (ny + 2) * 5);

	pressureRightHandSide.resize((nx + 2) * (ny + 2));
	pressureSolution.resize((nx + 2) * (ny + 2));


	int diagonalIndex = 0;

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {



			//if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {

			//	col.push_back(diagonalIndex);
			//	val.push_back(dummyValue);

			//}


			if (i == 0 && j == 0) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

			}
			else if (i == nx + 1 && j == 0) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

			}
			else if (j == ny + 1 && i == nx + 1) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

			}
			else if (j == ny + 1 && i == 0) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

			}


			else if (i == 0) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex + (ny + 2));

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (i == nx + 1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex - (ny + 2));

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (j == 0) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex + 1);

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}
			else if (j == ny + 1) {

				col.push_back(diagonalIndex);
				col.push_back(diagonalIndex - 1);

				val.push_back(dummyValue);
				val.push_back(dummyValue);

			}

			else {


				col.push_back(diagonalIndex - (ny + 2));
				val.push_back(dummyValue);


				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);


				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex + (ny + 2));
				val.push_back(dummyValue);


			}

			diagonalIndex++;

			ptr.push_back(col.size());

		}
	}


	for (int i = 0; i < (nx + 2) * (ny + 2); i++) {

		int const row_start = ptr[i];
		int const row_end = ptr[i + 1];

		for (int col_index = row_start; col_index < row_end; col_index++) {

			int const j = col[col_index];
			double const value = val[col_index];

			triplet.push_back(Eigen::Triplet<double>(i, j, value));

		}
	}



	A.setFromTriplets(triplet.begin(), triplet.end());

	triplet.clear();

	poissonSystemSolverBiCGSTAB.analyzePattern(A);


};
void PoissonSystemAssembly() {


	double const idx2 = 1.0 / (dx * dx);
	double const idy2 = 1.0 / (dy * dy);


	val.clear();
	triplet.clear();



	int k = 0;

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {


			//if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {

			//	val.push_back(1.0);

			//	pressureRightHandSide[k] = 0.0;

			//}


			if (i == 0 && j == 0) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}
			else if (i == nx + 1 && j == 0) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}
			else if (j == ny + 1 && i == nx + 1) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}
			else if (j == ny + 1 && i == 0) {

				val.push_back(1.0);
				val.push_back(-1.0);
				pressureRightHandSide[k] = 0.0;

			}


			else if (i == 0) {

				val.push_back(1.0 / dx);
				val.push_back(-1.0 / dx);

				pressureRightHandSide[k] = 0.0;

			}
			else if (i == nx + 1) {

				val.push_back(1.0 / dx);
				val.push_back(-1.0 / dx);

				pressureRightHandSide[k] = 0.0;

			}
			else if (j == 0) {

				val.push_back(1.0 / dy);
				val.push_back(-1.0 / dy);

				pressureRightHandSide[k] = 0.0;

			}
			else if (j == ny + 1) {

				val.push_back(1.0 / dx);
				val.push_back(-1.0 / dx);

				pressureRightHandSide[k] = 0.0;

			}


			else {

				double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
				double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
				double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
				double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

				double const irho_l = 1.0 / rho(phi_l);
				double const irho_d = 1.0 / rho(phi_d);
				double const irho_u = 1.0 / rho(phi_u);
				double const irho_r = 1.0 / rho(phi_r);

				val.push_back(idx2 * irho_l);



				val.push_back(idy2 * irho_d);



				double const val_c = idx2 * irho_l + idx2 * irho_r + idy2 * irho_u + idy2 * irho_d;

				val.push_back(-val_c);


				val.push_back(idy2 * irho_u);


				val.push_back(idx2 * irho_r);


				double const dudx = (ustar(i, j) - ustar(i - 1, j)) / dx;
				double const dvdy = (vstar(i, j) - vstar(i, j - 1)) / dy;

				pressureRightHandSide[k] = dudx + dvdy;

			}

			k++;

		}
	}


	for (int i = 0; i < (nx + 2)*(ny + 2); i++) {

		int const row_start = ptr[i];
		int const row_end = ptr[i + 1];

		for (int col_index = row_start; col_index < row_end; col_index++) {

			int const j = col[col_index];
			double const value = val[col_index];

			triplet.push_back(Eigen::Triplet<double>(i, j, value));

		}
	}


	A.setFromTriplets(triplet.begin(), triplet.end());

	//std::cout << A.toDense() << std::endl;


};
*/

/*
		void PoissonSystemInit() {


			double const dummyValue = 0.0;


			ptr.clear();
			col.clear();
			val.clear();


			ptr.reserve((nx + 2) * (ny + 2) + 1);
			col.reserve((nx + 2) * (ny + 2) * 5);
			val.reserve((nx + 2) * (ny + 2) * 5);

			ptr.push_back(0);


			A.resize((nx + 2) * (ny + 2), (nx + 2) * (ny + 2));
			A.reserve((nx + 2) * (ny + 2) * 5);

			pressureRightHandSide.resize((nx + 2) * (ny + 2));
			pressureSolution.resize((nx + 2) * (ny + 2));


			int diagonalIndex = 0;

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {


					//if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {
					//
					//	col.push_back(diagonalIndex);
					//	val.push_back(dummyValue);
					//
					//}

					if ((i == 0 && j == 0) || ( i == 0 && j == ny + 1) || (i == nx + 1 && j == 0) || (i == nx + 1 && j == ny + 1)) {

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);

					}
					else if (i == 0) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (i == nx + 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - (ny + 2));

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == 0) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex + 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else if (j == ny + 1) {

						col.push_back(diagonalIndex);
						col.push_back(diagonalIndex - 1);

						val.push_back(dummyValue);
						val.push_back(dummyValue);

					}
					else {


						col.push_back(diagonalIndex - (ny + 2));
						val.push_back(dummyValue);


						col.push_back(diagonalIndex - 1);
						val.push_back(dummyValue);

						col.push_back(diagonalIndex);
						val.push_back(dummyValue);


						col.push_back(diagonalIndex + 1);
						val.push_back(dummyValue);

						col.push_back(diagonalIndex + (ny + 2));
						val.push_back(dummyValue);


					}

					diagonalIndex++;

					ptr.push_back(col.size());

				}
			}


			for (int i = 0; i < (nx + 2) * (ny + 2); i++) {

				int const row_start = ptr[i];
				int const row_end = ptr[i + 1];

				for (int col_index = row_start; col_index < row_end; col_index++) {

					int const j = col[col_index];
					double const value = val[col_index];

					triplet.push_back(Eigen::Triplet<double>(i, j, value));

				}
			}



			A.setFromTriplets(triplet.begin(), triplet.end());

			triplet.clear();

			poissonSystemSolverBiCGSTAB.analyzePattern(A);


		};
		void PoissonSystemAssembly() {


			double const idx2 = 1.0 / (dx * dx);
			double const idy2 = 1.0 / (dy * dy);


			val.clear();
			triplet.clear();



			int k = 0;

			for (int i = 0; i < nx + 2; i++) {
				for (int j = 0; j < ny + 2; j++) {


					if ((i == 0 && j == 0) || (i == 0 && j == ny + 1) || (i == nx + 1 && j == 0) || (i == nx + 1 && j == ny + 1)) {

						val.push_back(1.0);
						pressureRightHandSide[k] = 0.0;

					}
					else if (i == 0) {

						//val.push_back(1.0);
						//val.push_back(1.0);

						val.push_back(1.0 / dx);
						val.push_back(-1.0 / dx);

						pressureRightHandSide[k] = 0.0;

					}
					else if (i == nx + 1) {

						//val.push_back(1.0);
						//val.push_back(1.0);

						val.push_back(1.0 / dx);
						val.push_back(-1.0 / dx);

						pressureRightHandSide[k] = 0.0;

					}
					else if (j == 0) {

						//val.push_back(1.0);
						//val.push_back(1.0);

						val.push_back(1.0 / dy);
						val.push_back(-1.0 / dy);

						pressureRightHandSide[k] = 0.0;

					}
					else if (j == ny + 1) {

						//val.push_back(1.0);
						//val.push_back(1.0);

						val.push_back(1.0 / dy);
						val.push_back(-1.0 / dy);

						pressureRightHandSide[k] = 0.0;

					}
					else {

						double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
						double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
						double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
						double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

						double const irho_l = 1.0 / rho(phi_l);
						double const irho_d = 1.0 / rho(phi_d);
						double const irho_u = 1.0 / rho(phi_u);
						double const irho_r = 1.0 / rho(phi_r);

						val.push_back(idx2 * irho_l);
						val.push_back(idy2 * irho_d);

						double const val_c = idx2 * irho_l + idx2 * irho_r + idy2 * irho_u + idy2 * irho_d;
						val.push_back(-val_c);

						val.push_back(idy2 * irho_u);
						val.push_back(idx2 * irho_r);


						double const dudx = (ustar(i, j) - ustar(i - 1, j)) / dx;
						double const dvdy = (vstar(i, j) - vstar(i, j - 1)) / dy;

						pressureRightHandSide[k] = dudx + dvdy;

					}

					k++;

				}
			}


			for (int i = 0; i < (nx + 2)*(ny + 2); i++) {

				int const row_start = ptr[i];
				int const row_end = ptr[i + 1];

				for (int col_index = row_start; col_index < row_end; col_index++) {

					int const j = col[col_index];
					double const value = val[col_index];

					triplet.push_back(Eigen::Triplet<double>(i, j, value));

				}
			}


			A.setFromTriplets(triplet.begin(), triplet.end());

			//std::cout << A.toDense() << std::endl;


		};
*/

/*****************************************************************************/
/*                                                                           */
/*		OLD ASSEMBLING OF THE POISSON MATRIX								 */
/*																			 */
/*****************************************************************************/
/*
void PoissonSystemInit() {


	double const dummyValue = 0.0;


	ptr.clear();
	col.clear();
	val.clear();




	ptr.reserve((nx + 2) * (ny + 2) + 1);
	col.reserve((nx + 2) * (ny + 2) * 5);
	val.reserve((nx + 2) * (ny + 2) * 5);

	ptr.push_back(0);



	A.resize((nx + 2) * (ny + 2), (nx + 2) * (ny + 2));
	A.reserve((nx + 2) * (ny + 2) * 5);


	pressureRightHandSide.resize((nx + 2) * (ny + 2));
	pressureSolution.resize((nx + 2) * (ny + 2));


	int diagonalIndex = 0;

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {



			if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);

			}
			//else if (i == 1) {

			//	col.push_back(diagonalIndex);
			//	col.push_back(diagonalIndex + (NumberOfPhysicalCellY + 2));

			//	val.push_back(dummyValue);
			//	val.push_back(dummyValue);

			//}
			//else if (i == (NumberOfPhysicalCellX + 2) - 1-1) {

			//	col.push_back(diagonalIndex);
			//	col.push_back(diagonalIndex - (NumberOfPhysicalCellY + 2));

			//	val.push_back(dummyValue);
			//	val.push_back(dummyValue);

			//}
			//else if (j == 1) {

			//	col.push_back(diagonalIndex);
			//	col.push_back(diagonalIndex + 1);

			//	val.push_back(dummyValue);
			//	val.push_back(dummyValue);

			//}
			//else if (j == (NumberOfPhysicalCellY + 2) - 1-1) {

			//	col.push_back(diagonalIndex);
			//	col.push_back(diagonalIndex - 1);

			//	val.push_back(dummyValue);
			//	val.push_back(dummyValue);

			//}
			else {



				col.push_back(diagonalIndex - (ny + 2));
				val.push_back(dummyValue);


				col.push_back(diagonalIndex - 1);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex);
				val.push_back(dummyValue);


				col.push_back(diagonalIndex + 1);
				val.push_back(dummyValue);

				col.push_back(diagonalIndex + (ny + 2));
				val.push_back(dummyValue);


			}

			diagonalIndex++;

			ptr.push_back(col.size());

		}
	}

	for (int i = 0; i < (nx + 2) * (ny + 2); i++) {

		int const row_start = ptr[i];
		int const row_end = ptr[i + 1];

		for (int col_index = row_start; col_index < row_end; col_index++) {

			int const j = col[col_index];
			double const value = val[col_index];

			triplet.push_back(Eigen::Triplet<double>(i, j, value));

		}
	}



	A.setFromTriplets(triplet.begin(), triplet.end());

	triplet.clear();


	poissonSystemSolverBiCGSTAB.analyzePattern(A);


};
void PoissonSystemAssembly() {


	double const idx2 = 1.0 / (dx * dx);
	double const idy2 = 1.0 / (dy * dy);


	val.clear();
	triplet.clear();



	int k = 0;

	for (int i = 0; i < nx + 2; i++) {
		for (int j = 0; j < ny + 2; j++) {

			if (i == 0 || i == nx + 1 || j == 0 || j == ny + 1) {

				val.push_back(1.0);

				pressureRightHandSide[k] = 0.0;

			}

			//else if (i == 1) {

			//	val.push_back(+1.0 / dx);
			//	val.push_back(-1.0 / dx);

			//	Rhs[k] = 0.0;

			//}
			//else if (i == (NumberOfPhysicalCellX + 2) - 1-1) {

			//	val.push_back(+1.0 / dx);
			//	val.push_back(-1.0 / dx);

			//	Rhs[k] = 0.0;

			//}
			//else if (j == 1) {

			//	val.push_back(+1.0 / dy);
			//	val.push_back(-1.0 / dy);

			//	Rhs[k] = 0.0;

			//}
			//else if (j == (NumberOfPhysicalCellY + 2) - 1-1) {

			//	val.push_back(+1.0 / dy);
			//	val.push_back(-1.0 / dy);

			//	Rhs[k] = 0.0;

			//}
			else {

				double const phi_l = 0.5*(phi(i, j) + phi(i - 1, j));
				double const phi_d = 0.5*(phi(i, j) + phi(i, j - 1));
				double const phi_u = 0.5*(phi(i, j) + phi(i, j + 1));
				double const phi_r = 0.5*(phi(i, j) + phi(i + 1, j));

				double const irho_l = 1.0 / rho(phi_l);
				double const irho_d = 1.0 / rho(phi_d);
				double const irho_u = 1.0 / rho(phi_u);
				double const irho_r = 1.0 / rho(phi_r);


				val.push_back(idx2 * irho_l);



				val.push_back(idy2 * irho_d);



				double const val_c = idx2 * irho_l + idx2 * irho_r + idy2 * irho_u + idy2 * irho_d;

				val.push_back(-val_c);


				val.push_back(idy2 * irho_u);


				val.push_back(idx2 * irho_r);


				double const dudx = (ustar(i, j) - ustar(i - 1, j)) / dx;
				double const dvdy = (vstar(i, j) - vstar(i, j - 1)) / dy;

				pressureRightHandSide[k] = dudx + dvdy;

			}

			k++;

		}
	}



	for (int i = 0; i < (nx + 2)*(ny + 2); i++) {

		int const row_start = ptr[i];
		int const row_end = ptr[i + 1];

		for (int col_index = row_start; col_index < row_end; col_index++) {

			int const j = col[col_index];
			double const value = val[col_index];

			triplet.push_back(Eigen::Triplet<double>(i, j, value));

		}
	}



	A.setFromTriplets(triplet.begin(), triplet.end());


};
*/


/*****************************************************************************/
/*                                                                           */
/*		5-TH ORDER WENO UPWIND SCHEME FOR THE LEVEL SET ADVECTION			 */
/*																			 */
/*****************************************************************************/
/*
double UDotDivPhi(int const i, int const j) {



			double uij;
			double vij;

			//uij = 0.5 * (u(i, j) + u(i - 1, j));
			//vij = 0.5 * (v(i, j) + v(i, j - 1));

			if (i == 1)
				uij = (5.0 * u(i - 1, j) + 15.0* u(i, j) - 5.0*u(i + 1, j) + u(i + 2, j)) / 16.0;
			else if (i == (NumberOfPhysicalCellX + 1) - 1)
				uij = (u(i - 3, j) - 5.0* u(i - 2, j) + 15.0*u(i - 1, j) + 5.0*u(i, j)) / 16.0;
			else
				uij = -(u(i - 2, j) + 9.0* u(i - 1, j) + 9.0*u(i, j) - u(i + 1, j)) / 16.0;


			if (j == 1)
				vij = (5.0 * v(i, j - 1) + 15.0* v(i, j) - 5.0*v(i, j + 1) + v(i, j + 2)) / 16.0;
			else if (j == (NumberOfPhysicalCellY + 1) - 1)
				vij = (v(i, j - 3) - 5.0* v(i, j - 2) + 15.0*v(i, j - 1) + 5.0*v(i, j)) / 16.0;
			else
				vij = -(v(i, j - 2) + 9.0* v(i, j - 1) + 9.0*v(i, j) - v(i, j + 1)) / 16.0;


			double const dphidx = UpwindPhi_x(uij, i, j);
			double const dphidy = UpwindPhi_y(vij, i, j);

			return uij * dphidx + vij * dphidy;

		};

double UpwindPhi_x(double const uij, int const i, int const j) {


	//if (uij > 0.0)
	//	return (phi(i, j) - phi(i - 1, j)) / dx;
	//if (uij < 0.0)
	//	return (phi(i + 1, j) - phi(i, j)) / dx;
	//
	//return 0.0;

	if (uij > 0.0) {

		if (i == 1)
			return q3minus_x(i, j) / 3.0 + 5.0 * q4minus_x(i, j) / 6.0 - q5minus_x(i, j) / 6.0;
		else if (i == 2 || i == (NumberOfPhysicalCellX + 1) - 2)
			return -q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0;
		else if (i == (NumberOfPhysicalCellX + 1) - 1)
			return q1minus_x(i, j) / 3.0 - 7.0 * q2minus_x(i, j) / 6.0 + 11.0 * q3minus_x(i, j) / 6.0;


		double const Value1 = omega1minus_x(i, j)*(q1minus_x(i, j) / 3.0 - 7.0 * q2minus_x(i, j) / 6.0 + 11.0 * q3minus_x(i, j) / 6.0);
		double const Value2 = omega2minus_x(i, j)*(-q2minus_x(i, j) / 6.0 + 5.0 * q3minus_x(i, j) / 6.0 + q4minus_x(i, j) / 3.0);
		double const Value3 = omega3minus_x(i, j)*(q3minus_x(i, j) / 3.0 + 5.0 * q4minus_x(i, j) / 6.0 - q5minus_x(i, j) / 6.0);

		return Value1 + Value2 + Value3;

	}
	else if (uij < 0.0) {

		if (i == 1)
			return q1plus_x(i, j) / 3.0 - 7.0 * q2plus_x(i, j) / 6.0 + 11.0 * q3plus_x(i, j) / 6.0;
		else if (i == 2 || i == (NumberOfPhysicalCellX + 1) - 2)
			return -q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0;
		else if (i == (NumberOfPhysicalCellX + 1) - 1)
			return q3plus_x(i, j) / 3.0 + 5.0 * q4plus_x(i, j) / 6.0 - q5plus_x(i, j) / 6.0;

		double const Value1 = omega1plus_x(i, j)*(q1plus_x(i, j) / 3.0 - 7.0 * q2plus_x(i, j) / 6.0 + 11.0 * q3plus_x(i, j) / 6.0);
		double const Value2 = omega2plus_x(i, j)*(-q2plus_x(i, j) / 6.0 + 5.0 * q3plus_x(i, j) / 6.0 + q4plus_x(i, j) / 3.0);
		double const Value3 = omega3plus_x(i, j)*(q3plus_x(i, j) / 3.0 + 5.0 * q4plus_x(i, j) / 6.0 - q5plus_x(i, j) / 6.0);

		return Value1 + Value2 + Value3;

	}

	return 0.0;

};
double UpwindPhi_y(double const vij, int const i, int const j) {

	//if (vij > 0.0)
	//	return (phi(i, j) - phi(i, j - 1)) / dy;
	//if (vij < 0.0)
	//	return (phi(i, j + 1) - phi(i, j)) / dy;
	//
	//return 0.0;

	if (vij > 0.0) {

		if (j == 1)
			return q3minus_y(i, j) / 3.0 + 5.0 * q4minus_y(i, j) / 6.0 - q5minus_y(i, j) / 6.0;
		else if (j == 2 || j == (NumberOfPhysicalCellY + 1) - 2)
			return -q2minus_y(i, j) / 6.0 + 5.0 * q3minus_y(i, j) / 6.0 + q4minus_y(i, j) / 3.0;
		else if (j == (NumberOfPhysicalCellY + 1) - 1)
			return q1minus_y(i, j) / 3.0 - 7.0 * q2minus_y(i, j) / 6.0 + 11.0 * q3minus_y(i, j) / 6.0;

		double const Value1 = omega1minus_y(i, j)*(q1minus_y(i, j) / 3.0 - 7.0 * q2minus_y(i, j) / 6.0 + 11.0 * q3minus_y(i, j) / 6.0);
		double const Value2 = omega2minus_y(i, j)*(-q2minus_y(i, j) / 6.0 + 5.0 * q3minus_y(i, j) / 6.0 + q4minus_y(i, j) / 3.0);
		double const Value3 = omega3minus_y(i, j)*(q3minus_y(i, j) / 3.0 + 5.0 * q4minus_y(i, j) / 6.0 - q5minus_y(i, j) / 6.0);

		return Value1 + Value2 + Value3;

	}
	else if (vij < 0.0) {

		if (j == 1)
			return q1plus_y(i, j) / 3.0 - 7.0 * q2plus_y(i, j) / 6.0 + 11.0 * q3plus_y(i, j) / 6.0;
		else if (j == 2 || j == (NumberOfPhysicalCellY + 1) - 2)
			return -q2plus_y(i, j) / 6.0 + 5.0 * q3plus_y(i, j) / 6.0 + q4plus_y(i, j) / 3.0;
		else if (j == (NumberOfPhysicalCellY + 1) - 1)
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
*/