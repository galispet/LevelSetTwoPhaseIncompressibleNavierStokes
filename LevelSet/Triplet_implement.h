#pragma once


namespace Solver {

	IndexPair::IndexPair() {

		indeces[0] = 0;
		indeces[1] = 0;

	};
	IndexPair::IndexPair(int const i, int const j) {

		indeces[0] = i;
		indeces[1] = j;

	};
	IndexPair::~IndexPair() {

	};
	   
}