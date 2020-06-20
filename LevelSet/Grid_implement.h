#pragma once


namespace Solver {


	Grid::Grid() {

		NumberOfPoints = 0;

	};
	Grid::~Grid() {

		this->clear();

	};

	int Grid::get_number_of_points() const {

		return NumberOfPoints;

	};

	void Grid::clear() {

		points.clear();
		points.shrink_to_fit();

		NumberOfPoints = 0;

	};

	template<PointMarker mark>
	void Grid::insert_point(Point const & p) {

		points	.push_back(p);
		markers	.push_back(mark);
	
		NumberOfPoints++;

	};

}











/*


namespace Solver {

class Grid {


	private:

		std::vector<Point> points;

		std::vector<Point> PointsInner;
		std::vector<Point> PointsNeumann;
		std::vector<Point> PointsDirichlet;

		int NumberOfPoints;

		int NumberOfPointsInner;
		int NumberOfPointsNeumann;
		int NumberOfPointsDirichlet;


		std::vector<IndexPair> indeces;


	public:

		Grid();
		Grid(int const size);
		~Grid();

		int get_number_of_points() const;
		int get_number_of_points_inner() const;
		int get_number_of_points_neumann() const;
		int get_number_of_points_dirichlet() const;

		void clear();

		template<PointMarker mark>
		void insert_point(int const i, int const j, Point const & p);

	};


	Grid::Grid() {

	};
	Grid::Grid(int const size) {

		int const SqrtSize = (int) floor(sqrt(size));

		points			.reserve(size);
		PointsInner		.reserve(size);
		PointsNeumann	.reserve(SqrtSize);
		PointsDirichlet	.reserve(SqrtSize);

	};
	Grid::~Grid() {

		this->clear();

	};

	int Grid::get_number_of_points() const {

		return NumberOfPoints;

	};
	int Grid::get_number_of_points_inner() const {

		return NumberOfPointsInner;

	};
	int Grid::get_number_of_points_neumann() const {

		return NumberOfPointsNeumann;

	};
	int Grid::get_number_of_points_dirichlet() const {

		return NumberOfPointsDirichlet;

	};

	void Grid::clear() {

		points			.clear();
		PointsInner		.clear();
		PointsNeumann	.clear();
		PointsDirichlet	.clear();


		points			.shrink_to_fit();
		PointsInner		.shrink_to_fit();
		PointsNeumann	.shrink_to_fit();
		PointsDirichlet	.shrink_to_fit();

		NumberOfPoints			= 0;
		NumberOfPointsInner		= 0;
		NumberOfPointsNeumann	= 0;
		NumberOfPointsDirichlet = 0;

	};

	template<PointMarker mark>
	void Grid::insert_point(int const i, int const j, Point const & p) {


		switch (mark) {

		case PointMarker::Inner:

			PointsInner.push_back(p);
			NumberOfPointsInner++;
			break;

		case PointMarker::Neumann:

			PointsNeumann.push_back(p);
			NumberOfPointsNeumann++;
			break;

		case PointMarker::Dirichlet:

			PointsDirichlet.push_back(p);
			NumberOfPointsDirichlet++;
			break;

		}

		indeces	.push_back(IndexPair(i, j));
		points	.push_back(p);

		NumberOfPoints++;

	};

}


*/