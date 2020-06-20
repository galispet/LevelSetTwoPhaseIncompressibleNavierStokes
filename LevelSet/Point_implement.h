#pragma once


namespace Solver {

	
		Point::Point() {
		
			coordinates[0] = 0.0;
			coordinates[1] = 0.0;

		};
		Point::Point(double const x, double const y) {

			coordinates[0] = x;
			coordinates[1] = y;

		};
		Point::Point(Point const & p) {

			coordinates[0] = p.coordinates[0];
			coordinates[1] = p.coordinates[1];

		};
		Point::~Point() {


		};


		double Point::get_x() const {

			return this->coordinates[0];

		};
		double Point::get_y() const {

			return this->coordinates[1];

		};

		void Point::set_x(double const x) {

			this->coordinates[0] = x;
		
		};
		void Point::set_y(double const y) {
		
			this->coordinates[1] = y;

		};

		double & Point::set_x() {

			return this->coordinates[0];
		
		};
		double & Point::set_y() {

			return this->coordinates[1];
		
		};


}