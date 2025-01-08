#include "main.h"


#include <Eigen/Dense>
using namespace Eigen;
//
//// Function to fit an ellipse using Gauss method
//VectorXd fitellipse(const std::vector<double>& x, const std::vector<double>& y) {
//	size_t n = x.size();
//	if (n < 5) {
//		std::cerr << "Need at least 5 points to fit an ellipse." << std::endl;
//		return VectorXd(0);
//	}
//
//	// Step 1: Construct the design matrix D
//	MatrixXd D(n, 6);
//	for (size_t i = 0; i < n; ++i) {
//		D(i, 0) = x[i] * x[i];
//		D(i, 1) = x[i] * y[i];
//		D(i, 2) = y[i] * y[i];
//		D(i, 3) = x[i];
//		D(i, 4) = y[i];
//		D(i, 5) = 1.0;
//	}
//
//	// Step 2: Compute the scatter matrix S = D^T * D
//	MatrixXd S = D.transpose() * D;
//
//	// Step 3: Define the constraint matrix C
//	MatrixXd C(6, 6);
//	C << 0, 0, 2, 0, 0, 0,
//		0, -1, 0, 0, 0, 0,
//		2, 0, 0, 0, 0, 0,
//		0, 0, 0, 0, 0, 0,
//		0, 0, 0, 0, 0, 0,
//		0, 0, 0, 0, 0, 0;
//
//	// Step 4: Solve the generalized eigenvalue problem S * v = ? * C * v
//	EigenSolver<MatrixXd> solver(S);
//	MatrixXd eigenvalues = solver.eigenvalues().real();
//	MatrixXd eigenvectors = solver.eigenvectors().real();
//
//	// Find the solution satisfying the ellipse constraint
//	int bestIndex = -1;
//	for (int i = 0; i < eigenvalues.size(); ++i) {
//		VectorXd v = eigenvectors.col(i);
//		if (v.transpose() * C * v > 0) {
//			bestIndex = i;
//			break;
//		}
//	}
//
//	if (bestIndex == -1) {
//		std::cerr << "No valid ellipse found." << std::endl;
//		return VectorXd(0);
//	}
//
//	return eigenvectors.col(bestIndex).transpose();
//}
//
//
//
//class GaussellipseFit {
//public:
//	static VectorXd fitellipse(const std::vector<double>& x, const std::vector<double>& y) {
//		if (x.size() != y.size() || x.size() < 5) {
//			throw std::runtime_error("Need at least 5 points to fit an ellipse");
//		}
//
//		// Build design matrix D
//		Eigen::MatrixXd D(x.size(), 6);
//		for (size_t i = 0; i < x.size(); i++) {
//			D(i, 0) = x[i] * x[i];
//			D(i, 1) = x[i] * y[i];
//			D(i, 2) = y[i] * y[i];
//			D(i, 3) = x[i];
//			D(i, 4) = y[i];
//			D(i, 5) = 1.0;
//		}
//
//		// Build scatter matrix S
//		Eigen::MatrixXd S = D.transpose() * D;
//
//		// Build constraint matrix C
//		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6);
//		C(0, 2) = 2;
//		C(2, 0) = 2;
//		C(1, 1) = -1;
//
//		// Solve generalized eigensystem
//		Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(S, C);
//
//		// Find the positive eigenvalue
//		Eigen::VectorXd eigenvalues = ges.eigenvalues().real();
//		int positiveIndex = -1;
//		for (int i = 0; i < eigenvalues.size(); i++) {
//			if (eigenvalues(i) > 0) {
//				positiveIndex = i;
//				break;
//			}
//		}
//
//		if (positiveIndex == -1) {
//			throw std::runtime_error("No valid ellipse found");
//		}
//
//		// Get the coefficients
//		Eigen::VectorXd coefficients = ges.eigenvectors().real().col(positiveIndex);
//
//		// Extract ellipse parameters
//		return coefficients;// extractellipseParameters(coefficients);
//	}
//
//	ellipseParameters extractellipseParameters(const Eigen::VectorXd& coefficients)
//	{
//		double a = coefficients(0);
//		double b = coefficients(1) / 2;
//		double c = coefficients(2);
//		double d = coefficients(3);
//		double e = coefficients(4);
//		double f = coefficients(5);
//
//		// Calculate center
//		double centerX = (2 * c * d - b * e) / (b * b - 4 * a * c);
//		double centerY = (2 * a * e - b * d) / (b * b - 4 * a * c);
//
//		// Calculate rotation angle
//		double theta = 0.5 * atan2(b, (a - c));
//
//		// Calculate semi-axes
//		double ct = cos(theta);
//		double st = sin(theta);
//		double ct2 = ct * ct;
//		double st2 = st * st;
//		double a2 = a * ct2 + b * ct * st + c * st2;
//		double c2 = a * st2 - b * ct * st + c * ct2;
//
//		// Calculate constants
//		double term = 2 * (a * centerX * centerX + b * centerX * centerY +
//			c * centerY * centerY + d * centerX + e * centerY + f);
//
//		double semiMajor = sqrt(abs(term / (2 * std::min(a2, c2))));
//		double semiMinor = sqrt(abs(term / (2 * std::max(a2, c2))));
//
//		if (a2 > c2) {
//			std::swap(semiMajor, semiMinor);
//			theta += pi / 2;
//		}
//
//		ellipseParameters params;
//		params.centerX = centerX;
//		params.centerY = centerY;
//		params.semiMajor = semiMajor;
//		params.semiMinor = semiMinor;
//
//		return params;
//	}
//};
//
//
//struct ellipseFoci {
//	double focus1X;
//	double focus1Y;
//	double focus2X;
//	double focus2Y;
//};
//
//
//
//
//ellipseFoci calculateFoci(const ellipseParameters& params) {
//	// c² = a² - b²
//	double c = sqrt(params.semiMajor * params.semiMajor -
//		params.semiMinor * params.semiMinor);
//
//	// Calculate foci before rotation
//	double fx = c;
//	double fy = 0;
//
//	// Rotate foci
//	double fx1 = fx * cos(0) - fy * sin(0);
//	double fy1 = fx * sin(0) + fy * cos(0);
//
//	// Return both foci (translated to ellipse center)
//	ellipseFoci foci;
//	foci.focus1X = params.centerX + fx1;
//	foci.focus1Y = params.centerY + fy1;
//	foci.focus2X = params.centerX - fx1;
//	foci.focus2Y = params.centerY - fy1;
//
//	return foci;
//}
//
//
void DrawEllipse(double cx, double cy, double rx, double ry, int num_segments)
{
	double theta = 2 * pi / double(num_segments);
	double c = cos(theta);
	double s = sin(theta);
	double t;

	double x = 1;//we start at angle = 0 
	double y = 0;

	glBegin(GL_LINE_LOOP);
	for (int ii = 0; ii < num_segments; ii++)
	{
		//apply radius and offset
		glVertex2d(x * rx + cx, y * ry + cy);//output vertex 

		//apply the rotation matrix
		t = x;
		x = c * x - s * y;
		y = s * t + c * y;
	}
	glEnd();
}




int main(int argc, char** argv)
{
	cout << setprecision(20) << endl;


	//cout << sqrt((grav_constant * 1.988435e30 / 69817079000.0) * (1 - 0.20563069) / (1 + 0.20563069)) << endl;







//	return 0;




	glutInit(&argc, argv);
	init_opengl(win_x, win_y);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
	glutKeyboardFunc(keyboard_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutPassiveMotionFunc(passive_motion_func);
	//glutIgnoreKeyRepeat(1);
	glutMainLoop();
	glutDestroyWindow(win_id);

	return 0;
}

custom_math::vector_3 grav_acceleration(const custom_math::vector_3& pos, const custom_math::vector_3& vel, const double G)
{
	custom_math::vector_3 grav_dir = sun_pos - pos;

	double distance = grav_dir.length();

	grav_dir.normalize();
	custom_math::vector_3 accel = grav_dir * (G * sun_mass / pow(distance, 2.0));

	return accel;
}



void proceed_Euler(custom_math::vector_3& pos, custom_math::vector_3& vel, const double G, const double dt)
{
	custom_math::vector_3 accel = grav_acceleration(pos, vel, G);

	vel += accel * dt;
	pos += vel * dt;
}

void proceed_RK4(custom_math::vector_3& pos, custom_math::vector_3& vel, const double G, const double dt)
{
	static const double one_sixth = 1.0 / 6.0;

	custom_math::vector_3 k1_velocity = vel;
	custom_math::vector_3 k1_acceleration = grav_acceleration(pos, k1_velocity, G);
	custom_math::vector_3 k2_velocity = vel + k1_acceleration * dt * 0.5;
	custom_math::vector_3 k2_acceleration = grav_acceleration(pos + k1_velocity * dt * 0.5, k2_velocity, G);
	custom_math::vector_3 k3_velocity = vel + k2_acceleration * dt * 0.5;
	custom_math::vector_3 k3_acceleration = grav_acceleration(pos + k2_velocity * dt * 0.5, k3_velocity, G);
	custom_math::vector_3 k4_velocity = vel + k3_acceleration * dt;
	custom_math::vector_3 k4_acceleration = grav_acceleration(pos + k3_velocity * dt, k4_velocity, G);

	vel += (k1_acceleration + (k2_acceleration + k3_acceleration) * 2.0 + k4_acceleration) * one_sixth * dt;
	pos += (k1_velocity + (k2_velocity + k3_velocity) * 2.0 + k4_velocity) * one_sixth * dt;
}


void proceed_symplectic4(custom_math::vector_3& pos, custom_math::vector_3& vel, const double G, const double dt)
{
	static double const cr2 = pow(2.0, 1.0 / 3.0);

	static const double c[4] =
	{
		1.0 / (2.0 * (2.0 - cr2)),
		(1.0 - cr2) / (2.0 * (2.0 - cr2)),
		(1.0 - cr2) / (2.0 * (2.0 - cr2)),
		1.0 / (2.0 * (2.0 - cr2))
	};

	static const double d[4] =
	{
		1.0 / (2.0 - cr2),
		-cr2 / (2.0 - cr2),
		1.0 / (2.0 - cr2),
		0.0
	};

	pos += vel * c[0] * dt;
	vel += grav_acceleration(pos, vel, G) * d[0] * dt;

	pos += vel * c[1] * dt;
	vel += grav_acceleration(pos, vel, G) * d[1] * dt;

	pos += vel * c[2] * dt;
	vel += grav_acceleration(pos, vel, G) * d[2] * dt;

	pos += vel * c[3] * dt;
	//	vel += grav_acceleration(pos, vel, G) * d[3] * dt; // last element d[3] is always 0
}


#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Assuming you have a library like Eigen for matrix operations or implement your own.

struct Observation {
	double time;  // Time of observation
	double azimuth; // In radians
	double elevation; // In radians
};

struct Orbit {
	double a; // Semi-major axis
	// Other orbital parameters if needed (e, i, etc.)
};

class OrbitEstimator {
private:
	std::vector<Observation> observations;

	// Convert azimuth and elevation to Cartesian coordinates
	Eigen::Vector3d azElToCartesian(double az, double el, double r = 1.0) {
		return Eigen::Vector3d(
			r * cos(az) * cos(el),
			r * sin(az) * cos(el),
			r * sin(el)
		);
	}

	// Fit an ellipse to a set of points (simplified)
	Orbit fitEllipse(const std::vector<Eigen::Vector3d>& points) {
		// This is a placeholder; actual implementation would use something like:
		// - Least Squares Ellipse Fitting (LSEF)
		// - Direct least square fitting of ellipses
		// Here we assume we can directly compute 'a' from points.
		double a = 0; // Simplified calculation, real method would involve more complex math
		return Orbit{ a };
	}

public:
	void addObservation(const Observation& obs) {
		observations.push_back(obs);
	}

	Orbit estimateOrbit() {
		std::vector<Eigen::Vector3d> points;
		for (const auto& obs : observations) {
			points.push_back(azElToCartesian(obs.azimuth, obs.elevation));
		}

		// Here we assume we can fit an ellipse in 3D space, which is complex. 
		// For simplicity, we might project onto the plane of the orbit if known.
		Orbit result = fitEllipse(points);
		return result;
	}
};



struct Point {
	double x, y;
};

class EllipseParameters
{
public:
	double centerX = 0;
	double centerY = 0;
	double semiMajor = 0;
	double semiMinor = 0;
	double angle = 0;
};


EllipseParameters extractEllipseParameters(const Eigen::VectorXd& coefficients)
{
	double a = coefficients(0);
	double b = coefficients(1) / 2;
	double c = coefficients(2);
	double d = coefficients(3);
	double e = coefficients(4);
	double f = 1;

	// Calculate center
	double centerX = (2 * c * d - b * e) / (b * b - 4 * a * c);
	double centerY = (2 * a * e - b * d) / (b * b - 4 * a * c);

	// Calculate rotation angle
	double theta = 0.5 * atan2(b, (a - c));

	// Calculate semi-axes
	double ct = cos(theta);
	double st = sin(theta);
	double ct2 = ct * ct;
	double st2 = st * st;
	double a2 = a * ct2 + b * ct * st + c * st2;
	double c2 = a * st2 - b * ct * st + c * ct2;

	// Calculate constants
	double term = 2 * (a * centerX * centerX + b * centerX * centerY +
		c * centerY * centerY + d * centerX + e * centerY + f);

	double semiMajor = sqrt(abs(term / (2 * std::min(a2, c2))));
	double semiMinor = sqrt(abs(term / (2 * std::max(a2, c2))));

	if (a2 > c2) {
		std::swap(semiMajor, semiMinor);
		theta += pi / 2;
	}

	EllipseParameters params;
	params.centerX = centerX;
	params.centerY = centerY;
	params.semiMajor = semiMajor;
	params.semiMinor = semiMinor;
	params.angle = theta;

	return params;
}









EllipseParameters fitEllipse(const std::vector<Point>& points) 
{
	const size_t num_points = 5;

	if (points.size() != num_points) {
		std::cerr << "Error: Exactly " << num_points << " points are required.\n";
		return EllipseParameters();
	}

	Eigen::MatrixXd A(num_points, 6);
	Eigen::VectorXd b(num_points);

	// Fill the matrix A and vector b with the equations from the points
	for (size_t i = 0; i < num_points; ++i)
	{
		double x = points[i].x;
		double y = points[i].y;
		A(i, 0) = x * x;       // Coefficient for x^2
		A(i, 1) = x * y;       // Coefficient for xy
		A(i, 2) = y * y;       // Coefficient for y^2
		A(i, 3) = x;           // Coefficient for x
		A(i, 4) = y;           // Coefficient for y
		A(i, 5) = 1;           // Constant term
		b(i) = -1;             // Right-hand side is -1. This is important!
	}

	//cout << A << endl;

	//cout << b << endl;

	// Solve for the ellipse parameters
	Eigen::VectorXd ellipseParams = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

	//// Extract parameters
	//double A_ = ellipseParams(0);
	//double B_ = ellipseParams(1);
	//double C_ = ellipseParams(2);
	//double D_ = ellipseParams(3);
	//double E_ = ellipseParams(4);
	//double F_ = ellipseParams(5);

	//cout << A_ << ' ' << B_ << ' ' << C_ << ' ' << D_ << ' ' << E_ << ' ' << F_ << endl;

	return extractEllipseParameters(ellipseParams);
}



#include <cmath>
#include <array>
//
//class OrbitalParameters {
//private:
//	double semiMajorAxis;    // a: semi-major axis in kilometers
//	double eccentricity;     // e: orbital eccentricity (dimensionless)
//	double inclination;      // i: orbital inclination in radians
//	double raan;            // Ω: right ascension of ascending node in radians
//	double argPerigee;      // ω: argument of perigee in radians
//	double meanAnomaly;     // M: mean anomaly in radians
//	double mu;              // μ: standard gravitational parameter (GM) in km³/s²
//
//public:
//	OrbitalParameters(double a, double e, double i, double omega, double w, double M, double gravitationalParam = 398600.4418)
//		: semiMajorAxis(a), eccentricity(e), inclination(i),
//		raan(omega), argPerigee(w), meanAnomaly(M), mu(gravitationalParam) {}
//
//	// Solve Kepler's equation using Newton-Raphson method
//	double solveKepler(double M, double e, double tolerance = 1e-8) const {
//		double E = M; // Initial guess
//		double delta;
//
//		do {
//			delta = (E - e * sin(E) - M) / (1 - e * cos(E));
//			E -= delta;
//		} while (std::abs(delta) > tolerance);
//
//		return E;
//	}
//
//	// Get position vector in orbital plane
//	std::array<double, 2> getPositionInOrbitalPlane(double trueAnomaly) const {
//		double r = semiMajorAxis * (1 - eccentricity * eccentricity) /
//			(1 + eccentricity * cos(trueAnomaly));
//
//		return { r * cos(trueAnomaly), r * sin(trueAnomaly) };
//	}
//
//	// Convert orbital position to cartesian coordinates (ECI frame)
//	std::array<double, 3> getCartesianPosition(double time) const {
//		// Solve Kepler's equation to get eccentric anomaly
//		double E = solveKepler(meanAnomaly +
//			sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis)) * time,
//			eccentricity);
//
//		// Calculate true anomaly
//		double trueAnomaly = 2 * atan(sqrt((1 + eccentricity) / (1 - eccentricity)) *
//			tan(E / 2));
//
//		// Get position in orbital plane
//		auto pos = getPositionInOrbitalPlane(trueAnomaly);
//
//		// Transform to ECI frame using rotation matrices
//		double x = pos[0] * (cos(raan) * cos(argPerigee) -
//			sin(raan) * sin(argPerigee) * cos(inclination)) -
//			pos[1] * (cos(raan) * sin(argPerigee) +
//				sin(raan) * cos(argPerigee) * cos(inclination));
//
//		double y = pos[0] * (sin(raan) * cos(argPerigee) +
//			cos(raan) * sin(argPerigee) * cos(inclination)) +
//			pos[1] * (cos(raan) * cos(argPerigee) * cos(inclination) -
//				sin(raan) * sin(argPerigee));
//
//		double z = pos[0] * sin(argPerigee) * sin(inclination) +
//			pos[1] * cos(argPerigee) * sin(inclination);
//
//		return { x, y, z };
//	}
//
//	// Get orbital period in seconds
//	double getPeriod() const {
//		return 2 * pi * sqrt(semiMajorAxis * semiMajorAxis * semiMajorAxis / mu);
//	}
//
//	// Get specific orbital energy
//	double getOrbitalEnergy() const {
//		return -mu / (2 * semiMajorAxis);
//	}
//
//	// Getters for orbital elements
//	double getSemiMajorAxis() const { return semiMajorAxis; }
//	double getEccentricity() const { return eccentricity; }
//	double getInclination() const { return inclination; }
//	double getRAAN() const { return raan; }
//	double getArgPerigee() const { return argPerigee; }
//	double getMeanAnomaly() const { return meanAnomaly; }
//};


class OrbitalEllipse {
private:
	// Orbital parameters
	double semiMajorAxis;
	double eccentricity;
	double inclination;
	double raan;         // Right ascension of ascending node
	double argPerigee;   // Argument of periapsis
	int numPoints;       // Number of points to generate for the ellipse

public:
	OrbitalEllipse(double a, double e, double i, double omega, double w, int points = 20)
		: semiMajorAxis(a), eccentricity(e), inclination(i),
		raan(omega), argPerigee(w), numPoints(points) {}

	// Calculate semi-minor axis
	double getSemiMinorAxis() const {
		return semiMajorAxis * sqrt(1.0 - eccentricity * eccentricity);
	}

	// Calculate focus distance from center
	double getFocalDistance() const {
		return semiMajorAxis * eccentricity;
	}

	// Generate points on the ellipse in the orbital plane
	std::vector<std::array<double, 3>> generateEllipsePoints() const {
		std::vector<std::array<double, 3>> points;
		points.reserve(numPoints);

		double semiMinor = getSemiMinorAxis();
		double c = getFocalDistance();

		for (int i = 0; i < numPoints; i++) {
			double theta = 2.0 * pi * i / numPoints;
			double r = semiMajorAxis * (1.0 - eccentricity * eccentricity) /
				(1.0 + eccentricity * cos(theta));

			// Calculate point in orbital plane (perifocal coordinates)
			std::array<double, 3> point = {
				r * cos(theta),
				r * sin(theta),
				0.0
			};

			// Transform to reference coordinates
			points.push_back(transformToReferenceFrame(point));
		}

		return points;
	}

	// Transform a point from the orbital plane to the reference frame
	std::array<double, 3> transformToReferenceFrame(const std::array<double, 3>& point) const {
		// First rotation: argument of periapsis (w)
		double x1 = point[0] * cos(argPerigee) - point[1] * sin(argPerigee);
		double y1 = point[0] * sin(argPerigee) + point[1] * cos(argPerigee);
		double z1 = point[2];

		// Second rotation: inclination (i)
		double x2 = x1;
		double y2 = y1 * cos(inclination) - z1 * sin(inclination);
		double z2 = y1 * sin(inclination) + z1 * cos(inclination);

		// Third rotation: RAAN (Ω)
		double x3 = x2 * cos(raan) - y2 * sin(raan);
		double y3 = x2 * sin(raan) + y2 * cos(raan);
		double z3 = z2;

		return { x3, y3, z3 };
	}

	// Get points at specific true anomalies
	std::array<double, 3> getPointAtTrueAnomaly(double trueAnomaly) const {
		double r = semiMajorAxis * (1.0 - eccentricity * eccentricity) /
			(1.0 + eccentricity * cos(trueAnomaly));

		std::array<double, 3> point = {
			r * cos(trueAnomaly),
			r * sin(trueAnomaly),
			0.0
		};

		return transformToReferenceFrame(point);
	}

	// Get periapsis point
	std::array<double, 3> getPeriapsisPoint() const {
		return getPointAtTrueAnomaly(0.0);
	}

	// Get apoapsis point
	std::array<double, 3> getApoapsisPoint() const {
		return getPointAtTrueAnomaly(pi);
	}

	// Get nodes (points where orbit crosses reference plane)
	std::pair<std::array<double, 3>, std::array<double, 3>> getNodes() const {
		// Ascending node (true anomaly = -w)
		auto ascending = getPointAtTrueAnomaly(-argPerigee);
		// Descending node (true anomaly = π-w)
		auto descending = getPointAtTrueAnomaly(pi - argPerigee);
		return { ascending, descending };
	}

	// Calculate orbital period (if provided with gravitational parameter)
	double getPeriod(double mu) const {
		return 2.0 * pi * sqrt(pow(semiMajorAxis, 3) / mu);
	}

	// Get ellipse parameters for visualization
	struct EllipseParameters {
		double a;        // semi-major axis
		double b;        // semi-minor axis
		double c;        // focal distance
		double theta;    // rotation angle in xy-plane
		double phi;      // rotation angle from z-axis
	};

	EllipseParameters getVisualizationParameters() const {
		return {
			semiMajorAxis,
			getSemiMinorAxis(),
			getFocalDistance(),
			raan,
			inclination
		};
	}
};







struct EllipseParams {
	Eigen::Vector2d center;
	double a;  // semi-major axis
	double b;  // semi-minor axis
	double theta;  // rotation angle
};



EllipseParams findEllipse(const Eigen::Vector2d& p1,
	const Eigen::Vector2d& p2,
	const Eigen::Vector2d& p3,
	const Eigen::Vector2d& focus) {
	// Form the system of equations based on the definition of an ellipse
	// For each point: |r1| + |r2| = 2a where r1, r2 are distances to foci

	// Calculate distances from points to the given focus
	double d1 = (p1).norm();
	double d2 = (p2).norm();
	double d3 = (p3).norm();

	// Form matrix A and vector b for the system Ax = b
	Eigen::Matrix3d A;
	A << 1, p1.x(), p1.y(),
		1, p2.x(), p2.y();
	1, p3.x(), p3.y();

	Eigen::Vector3d b(d1,d2, d3);

	// Solve the system using Jacobi SVD
	Eigen::Vector3d x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

	// Extract the center and semi-axes
	EllipseParams params;
	params.center = focus + Eigen::Vector2d(x[1] / 2, x[2] / 2);
	params.a = x[0] / 2;  // semi-major axis

	// Calculate semi-minor axis using the focus
	double c = (focus - params.center).norm();
	params.b = sqrt(params.a * params.a - c * c);

	// Calculate rotation angle
	params.theta = atan2(focus.y() - params.center.y(),
		focus.x() - params.center.x());

	return params;
}





EllipseParameters global_ep;
vector<array<double, 3Ui64>> double_check_ellipse_points;


double deg_to_rad(double degree)
{
	return degree * (pi/180.0);
}

double hours_to_seconds(double hours)
{
	return hours * 3600.0;
}






struct TrackingData {
	double timestamp;  // Time in seconds
	double azimuth;    // Azimuth in degrees
};

struct OrbitalData {
	double radius;
	double velocity;
};



class cartesian_point {

public:

	double length(void) const
	{
		return sqrt(x * x + y * y);
	}

	double x;
	double y;
};

cartesian_point polarToCartesian(double radius, double azimuth)
{
	cartesian_point result;

	double radians = azimuth;

	// Convert polar to Cartesian
	result.x = radius * std::cos(radians);
	result.y = radius * std::sin(radians);

	return result;
}


struct Vector3D {
	double x, y, z;
};

double distanceBetweenPoints(const Vector3D& p1, const Vector3D& p2) {
	return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2) + std::pow(p2.z - p1.z, 2));
}

void calculateOrbitalParameters(const Vector3D& centralBody, const Vector3D& point1, const Vector3D& point2, double centralMass, double gravitationalConstant)
{
	double r1 = distanceBetweenPoints(centralBody, point1);
	double r2 = distanceBetweenPoints(centralBody, point2);

	// Semi-major axis
	double a = (r1 + r2) / 2.0;

	// Eccentricity
	double e = std::abs(r1 - r2) / (r1 + r2);

	// Focal parameter (distance from center to focus)
	double c = e * a;

	// Orbital period (T) calculation using Kepler's Third Law (T^2 = (4π^2/GM) * a^3)
	double T = std::sqrt((4.0 * pi * pi * std::pow(a, 3)) / (grav_constant * sun_mass));

	std::cout << "Semi-major axis (a): " << a << endl;
	std::cout << "Eccentricity (e): " << e << endl;
	std::cout << "Focal parameter (c): " << c << endl;
	std::cout << "Orbital period (T): " << T << endl;
}



		
void idle_func(void)
{
	static size_t frame_count = 0;

	const double dt = 10000; // 10000 seconds == 2.77777 hours

	proceed_symplectic4(mercury_pos, mercury_vel, grav_constant, dt);

	positions.push_back(mercury_pos);

	static bool calculated_ellipse = false;

	if (calculated_ellipse == false && frame_count % 200 == 0)
		ellipse_positions.push_back(positions[positions.size() - 1]);

	if (false == calculated_ellipse && ellipse_positions.size() == 3)
	{
		calculated_ellipse = true;

		std::vector<TrackingData> measurements = 
		{
			{hours_to_seconds(0),  deg_to_rad(0)},
			{hours_to_seconds(24), deg_to_rad(1)},
			{hours_to_seconds(52), deg_to_rad(8)}
		};

		std::vector<OrbitalData> dataPoints(2);

		for (size_t i = 0; i < measurements.size() - 1; i++)
		{
			//double omega = 4.31e-8; // Ceres average angular velocity
			//double omega = 1.99e-7; // Earth average angular velocity

			// Variable angular velocity
			double omega = (measurements[i + 1].azimuth - measurements[i].azimuth) / (measurements[i + 1].timestamp - measurements[i].timestamp);

			cout << "Angular velocity (omega): " << omega << endl;

			double T = 2 * pi / omega; // Orbital period in seconds
			double r = std::cbrt((sun_mass * grav_constant * T * T) / (4 * pi * pi));
			double v = omega * r;

			// Output results
			std::cout << "Orbital radius (r): " << r << " meters" << std::endl;
			std::cout << "Orbital linear velocity (v): " << v << " meters per second" << std::endl;

			dataPoints[i].radius = r;
			dataPoints[i].velocity = v;
		}


		double angle1 = measurements[1].azimuth;
		double r1 = dataPoints[0].radius;
		double v1 = dataPoints[0].velocity;

		double angle2 = measurements[2].azimuth;
		double r2 = dataPoints[1].radius;
		double v2 = dataPoints[1].velocity;

		cartesian_point cart1 = polarToCartesian(r1, angle1);
		cartesian_point cart2 = polarToCartesian(r2, angle2);

		Vector3D centralBody = { 0, 0, 0 };  // Central body at origin for simplicity
		Vector3D point1 = { cart1.x, cart1.y, 0 };     // First point on the orbit
		Vector3D point2 = { cart2.x, cart2.y, 0 };    // Second point on the orbit
		double centralMass = 1.989e30;     // Mass of sun in kg, for example
		double gravitationalConstant = 6.67430e-11; // m^3 kg^-1 s^-2

		calculateOrbitalParameters(centralBody, point1, point2, centralMass, gravitationalConstant);
	}

	frame_count++;
	glutPostRedisplay();
}

void init_opengl(const int& width, const int& height)
{
	win_x = width;
	win_y = height;

	if (win_x < 1)
		win_x = 1;

	if (win_y < 1)
		win_y = 1;

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("orbit");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	glClearColor((float)background_colour.x, (float)background_colour.y, (float)background_colour.z, 1);
	glClearDepth(1.0f);

	main_camera.Set(0, 0, camera_w, camera_fov, win_x, win_y, camera_near, camera_far);
}

void reshape_func(int width, int height)
{
	win_x = width;
	win_y = height;

	if (win_x < 1)
		win_x = 1;

	if (win_y < 1)
		win_y = 1;

	glutSetWindow(win_id);
	glutReshapeWindow(win_x, win_y);
	glViewport(0, 0, win_x, win_y);

	main_camera.Set(main_camera.u, main_camera.v, main_camera.w, main_camera.fov, win_x, win_y, camera_near, camera_far);
}

// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
void render_string(int x, const int y, void* font, const string& text)
{
	for (size_t i = 0; i < text.length(); i++)
	{
		glRasterPos2i(x, y);
		glutBitmapCharacter(font, text[i]);
		x += glutBitmapWidth(font, text[i]) + 1;
	}
}
// End text drawing code.

void draw_objects(void)
{
	glDisable(GL_LIGHTING);

	glPushMatrix();


	glPointSize(6.0);
	glLineWidth(1.0f);


	glBegin(GL_POINTS);
	glVertex3d(sun_pos.x, sun_pos.y, sun_pos.z);

	glColor3f(1.0, 1.0, 1.0);

	for (size_t i = 0; i < ellipse_positions.size(); i++)
		glVertex3d(ellipse_positions[i].x, ellipse_positions[i].y, ellipse_positions[i].z);

	glEnd();





	glColor3f(1.0, 0.5, 0.0);

	glPushMatrix();

	glTranslated(global_ep.centerX, global_ep.centerY, 0);
	glRotated(global_ep.angle / (2 * pi) * 360.0, 0, 0, 1);
	glTranslated(-global_ep.centerX, -global_ep.centerY, 0);

	DrawEllipse(global_ep.centerX, global_ep.centerY, global_ep.semiMinor, global_ep.semiMajor, 100);

	glPopMatrix();



	glBegin(GL_POINTS);

	glColor3f(1.0, 0.0, 0.0);

	for (size_t i = 0; i < double_check_ellipse_points.size(); i++)
		glVertex3d(double_check_ellipse_points[i][0], double_check_ellipse_points[i][1], double_check_ellipse_points[i][2]);

	glEnd();



	//glBegin(GL_TRIANGLES);

	//glColor3f(1.0, 1.0, 1.0);

	//for (size_t i = 0; i < triangles.size(); i++)
	//{
	//	glVertex3d(triangles[i].vertex[0].x, triangles[i].vertex[0].y, 0);
	//	glVertex3d(triangles[i].vertex[1].x, triangles[i].vertex[1].y, 0);
	//	glVertex3d(triangles[i].vertex[2].x, triangles[i].vertex[2].y, 0);

	//}

	//glEnd();


	//glLineWidth(1.0f);

	//glBegin(GL_LINES);
	//glColor3f(1.0, 0.5, 0.0);

	//for (size_t i = 0; i < line_segments.size(); i++)
	//{
	//	glVertex3d(line_segments[i].vertex[0].x, line_segments[i].vertex[0].y, 0);
	//	glVertex3d(line_segments[i].vertex[1].x, line_segments[i].vertex[1].y, 0);
	//}

	//glEnd();


 //   
 //   
	//// If we do draw the axis at all, make sure not to draw its outline.
	//if(true == draw_axis)
	//{
	//	glBegin(GL_LINES);

	//	glColor3f(1, 0, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(1, 0, 0);
	//	glColor3f(0, 1, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 1, 0);
	//	glColor3f(0, 0, 1);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 0, 1);

	//	glColor3f(0.5, 0.5, 0.5);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(-1, 0, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, -1, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 0, -1);

	//	glEnd();
	//}

	glPopMatrix();
}




void display_func(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw the model's components using OpenGL/GLUT primitives.
	draw_objects();

	if (true == draw_control_list)
	{
		// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
		// http://www.lighthouse3d.com/opengl/glut/index.php?bmpfontortho
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, win_x, 0, win_y);
		glScaled(1, -1, 1); // Neat. :)
		glTranslated(0, -win_y, 0); // Neat. :)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor3d(control_list_colour.x, control_list_colour.y, control_list_colour.z);

		size_t break_size = 22;
		size_t start = 20;
		ostringstream oss;

		render_string(10, static_cast<int>(start), GLUT_BITMAP_HELVETICA_18, string("Mouse controls:"));
		render_string(10, static_cast<int>(start + 1 * break_size), GLUT_BITMAP_HELVETICA_18, string("  LMB + drag: Rotate camera"));
		render_string(10, static_cast<int>(start + 2 * break_size), GLUT_BITMAP_HELVETICA_18, string("  RMB + drag: Zoom camera"));

		render_string(10, static_cast<int>(start + 4 * break_size), GLUT_BITMAP_HELVETICA_18, string("Keyboard controls:"));
		render_string(10, static_cast<int>(start + 5 * break_size), GLUT_BITMAP_HELVETICA_18, string("  w: Draw axis"));
		render_string(10, static_cast<int>(start + 6 * break_size), GLUT_BITMAP_HELVETICA_18, string("  e: Draw text"));
		render_string(10, static_cast<int>(start + 7 * break_size), GLUT_BITMAP_HELVETICA_18, string("  u: Rotate camera +u"));
		render_string(10, static_cast<int>(start + 8 * break_size), GLUT_BITMAP_HELVETICA_18, string("  i: Rotate camera -u"));
		render_string(10, static_cast<int>(start + 9 * break_size), GLUT_BITMAP_HELVETICA_18, string("  o: Rotate camera +v"));
		render_string(10, static_cast<int>(start + 10 * break_size), GLUT_BITMAP_HELVETICA_18, string("  p: Rotate camera -v"));



		custom_math::vector_3 eye = main_camera.eye;
		custom_math::vector_3 eye_norm = eye;
		eye_norm.normalize();

		oss.clear();
		oss.str("");
		oss << "Camera position: " << eye.x << ' ' << eye.y << ' ' << eye.z;
		render_string(10, static_cast<int>(win_y - 2 * break_size), GLUT_BITMAP_HELVETICA_18, oss.str());

		oss.clear();
		oss.str("");
		oss << "Camera position (normalized): " << eye_norm.x << ' ' << eye_norm.y << ' ' << eye_norm.z;
		render_string(10, static_cast<int>(win_y - break_size), GLUT_BITMAP_HELVETICA_18, oss.str());

		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		// End text drawing code.
	}

	glFlush();
	glutSwapBuffers();
}

void keyboard_func(unsigned char key, int x, int y)
{
	switch (tolower(key))
	{
	case 'w':
	{
		draw_axis = !draw_axis;
		break;
	}
	case 'e':
	{
		draw_control_list = !draw_control_list;
		break;
	}
	case 'u':
	{
		main_camera.u -= u_spacer;
		main_camera.Set();
		break;
	}
	case 'i':
	{
		main_camera.u += u_spacer;
		main_camera.Set();
		break;
	}
	case 'o':
	{
		main_camera.v -= v_spacer;
		main_camera.Set();
		break;
	}
	case 'p':
	{
		main_camera.v += v_spacer;
		main_camera.Set();
		break;
	}

	default:
		break;
	}
}

void mouse_func(int button, int state, int x, int y)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		if (GLUT_DOWN == state)
			lmb_down = true;
		else
			lmb_down = false;
	}
	else if (GLUT_MIDDLE_BUTTON == button)
	{
		if (GLUT_DOWN == state)
			mmb_down = true;
		else
			mmb_down = false;
	}
	else if (GLUT_RIGHT_BUTTON == button)
	{
		if (GLUT_DOWN == state)
			rmb_down = true;
		else
			rmb_down = false;
	}
}

void motion_func(int x, int y)
{
	int prev_mouse_x = mouse_x;
	int prev_mouse_y = mouse_y;

	mouse_x = x;
	mouse_y = y;

	int mouse_delta_x = mouse_x - prev_mouse_x;
	int mouse_delta_y = prev_mouse_y - mouse_y;

	if (true == lmb_down && (0 != mouse_delta_x || 0 != mouse_delta_y))
	{
		main_camera.u -= static_cast<float>(mouse_delta_y) * u_spacer;
		main_camera.v += static_cast<float>(mouse_delta_x) * v_spacer;
	}
	else if (true == rmb_down && (0 != mouse_delta_y))
	{
		main_camera.w -= static_cast<float>(mouse_delta_y) * w_spacer;

		if (main_camera.w < 1.1f)
			main_camera.w = 1.1f;

	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}



