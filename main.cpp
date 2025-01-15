#include "main.h"


#include <Eigen/Dense>
using namespace Eigen;



class cartesian_point
{
public:

	double x;
	double y;

	double length(void) const
	{
		return sqrt(x * x + y * y);
	}

	cartesian_point operator-(const cartesian_point& rhs) const
	{
		cartesian_point ret;
		ret.x = x - rhs.x;
		ret.y = y - rhs.y;

		return ret;
	}
};

//struct Point {
//	double x, y;
//};

class EllipseParameters
{
public:
	double centerX = 0;
	double centerY = 0;
	double semiMajor = 0;
	double semiMinor = 0;
	double angle = 0;
};


EllipseParameters global_ep;
vector<array<double, 3Ui64>> double_check_ellipse_points;







EllipseParameters extractEllipseParameters(const Eigen::VectorXd& coefficients)
{
	double a = coefficients(0);
	double b = coefficients(1);
	double c = coefficients(2);
	double d = coefficients(3);
	double e = coefficients(4);
	double f = 1;// coefficients(5);

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






EllipseParameters fitEllipse(const std::vector<cartesian_point>& points, const cartesian_point& focus)
{
	if (points.size() != 5) {
		std::cerr << "Error: Exactly 5 points are required.\n";
		return EllipseParameters();
	}

	Eigen::MatrixXd A(5, 6);
	Eigen::VectorXd b(5);

	// Fill the matrix A and vector b with the equations from the points
	for (size_t i = 0; i < 5; ++i)
	{
		double x = points[i].x;
		double y = points[i].y;
		A(i, 0) = x * x;       // Coefficient for x^2
		A(i, 1) = x * y;       // Coefficient for xy
		A(i, 2) = y * y;       // Coefficient for y^2
		A(i, 3) = x;           // Coefficient for x
		A(i, 4) = y;           //  Coefficient for y
		A(i, 5) = 1;           // Constant term
		b(i) = -1;             // Right-hand side is -1. This is important!
	}

	// Solve for the ellipse parameters
	Eigen::VectorXd ellipseParams = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

	// Extract parameters
	double A_ = ellipseParams(0);
	double B_ = ellipseParams(1);
	double C_ = ellipseParams(2);
	double D_ = ellipseParams(3);
	double E_ = ellipseParams(4);
	double F_ = ellipseParams(5);




	// Compute center of ellipse
	EllipseParameters ep = extractEllipseParameters(ellipseParams);


	global_ep.angle = ep.angle;
	global_ep.centerX = ep.centerX;
	global_ep.centerY = ep.centerY;
	global_ep.semiMajor = ep.semiMajor;
	global_ep.semiMinor = ep.semiMinor;

	//cout << global_ep.angle << endl;
	//cout << global_ep.centerX << endl;
	//cout << global_ep.centerY << endl;
	//cout << global_ep.semiMajor << endl;
	//cout << global_ep.semiMinor << endl;


	return ep;
}
// https://chatgpt.com/c/67809363-9d0c-8003-b1cc-ad9a37e10e54


struct EllipseParams_min {
	double semiMajorAxis;
	double eccentricity;
	double angle; // in radians
};

void gatherConstraints(
	const std::vector<cartesian_point>& points,
	const std::vector<cartesian_point>& velocities,
	Eigen::MatrixXd& A,
	Eigen::VectorXd& b)
{
	size_t n = points.size();
	A = Eigen::MatrixXd(2 * n, 6);
	b = Eigen::VectorXd::Ones(2 * n);
	b = -b; // Make them all -1s

	for (int i = 0; i < n; ++i) {
		double x = points[i].x, y = points[i].y;
		double vx = velocities[i].x, vy = velocities[i].y;

		// Ellipse equation constraint
		A(i, 0) = x * x;      // A
		A(i, 1) = x * y;      // B
		A(i, 2) = y * y;      // C
		A(i, 3) = x;          // D
		A(i, 4) = y;          // E
		A(i, 5) = 1;          // F

		// Tangent velocity constraint
		A(n + i, 0) = 2 * x * vx; // 2A x vx
		A(n + i, 1) = x * vy + y * vx; // B (x vy + y vx)
		A(n + i, 2) = 2 * y * vy; // 2C y vy
		A(n + i, 3) = vx;         // D vx
		A(n + i, 4) = vy;         // E vy
		A(n + i, 5) = 0;          // F
	}
}


std::vector<double> solveEllipseCoefficients(
	const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b)
{
	Eigen::VectorXd coefficients = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	return std::vector<double>(coefficients.data(), coefficients.data() + coefficients.size());
}

EllipseParams_min extractEllipseParams(
	const std::vector<double>& coefficients,
	const cartesian_point& focus)
{
	double A = coefficients[0], B = coefficients[1], C = coefficients[2];
	double D = coefficients[3], E = coefficients[4], F = coefficients[5];

	// Calculate orientation angle (theta)
	double theta = 0.5 * atan2(B, A - C);

	// Transform coefficients to canonical form
	double cosTheta = cos(theta), sinTheta = sin(theta);
	double Ap = A * cosTheta * cosTheta + B * cosTheta * sinTheta + C * sinTheta * sinTheta;
	double Cp = A * sinTheta * sinTheta - B * cosTheta * sinTheta + C * cosTheta * cosTheta;

	// Calculate semi-major and semi-minor axes
	double a = sqrt(1 / fabs(Ap));
	double b = sqrt(1 / fabs(Cp));

	// Calculate eccentricity
	double e = sqrt(1 - (b * b) / (a * a));

	double centerX = (2 * C * D - B * E) / (B * B - 4 * A * C);
	double centerY = (2 * A * E - B * D) / (B * B - 4 * A * C);

	global_ep.angle = theta + pi / 2;
	global_ep.centerX = centerX;
	global_ep.centerY = centerY;
	global_ep.semiMajor = a;
	global_ep.semiMinor = b;// ep.semiMinor;

	//cout << global_ep.angle << endl;
	//cout << global_ep.centerX << endl;
	//cout << global_ep.centerY << endl;
	//cout << global_ep.semiMajor << endl;
	//cout << global_ep.semiMinor << endl;

	return { a, e, theta };
}


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








struct EllipseParameters2 {
	double semi_major_axis;
	double eccentricity;
	double angle;
	Eigen::Vector2d center;
};


Eigen::Vector2d calculateSecondFocus(
	const EllipseParameters2& params,
	const Eigen::Vector2d& focus1)
{

	// Validate parameters
	if (params.eccentricity < 0.0 || params.eccentricity >= 1.0) {
		return Vector2d();// throw std::invalid_argument("Eccentricity must be in [0,1)");
	}
	if (params.semi_major_axis <= 0.0) {
		return Vector2d();// throw std::invalid_argument("Semi-major axis must be positive");
	}

	// Calculate focal distance (distance from center to focus)
	double focal_distance = params.semi_major_axis * params.eccentricity;

	// Vector from center to first focus
	Eigen::Vector2d center_to_focus1 = focus1 - params.center;

	// Verify first focus is at correct distance
	//double actual_distance = center_to_focus1.norm();
	//if (std::abs(actual_distance - focal_distance) > 1e-10) {
	//	throw std::invalid_argument("First focus is not at correct distance from center");
	//}

	// Second focus is opposite to first focus through center
	Eigen::Vector2d focus2 = params.center - center_to_focus1;

	return focus2;
}



Eigen::VectorXd ellipseParamsToCoefficients(double h, double k, double a, double b, double theta) {
	Eigen::VectorXd coeffs(6);  // a, b, c, d, e, f

	double cosTheta = cos(theta);
	double sinTheta = sin(theta);
	double a2 = a * a;
	double b2 = b * b;

	// Coefficients from the implicit equation of an ellipse
	coeffs[0] = b2 * cosTheta * cosTheta + a2 * sinTheta * sinTheta;  // a
	coeffs[1] = 2 * (b2 - a2) * sinTheta * cosTheta;              // b
	coeffs[2] = b2 * sinTheta * sinTheta + a2 * cosTheta * cosTheta;  // c
	coeffs[3] = -2 * coeffs[0] * h - coeffs[1] * k;               // d
	coeffs[4] = -coeffs[1] * h - 2 * coeffs[2] * k;               // e
	coeffs[5] = coeffs[0] * h * h + coeffs[2] * k * k + coeffs[1] * h * k - a2 * b2;  // f

	return coeffs;
}












// Helper function for distance calculation
double distance(double x1, double y1, double x2, double y2)
{
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

// Objective function for optimization (least squares)
double objectiveFunction(
	const VectorXd& params,
	const vector<cartesian_point>& points,
	const vector<cartesian_point>& velocities,
	const cartesian_point& focus)
{
	double h = params[0], k = params[1], a = params[2], b = params[3], angle = params[4];
	double error = 0;

	EllipseParameters2 ep;
	ep.angle = angle;
	ep.center(0) = 0;
	ep.center(1) = k;
	ep.semi_major_axis = a;
	ep.eccentricity = sqrt(1 - (b * b) / (a * a));

	Eigen::Vector2d focus1;
	focus1(0) = focus.x;
	focus1(1) = focus.y;

	Eigen::Vector2d focus2 = calculateSecondFocus(ep, focus1);

	for (size_t i = 0; i < points.size(); i++)
	{
		cartesian_point p = points[i];
		cartesian_point v = velocities[i];

		double dist1 = distance(p.x, p.y, focus1(0), focus1(1));
		double dist2 = distance(p.x, p.y, focus2(0), focus2(1));

		error += pow(dist1 + dist2 - 2 * a, 2);

		cartesian_point f2(focus2(0), focus2(1));
		cartesian_point centre(h, k);
		double focal_dist = (f2 - centre).length();
		double c = std::sqrt(std::abs(a * a - b * b));  // focal distance from center

		error += std::pow(focal_dist - c, 2.0);

		// Since we're axis-aligned, we simplify velocity condition:
		// Velocity should be more in line with the axis of the ellipse
		//double velError = 0;

		//if (abs(v.x) > abs(v.y)) // Suggesting a is along x
		//	velError = pow(v.y / v.x - (k - p.y) / (h - p.x), 2); // Check alignment with y
		//else
		//	velError = pow(v.x / v.y - (h - p.x) / (k - p.y), 2); // Check alignment with x

		//error += velError;
	}

	return error;
}

// Simple solver function using gradient descent (for demonstration)
VectorXd solveEllipseParameters(const vector<cartesian_point>& points, const vector<cartesian_point>& velocities, const cartesian_point& focus)
{
	// Get max distance data
	vector<double> mvec;

	mvec.push_back(points[0].length());
	mvec.push_back(points[1].length());
	mvec.push_back(points[2].length());
	mvec.push_back(points[3].length());
	mvec.push_back(points[4].length());

	sort(mvec.begin(), mvec.end());

	// Use the maximum distance data	
	const double m = mvec[4];

	double d = (mvec[4] - mvec[0]) / mvec[4];
	d = pow(1 - d, 20.0);

	VectorXd params(5); // h, k, a, b, angle
	params << 1, 1, m* d, m* d, 0; // Initial guess

	int iterations = 100000;
	double stepSize = 0.0001;

	for (int i = 0; i < iterations; i++)
	{
		VectorXd gradient = VectorXd::Zero(5);

		for (int j = 0; j < 5; j++)
		{
			VectorXd paramsPlus = params;
			paramsPlus[j] += stepSize;
			VectorXd paramsMinus = params;
			paramsMinus[j] -= stepSize;

			gradient[j] = (objectiveFunction(paramsPlus, points, velocities, focus) - objectiveFunction(paramsMinus, points, velocities, focus)) / (2 * stepSize);
		}

		params -= stepSize * gradient;
	}

	return params;
}


//
//
//
//
//const double GM = sun_mass * grav_constant; // km^3/s^2, gravitational constant * mass of Earth
//
//namespace ns
//{
//
//
//
//	struct Vector3 {
//		double x, y, z;
//
//		Vector3 operator+(const Vector3& other) const {
//			return { x + other.x, y + other.y, z + other.z };
//		}
//
//		Vector3 operator-(const Vector3& other) const {
//			return { x - other.x, y - other.y, z - other.z };
//		}
//
//		Vector3 operator*(double scalar) const {
//			return { x * scalar, y * scalar, z * scalar };
//		}
//
//		Vector3 operator/(double scalar) const {
//			return { x / scalar, y / scalar, z / scalar };
//		}
//
//		double dot(const Vector3& other) const {
//			return x * other.x + y * other.y + z * other.z;
//		}
//
//		Vector3 cross(const Vector3& other) const {
//			return { y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x };
//		}
//
//		double magnitude() const {
//			return std::sqrt(x * x + y * y + z * z);
//		}
//	};
//}
//
//struct OrbitElements {
//	double semi_major_axis;
//	double eccentricity;
//	double inclination;
//	double raan; // Right ascension of the ascending node
//	double arg_periapsis;
//	double true_anomaly;
//};
//
//ns::Vector3 solve_gauss(const ns::Vector3& r1, const ns::Vector3& r2, const ns::Vector3& r3, double dt21, double dt31) {
//	double r1_mag = r1.magnitude();
//	double r2_mag = r2.magnitude();
//	double r3_mag = r3.magnitude();
//
//	double cos_theta21 = r1.dot(r2) / (r1_mag * r2_mag);
//	double cos_theta31 = r1.dot(r3) / (r1_mag * r3_mag);
//
//	double sin_theta21 = std::sqrt(1 - cos_theta21 * cos_theta21);
//	double sin_theta31 = std::sqrt(1 - cos_theta31 * cos_theta31);
//
//	double a1 = sin_theta31 / sin_theta31 - cos_theta31;
//	double a3 = sin_theta21 / sin_theta31 - cos_theta21;
//
//	ns::Vector3 d1 = r2 * a1;
//	ns::Vector3 d3 = r2 * a3;
//
//	ns::Vector3 d = d1 + d3;
//
//	double d_mag = d.magnitude();
//
//	ns::Vector3 v2 = (d / dt21 - r1 / dt21) * GM / (r2_mag * d_mag);
//
//	return v2;
//}
//
//OrbitElements calc_orbit_elements(const ns::Vector3& r, const ns::Vector3& v) {
//	ns::Vector3 h = r.cross(v);
//	double h_mag = h.magnitude();
//
//	ns::Vector3 e = ((v.cross(h)) / GM) - (r / r.magnitude());
//	double e_mag = e.magnitude();
//
//	double semi_major_axis = 1 / ((2 / r.magnitude()) - (v.magnitude() * v.magnitude() / GM));
//
//	double inclination = std::acos(h.z / h_mag);
//
//	double raan = std::atan2(h.x, -h.y);
//
//	double arg_periapsis = std::atan2(e.z / std::sin(inclination), e.x * std::cos(raan) + e.y * std::sin(raan));
//
//	double true_anomaly = std::atan2(r.dot(v) / std::sqrt(GM * semi_major_axis), 1 - r.magnitude() / semi_major_axis);
//
//	return { semi_major_axis, e_mag, inclination, raan, arg_periapsis, true_anomaly };
//}
//
//




//
//
//// Define constants (you would replace these with actual data)
//const double mu = sun_mass*grav_constant; // Earth's gravitational parameter in km^3/s^2
//
//// Struct to represent a vector in 2D space (since z is always zero)
//struct vector_2d {
//	double x;
//	double y;
//};
//
//// Struct to represent an observation
//struct observation {
//	double time;
//	vector_2d position;
//};
//
//// Function to compute the magnitude of a 2D vector
//double calculate_magnitude(const vector_2d& v) {
//	return std::sqrt(v.x * v.x + v.y * v.y);
//}
//
//// Function to compute the dot product of two 2D vectors
//double dot_product(const vector_2d& a, const vector_2d& b) {
//	return a.x * b.x + a.y * b.y;
//}
//
//// Function to calculate cross product in 2D (which gives a scalar in z direction)
//double cross_product(const vector_2d& a, const vector_2d& b) {
//	return a.x * b.y - a.y * b.x;
//}
//
//// Function to perform Gauss's method for orbit determination (velocity at t2)
//vector_2d gauss_method_orbit_determination(const std::vector<observation>& observations) {
//	if (observations.size() != 3) {
//		std::cerr << "Error: Need exactly three observations for Gauss's method." << std::endl;
//		return { 0, 0 }; // Return zero vector on error
//	}
//
//	// Time intervals
//	double tau1 = observations[1].time - observations[0].time;
//	double tau3 = observations[2].time - observations[1].time;
//	double tau = tau1 + tau3;
//
//	// Position vectors
//	vector_2d r1 = observations[0].position;
//	vector_2d r2 = observations[1].position;
//	vector_2d r3 = observations[2].position;
//
//	// Compute necessary constants
//	double D = cross_product(r1, r3);
//	double E = cross_product(r2, r3) + cross_product(r1, r2);
//	double F = 3 * dot_product(r2, r2) + mu * tau * tau;
//	double G = dot_product(r1, r3);
//
//	// Denominator for solving lambda
//	double d = 2 * F * E - D * (F + 2 * G);
//
//	// Solve for lambda (assuming positive solution for simplicity)
//	double lambda = (D * F * tau1 - 2 * mu * tau * tau * tau3 * E) / d;
//
//	// Solve for velocity at time t2 (middle observation)
//	vector_2d velocity = {
//		(2 * (cross_product(r3, r2) + lambda * r3.y - lambda * r2.y)) / (tau3 * (tau1 + tau3)),
//		(-2 * (cross_product(r3, r2) + lambda * r3.x - lambda * r2.x)) / (tau3 * (tau1 + tau3))
//	};
//
//	return velocity;
//}
//
//// Compute orbital parameters
//void compute_orbit_parameters(const vector_2d& position, const vector_2d& velocity, double& semi_major_axis, double& eccentricity) {
//	double r = calculate_magnitude(position);
//	double v = calculate_magnitude(velocity);
//	double vr = dot_product(position, velocity) / r; // Radial component of velocity
//
//	// Specific angular momentum
//	double h = r * std::sqrt(v * v - vr * vr);
//
//	// Semi-major axis calculation
//	semi_major_axis = 1.0 / (2.0 / r - v * v / mu);
//
//	// Eccentricity vector
//	vector_2d e_vector = {
//		(v * v / mu - 1.0 / r) * position.x - (r * vr / mu) * velocity.x,
//		(v * v / mu - 1.0 / r) * position.y - (r * vr / mu) * velocity.y
//	};
//
//	eccentricity = calculate_magnitude(e_vector);
//
//	// Other parameters could be calculated here like inclination, argument of perigee etc., but not needed for equatorial plane
//}
//
//const double gravitational_constant = grav_constant*sun_mass; // Earth's gravitational constant in km^3/s^2
//
//struct vector3d {
//	double x;
//	double y;
//	double z;
//};
//
//struct orbital_elements {
//	double semi_major_axis;
//	double eccentricity;
//	double inclination;
//	double longitude_of_ascending_node;
//	double argument_of_periapsis;
//	double true_anomaly;
//};
//
//vector3d vector_subtract(const vector3d& a, const vector3d& b) {
//	return { a.x - b.x, a.y - b.y, a.z - b.z };
//}
//
//vector3d vector_add(const vector3d& a, const vector3d& b) {
//	return { a.x + b.x, a.y + b.y, a.z + b.z };
//}
//
//double vector_dot(const vector3d& a, const vector3d& b) {
//	return a.x * b.x + a.y * b.y + a.z * b.z;
//}
//
//vector3d vector_cross(const vector3d& a, const vector3d& b) {
//	return {
//		a.y * b.z - a.z * b.y,
//		a.z * b.x - a.x * b.z,
//		a.x * b.y - a.y * b.x
//	};
//}
//
//double vector_magnitude(const vector3d& v) {
//	return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
//}
//
//vector3d vector_scale(const vector3d& v, double scalar) {
//	return { v.x * scalar, v.y * scalar, v.z * scalar };
//}
//
//vector3d vector_divide(const vector3d& v, double scalar) {
//	return { v.x / scalar, v.y / scalar, v.z / scalar };
//}
//
//orbital_elements gauss_method_orbital_elements(const vector3d& position, const vector3d& velocity) {
//	orbital_elements elements;
//
//	// Calculate specific angular momentum
//	vector3d h = vector_cross(position, velocity);
//	double h_magnitude = vector_magnitude(h);
//
//	// Node vector
//	vector3d n = { -h.y, h.x, 0.0 };
//	double n_magnitude = vector_magnitude(n);
//
//	// Eccentricity vector
//	vector3d r_unit = vector_scale(position, 1.0 / vector_magnitude(position));
//	vector3d v_cross_h = vector_cross(velocity, h);
//	vector3d e_vector = vector_subtract(vector_scale(v_cross_h, 1.0 / gravitational_constant), r_unit);
//	elements.eccentricity = vector_magnitude(e_vector);
//
//	// Semi-major axis
//	double r = vector_magnitude(position);
//	double v = vector_magnitude(velocity);
//	double energy = (v * v) / 2.0 - gravitational_constant / r;
//	elements.semi_major_axis = -gravitational_constant / (2.0 * energy);
//
//	// Inclination
//	elements.inclination = std::acos(h.z / h_magnitude);
//
//	// Longitude of ascending node
//	if (n_magnitude != 0) {
//		elements.longitude_of_ascending_node = std::acos(n.x / n_magnitude);
//		if (n.y < 0) {
//			elements.longitude_of_ascending_node = 2 * pi - elements.longitude_of_ascending_node;
//		}
//	}
//	else {
//		elements.longitude_of_ascending_node = 0.0;
//	}
//
//	// Argument of periapsis
//	if (n_magnitude != 0 && elements.eccentricity != 0) {
//		double arg_peri = vector_dot(n, e_vector) / (n_magnitude * elements.eccentricity);
//		arg_peri = std::acos(arg_peri);
//		if (e_vector.z < 0) {
//			arg_peri = 2 * pi - arg_peri;
//		}
//		elements.argument_of_periapsis = arg_peri;
//	}
//	else {
//		elements.argument_of_periapsis = 0.0;
//	}
//
//	// True anomaly
//	double true_anom = vector_dot(e_vector, position) / (elements.eccentricity * r);
//	true_anom = std::acos(true_anom);
//	if (vector_dot(position, velocity) < 0) {
//		true_anom = 2 * pi - true_anom;
//	}
//	elements.true_anomaly = true_anom;
//
//	return elements;
//}
//
//
//orbital_elements gauss_method_with_three_points(const vector3d& r1, const vector3d& r2, const vector3d& r3, double t1, double t2, double t3) {
//	orbital_elements elements;
//
//	// Time intervals
//	double tau1 = t2 - t1;
//	double tau3 = t3 - t2;
//	double tau = tau1 + tau3;
//
//	// Calculate the coefficients for Gauss's method
//	double f1 = 1 - (gravitational_constant * tau1 * tau1) / (2 * std::pow(vector_magnitude(r2), 3));
//	double f3 = 1 - (gravitational_constant * tau3 * tau3) / (2 * std::pow(vector_magnitude(r2), 3));
//	double g1 = tau1 - (gravitational_constant * std::pow(tau1, 3)) / (6 * std::pow(vector_magnitude(r2), 3));
//	double g3 = tau3 - (gravitational_constant * std::pow(tau3, 3)) / (6 * std::pow(vector_magnitude(r2), 3));
//
//	// Velocity at the second point
//	vector3d v2 = vector_divide(vector_subtract(vector_scale(r1, -f3), vector_scale(r3, f1)), g1 * f3 - g3 * f1);
//
//	// Use r2 and v2 to calculate the orbital elements
//	elements = gauss_method_orbital_elements(r2, v2);
//
//	return elements;
//}
//
//
//void print_orbital_elements(const orbital_elements& elements) {
//	std::cout << std::fixed << std::setprecision(6);
//	std::cout << "Semi-Major Axis (km): " << elements.semi_major_axis << std::endl;
//	std::cout << "Eccentricity: " << elements.eccentricity << std::endl;
//	std::cout << "Inclination (rad): " << elements.inclination << std::endl;
//	std::cout << "Longitude of Ascending Node (rad): " << elements.longitude_of_ascending_node << std::endl;
//	std::cout << "Argument of Periapsis (rad): " << elements.argument_of_periapsis << std::endl;
//	std::cout << "True Anomaly (rad): " << elements.true_anomaly << std::endl;
//}









struct vector3d {
	double x;
	double y;
	double z;

	vector3d(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}

	vector3d operator+(const vector3d& v) const {
		return vector3d(x + v.x, y + v.y, z + v.z);
	}

	vector3d operator-(const vector3d& v) const {
		return vector3d(x - v.x, y - v.y, z - v.z);
	}

	vector3d operator*(double scalar) const {
		return vector3d(x * scalar, y * scalar, z * scalar);
	}

	vector3d operator/(double scalar) const {
		return vector3d(x / scalar, y / scalar, z / scalar);
	}

	vector3d cross(const vector3d& v) const {
		return vector3d(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}
};

// Function to compute Lagrange coefficients
void compute_lagrange_coefficients(double t0, double t1, double t2, double mu,
	double& f0, double& f1, double& f2,
	double& g0, double& g1, double& g2) {
	double tau0 = t1 - t0;
	double tau1 = t2 - t1;
	double tau2 = t2 - t0;

	double r0 = std::sqrt(tau1 / tau2);
	double r1 = std::sqrt(tau0 / tau2);
	double r2 = std::sqrt(tau0 / tau1);

	f0 = 1.0 - (tau1 / tau2);
	f1 = 1.0 - (tau0 * tau2) / (tau1 * tau2);
	f2 = 1.0 - (tau0 / tau1);

	g0 = (tau1 * r0) / mu;
	g1 = (tau0 * tau2) / (mu * (r1 * r2));
	g2 = (tau0 * r2) / mu;
}

// Function to solve for position and velocity at intermediate time
vector3d solve_position(vector3d r0, vector3d r1, vector3d r2,
	double t0, double t1, double t2, double mu) {
	double f0, f1, f2, g0, g1, g2;
	compute_lagrange_coefficients(t0, t1, t2, mu, f0, f1, f2, g0, g1, g2);

	return (r0 * f0 + r1 * f1 + r2 * f2) / (f0 + f1 + f2);
}

vector3d solve_velocity(vector3d r0, vector3d r1, vector3d r2,
	double t0, double t1, double t2, double mu) {
	double f0, f1, f2, g0, g1, g2;
	compute_lagrange_coefficients(t0, t1, t2, mu, f0, f1, f2, g0, g1, g2);

	return (r0 * g0 + r1 * g1 + r2 * g2) / (g0 + g1 + g2);
}

// Main function to determine orbit parameters
void determine_orbit_parameters(vector3d r1, vector3d r2, vector3d r3,
	double t1, double t2, double t3, double mu) {
	vector3d r_mid = solve_position(r1, r2, r3, t1, t2, t3, mu);
	vector3d v_mid = solve_velocity(r1, r2, r3, t1, t2, t3, mu);

	// Compute angular momentum
	vector3d h = r_mid.cross(v_mid);

	// Compute eccentricity vector (assuming equatorial plane, z component of h should be zero)
	vector3d e = (v_mid.cross(h) / mu) - (r_mid / std::sqrt(r_mid.x * r_mid.x + r_mid.y * r_mid.y));

	// Semi-major axis
	double a = 1.0 / (2.0 / std::sqrt(r_mid.x * r_mid.x + r_mid.y * r_mid.y) - (v_mid.x * v_mid.x + v_mid.y * v_mid.y) / mu);

	// Eccentricity magnitude
	double e_mag = std::sqrt(e.x * e.x + e.y * e.y);

	// Output parameters
	std::cout << "Semi-major axis (a): " << a << std::endl;
	std::cout << "Eccentricity (e): " << e_mag << std::endl;
	std::cout << "Inclination (i): 0.0 (equatorial plane)" << std::endl;
}




//
//#include <vector>
//#include <cmath>
//#include <array>
//#include <iostream>
//#include <stdexcept>
//
//struct vector_3d {
//	double x, y, z;
//
//	vector_3d(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
//
//	vector_3d operator+(const vector_3d& v) const {
//		return vector_3d(x + v.x, y + v.y, z + v.z);
//	}
//
//	vector_3d operator-(const vector_3d& v) const {
//		return vector_3d(x - v.x, y - v.y, z - v.z);
//	}
//
//	vector_3d operator*(double scalar) const {
//		return vector_3d(x * scalar, y * scalar, z * scalar);
//	}
//
//	vector_3d operator/(double scalar) const {
//		return vector_3d(x / scalar, y / scalar, z / scalar);
//	}
//
//	double get_magnitude() const {
//		return std::sqrt(x * x + y * y + z * z);
//	}
//
//	vector_3d normalize() const {
//		double mag = get_magnitude();
//		return vector_3d(x / mag, y / mag, z / mag);
//	}
//};
//
//struct orbit_parameters {
//	double semi_major_axis;      // a
//	double eccentricity;         // e
//	double true_anomaly;         // ν
//	double mean_motion;          // n
//	double period;               // T
//	double mean_anomaly;         // M
//	double eccentric_anomaly;    // E
//	vector_3d angular_momentum;  // h
//	vector_3d eccentricity_vector; // e vector
//};
//
//class gauss_orbit_calculator {
//private:
//	const double gravitational_parameter = sun_mass*grav_constant; // Earth's GM in km³/s²
//	const double tolerance = 1e-12;
//	const int max_iterations = 1000;
//
//	struct internal_orbit_data {
//		double tau;          // time difference τ
//		double tau1;         // τ₁
//		double tau3;         // τ₃
//		double p;           // semi-parameter
//		std::array<double, 3> rho;  // range values ρ
//		std::array<vector_3d, 3> r;  // position vectors
//	};
//
//	double calculate_stumpff_c(double z) {
//		if (std::abs(z) < tolerance) {
//			return 1.0 / 2.0;
//		}
//		else if (z > 0) {
//			return (1.0 - std::cos(std::sqrt(z))) / z;
//		}
//		else {
//			double sqrt_neg_z = std::sqrt(-z);
//			return (std::cosh(sqrt_neg_z) - 1.0) / (-z);
//		}
//	}
//
//	double calculate_stumpff_s(double z) {
//		if (std::abs(z) < tolerance) {
//			return 1.0 / 6.0;
//		}
//		else if (z > 0) {
//			double sqrt_z = std::sqrt(z);
//			return (sqrt_z - std::sin(sqrt_z)) / (sqrt_z * sqrt_z * sqrt_z);
//		}
//		else {
//			double sqrt_neg_z = std::sqrt(-z);
//			return (std::sinh(sqrt_neg_z) - sqrt_neg_z) / (sqrt_neg_z * sqrt_neg_z * sqrt_neg_z);
//		}
//	}
//
//	vector_3d calculate_cross_product(const vector_3d& a, const vector_3d& b) {
//		return vector_3d(
//			a.y * b.z - a.z * b.y,
//			a.z * b.x - a.x * b.z,
//			a.x * b.y - a.y * b.x
//		);
//	}
//
//	double calculate_dot_product(const vector_3d& a, const vector_3d& b) {
//		return a.x * b.x + a.y * b.y + a.z * b.z;
//	}
//
//	double calculate_sector_triangle_ratio(const vector_3d& r1, const vector_3d& r2, double dt) {
//		double r1_mag = r1.get_magnitude();
//		double r2_mag = r2.get_magnitude();
//		double sqrt_mu = std::sqrt(gravitational_parameter);
//
//		// Initial guess using Kepler's equation
//		double chord = (r2 - r1).get_magnitude();
//		double s = (r1_mag + r2_mag + chord) / 2.0;
//		double area_triangle = std::sqrt(s * (s - r1_mag) * (s - r2_mag) * (s - chord));
//		double sector_area = area_triangle;
//
//		// Iterative improvement using universal variables
//		for (int i = 0; i < max_iterations; i++) {
//			double previous_sector = sector_area;
//			double z = chord * chord / (r1_mag * r2_mag);
//			double y = sector_area * sector_area / (r1_mag * r2_mag * r2_mag);
//
//			double c = calculate_stumpff_c(z);
//			double s = calculate_stumpff_s(z);
//
//			double new_sector = sqrt_mu * dt / 6.0 - r1_mag * r2_mag * std::sqrt(y / 2.0) *
//				(z * s - (1.0 - z * c));
//
//			if (std::abs(new_sector - previous_sector) < tolerance) {
//				return new_sector / area_triangle;
//			}
//
//			sector_area = new_sector;
//		}
//
//		throw std::runtime_error("Sector-triangle ratio calculation did not converge");
//	}
//
//	std::array<double, 3> solve_gauss_equation(const internal_orbit_data& data) {
//		// Set up the Gaussian equation matrices
//		double d11 = data.tau3 / data.tau;
//		double d21 = -1.0;
//		double d31 = -data.tau1 / data.tau;
//
//		double d12 = -data.tau3;
//		double d22 = data.tau;
//		double d32 = data.tau1;
//
//		// Calculate determinants
//		double det = d11 * (d22 * 1.0 - d32 * 0.0) - d21 * (d12 * 1.0 - d32 * 0.0) + d31 * (d12 * 0.0 - d22 * 0.0);
//
//		// Solve for ranges (ρ values)
//		std::array<double, 3> rho;
//		rho[0] = (data.r[0].get_magnitude() * d11 + data.r[1].get_magnitude() * d21 +
//			data.r[2].get_magnitude() * d31) / det;
//		rho[1] = (data.r[0].get_magnitude() * d12 + data.r[1].get_magnitude() * d22 +
//			data.r[2].get_magnitude() * d32) / det;
//		rho[2] = 0;  // Third equation constraint
//
//		return rho;
//	}
//
//	vector_3d calculate_velocity(const vector_3d& r1, const vector_3d& r2,
//		const vector_3d& r3, double f1, double f3) {
//		vector_3d v2;
//		double r2_mag = r2.get_magnitude();
//
//		// Use Gibbs method for velocity
//		vector_3d n = r1 * (calculate_cross_product(r2, r3).get_magnitude()) +
//			r2 * (calculate_cross_product(r3, r1).get_magnitude()) +
//			r3 * (calculate_cross_product(r1, r2).get_magnitude());
//
//		vector_3d d = calculate_cross_product(r1, r2) +
//			calculate_cross_product(r2, r3) +
//			calculate_cross_product(r3, r1);
//
//		vector_3d s = r1 * ((r2 - r3).get_magnitude()) +
//			r2 * ((r3 - r1).get_magnitude()) +
//			r3 * ((r1 - r2).get_magnitude());
//
//		double v2_mag = std::sqrt(gravitational_parameter / (n.get_magnitude() * d.get_magnitude()));
//		v2 = (calculate_cross_product(d, r2) / r2_mag + s) * v2_mag;
//
//		return v2;
//	}
//
//public:
//	orbit_parameters determine_orbit(double time_1, double time_2, double time_3,
//		const vector_3d& pos_1, const vector_3d& pos_2,
//		const vector_3d& pos_3) {
//		internal_orbit_data data;
//		data.tau = time_3 - time_1;
//		data.tau1 = time_1 - time_2;
//		data.tau3 = time_3 - time_2;
//		data.r = { pos_1, pos_2, pos_3 };
//
//		// Calculate sector/triangle ratios
//		double ratio_12 = calculate_sector_triangle_ratio(pos_1, pos_2, std::abs(data.tau1));
//		double ratio_23 = calculate_sector_triangle_ratio(pos_2, pos_3, std::abs(data.tau3));
//
//		// Solve for ranges using Gauss's equation
//		data.rho = solve_gauss_equation(data);
//
//		// Calculate position vectors
//		std::array<vector_3d, 3> r;
//		for (int i = 0; i < 3; i++) {
//			r[i] = data.r[i] * data.rho[i];
//		}
//
//		// Calculate f and g series coefficients
//		double r2_mag = r[1].get_magnitude();
//		double f1 = 1.0 - ((gravitational_parameter * data.tau1 * data.tau1) / (2.0 * r2_mag * r2_mag * r2_mag));
//		double f3 = 1.0 - ((gravitational_parameter * data.tau3 * data.tau3) / (2.0 * r2_mag * r2_mag * r2_mag));
//
//		// Calculate velocity
//		vector_3d velocity = calculate_velocity(r[0], r[1], r[2], f1, f3);
//
//		// Calculate orbital elements
//		orbit_parameters params;
//		vector_3d h = calculate_cross_product(r[1], velocity);
//		params.angular_momentum = h;
//
//		// Eccentricity vector
//		params.eccentricity_vector = calculate_cross_product(velocity, h) * (1.0 / gravitational_parameter) -
//			r[1] * (1.0 / r[1].get_magnitude());
//		params.eccentricity = params.eccentricity_vector.get_magnitude();
//
//		// Semi-major axis
//		params.semi_major_axis = r2_mag / (2.0 - r2_mag * velocity.get_magnitude() * velocity.get_magnitude() /
//			gravitational_parameter);
//
//		// True anomaly
//		double cos_true_anomaly = calculate_dot_product(params.eccentricity_vector, r[1]) /
//			(params.eccentricity * r2_mag);
//		params.true_anomaly = std::acos(cos_true_anomaly);
//		if (calculate_dot_product(r[1], velocity) < 0) {
//			params.true_anomaly = 2.0 * pi - params.true_anomaly;
//		}
//
//		// Mean motion and period
//		params.mean_motion = std::sqrt(gravitational_parameter /
//			(params.semi_major_axis * params.semi_major_axis * params.semi_major_axis));
//		params.period = 2.0 * pi / params.mean_motion;
//
//		// Eccentric anomaly
//		params.eccentric_anomaly = 2.0 * std::atan(std::sqrt((1.0 - params.eccentricity) /
//			(1.0 + params.eccentricity)) *
//			std::tan(params.true_anomaly / 2.0));
//
//		// Mean anomaly
//		params.mean_anomaly = params.eccentric_anomaly - params.eccentricity * std::sin(params.eccentric_anomaly);
//
//		return params;
//	}
//};

//
//// Define constants
//const double GM = 1.32712440018e20; // Gravitational constant * mass of the sun (m^3/s^2)
//
//// Struct to store a vector
//struct vector_3d {
//	double x, y, z;
//};
//
//// Function to calculate the dot product of two vectors
//double dot(const vector_3d& a, const vector_3d& b) {
//	return a.x * b.x + a.y * b.y + a.z * b.z;
//}
//
//// Function to calculate the cross product of two vectors
//vector_3d cross(const vector_3d& a, const vector_3d& b) {
//	return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
//}
//
//// Function to calculate the norm of a vector
//double norm(const vector_3d& a) {
//	return sqrt(dot(a, a));
//}
//
//// Function to calculate the orbital elements given three position vectors
//void gaussMethod(const vector<vector_3d>& R, const vector<double>& time) {
//	vector_3d r1 = R[0];
//	vector_3d r2 = R[1];
//	vector_3d r3 = R[2];
//
//	double tau1 = time[1] - time[0];
//	double tau3 = time[2] - time[1];
//
//	vector_3d r1xr3 = cross(r1, r3);
//
//	double D0 = dot(r1xr3, r2);
//	double D1 = dot(cross(r2, r3), r2);
//	double D2 = dot(r1xr3, r2);
//	double D3 = dot(cross(r1, r2), r2);
//
//	double A = tau3 / (tau3 - tau1);
//	double B = tau1 / (tau1 - tau3);
//	double E = (D2 * A + D3 * B) / (D0 * (D2 + D3));
//
//	vector_3d R2hat = { r2.x / norm(r2), r2.y / norm(r2), r2.z / norm(r2) };
//	vector_3d H = cross(r2, R2hat);
//
//	double h = norm(H);
//
//	vector_3d eVec = { D2 / (GM * norm(r2)) - r2.x, D2 / (GM * norm(r2)) - r2.y, D2 / (GM * norm(r2)) - r2.z };
//
//	//double a = norm(eVec);
//
//	//double e = sqrt(a*a*a * grav_constant*sun_mass - h*h) / (pow(a, 3.0/2.0)*sqrt(grav_constant * sun_mass));
//
//	double e = norm(eVec);
//
//	double a = pow((h * h / GM) / (1 - e * e), 1.0 / 3.0);
//	double P = 2 * pi * sqrt(a * a * a / GM);
//
//	cout << "Semi-major axis (a): " << a << " m\n";
//	cout << "Eccentricity (e): " << e << "\n";
//	cout << "Orbital period (P): " << P << " s\n";
//
//	swap(a, e);
//
//global_ep.angle = 0;
//global_ep.centerX = 0;
//global_ep.centerY = 0;
//global_ep.semiMajor = a;
//global_ep.semiMinor = a * sqrt(1 - e*e);
//
//}




// This almost works
//
//
//// Constants
//const double MU_SUN = 1.32712440042e20; // Gravitational parameter of the Sun (m^3/s^2)
//const double PI = 3.14159265358979323846;
//
//// Structure to represent a 3D vector
//struct vector_3d {
//	double x, y, z;
//};
//
//// Function to calculate the magnitude of a vector
//double magnitude(const vector_3d& v) {
//	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
//}
//
//// Function to calculate the cross product of two vectors
//vector_3d crossProduct(const vector_3d& v1, const vector_3d& v2) {
//	return { v1.y * v2.z - v1.z * v2.y,
//			v1.z * v2.x - v1.x * v2.z,
//			v1.x * v2.y - v1.y * v2.x };
//}
//
//// Function to calculate the dot product of two vectors
//double dotProduct(const vector_3d& v1, const vector_3d& v2) {
//	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
//}
//
//// Function to calculate the orbital elements using Gauss's method
//void gaussMethod(const vector_3d& r1, const vector_3d& r2, const vector_3d& r3, double t1, double t2, double t3) {
//	// Calculate the velocity vectors
//	//vector_3d v1 = { (r2.x - r1.x) / (t2 - t1), (r2.y - r1.y) / (t2 - t1), 0.0 };
//	//vector_3d v3 = { (r3.x - r2.x) / (t3 - t2), (r3.y - r2.y) / (t3 - t2), 0.0 };
//	vector_3d v1 = { (r2.x - r1.x) / (t2 - t1), (r2.y - r1.y) / (t2 - t1), 0.0 };
//	vector_3d v2 = { (r3.x - r1.x) / (t3 - t1), (r3.y - r1.y) / (t3 - t1), 0.0 };
//	vector_3d v3 = { (r3.x - r2.x) / (t3 - t2), (r3.y - r2.y) / (t3 - t2), 0.0 };
//
//
//	// Calculate the acceleration vector
//	vector_3d a = { (v3.x - v1.x) / (t3 - t1), (v3.y - v1.y) / (t3 - t1), 0.0 };
//
//	// Calculate the specific angular momentum vector
//	vector_3d h = crossProduct(r2, v2);
//
//	// Calculate the semi-major axis
//	double a_orbit = 1.0 / (2.0 / magnitude(r2) - dotProduct(v2, v2) / MU_SUN);
//
//	// Calculate the eccentricity
//	double e = sqrt(1.0 - (magnitude(h) * magnitude(h) / (MU_SUN * a_orbit)));
//
//	// Calculate the inclination (in the equatorial plane, inclination is 0)
//	double i = 0.0;
//
//	// Calculate the longitude of the ascending node (in the equatorial plane, longitude is 0)
//	double Omega = 0.0;
//
//	// Calculate the argument of periapsis
//	double omega = atan2(h.y, h.x);
//
//	// Calculate the mean anomaly
//	double M = sqrt(MU_SUN / (a_orbit * a_orbit * a_orbit)) * (t3 - t1);
//
//	// Print the orbital elements
//	std::cout << "Semi-major axis: " << a_orbit << " m" << std::endl;
//	std::cout << "Eccentricity: " << e << std::endl;
//	std::cout << "Inclination: " << i * 180.0 / PI << " deg" << std::endl;
//	std::cout << "Longitude of the ascending node: " << Omega * 180.0 / PI << " deg" << std::endl;
//	std::cout << "Argument of periapsis: " << omega * 180.0 / PI << " deg" << std::endl;
//	std::cout << "Mean anomaly: " << M * 180.0 / PI << " deg" << std::endl;
//
//
//
//	global_ep.angle = 0;
//	global_ep.centerX = 0;
//	global_ep.centerY = 0;
//	global_ep.semiMajor = a_orbit;
//	global_ep.semiMinor = a_orbit * sqrt(1 - e*e);// ep.semiMinor;
//
//
//}


// Constants
const double MU_SUN = 1.32712440042e20; // Gravitational parameter of the Sun (m^3/s^2)
const double PI = 3.14159265358979323846;

// Structure to represent a 3D vector
struct vector_3d {
	double x, y, z;
};

// Function to calculate the magnitude of a vector
double magnitude(const vector_3d& v) {
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Function to calculate the cross product of two vectors
vector_3d crossProduct(const vector_3d& v1, const vector_3d& v2) {
	return { v1.y * v2.z - v1.z * v2.y,
			v1.z * v2.x - v1.x * v2.z,
			v1.x * v2.y - v1.y * v2.x };
}

// Function to calculate the dot product of two vectors
double dotProduct(const vector_3d& v1, const vector_3d& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// Function to calculate the orbital elements using Gauss's method
void gaussMethod(const vector_3d& r1, const vector_3d& r2, const vector_3d& r3, double t1, double t2, double t3) {
	// Calculate the velocity vectors
	vector_3d v1 = { (r2.x - r1.x) / (t2 - t1), (r2.y - r1.y) / (t2 - t1), 0.0 };
	vector_3d v2 = { (r3.x - r1.x) / (t3 - t1), (r3.y - r1.y) / (t3 - t1), 0.0 };
	vector_3d v3 = { (r3.x - r2.x) / (t3 - t2), (r3.y - r2.y) / (t3 - t2), 0.0 };

	// Calculate the acceleration vector
	vector_3d a = { (v3.x - v1.x) / (t3 - t1), (v3.y - v1.y) / (t3 - t1), 0.0 };

	// Calculate the specific angular momentum vector
	vector_3d h = crossProduct(r2, v2);

	// Calculate the semi-major axis
	double a_orbit = 1.0 / (2.0 / magnitude(r2) - dotProduct(v2, v2) / MU_SUN);

	// Calculate the eccentricity
	double e = sqrt(1.0 - (magnitude(h) * magnitude(h) / (MU_SUN * a_orbit)));

	// Calculate the semi-minor axis
	double b_orbit = a_orbit * sqrt(1.0 - e * e);

	// Calculate the center of the ellipse
	vector_3d center;
	center.x = 0;// = { (r1.x + r2.x + r3.x) / 3.0, (r1.y + r2.y + r3.y) / 3.0, 0.0 };
	center.y = a_orbit * e;
	center.z = 0;

	// Print the orbital elements
	std::cout << "Center of the ellipse: (" << center.x << ", " << center.y << ") m" << std::endl;
	std::cout << "Semi-major axis: " << a_orbit << " m" << std::endl;
	std::cout << "Semi-minor axis: " << b_orbit << " m" << std::endl;
	std::cout << "Eccentricity: " << e << std::endl;

	global_ep.angle = 0;
	global_ep.centerX = center.x;
	global_ep.centerY = center.y;
	global_ep.semiMajor = a_orbit;
	global_ep.semiMinor = b_orbit;



}


int main(int argc, char** argv)
{
	cout << setprecision(20) << endl;

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

double deg_to_rad(double degree)
{
	return degree * (pi / 180.0);
}

double hours_to_seconds(double hours)
{
	return hours * 3600.0;
}

double hours_to_radians(double hours)
{
	return hours * pi / 12.0;
}

struct timestamp_azimuth_data
{
	double timestamp; // seconds 
	double azimuth; // radians
};

struct radius_velocity_data
{
	double angular_velocity;
	double radius;
	double velocity;
};

cartesian_point to_cartesian(double radius, double azimuth)
{
	cartesian_point result;

	result.x = radius * cos(azimuth);
	result.y = radius * sin(azimuth);

	return result;
}

cartesian_point to_spherical(double x, double y)
{
	cartesian_point result;

	result.x = sqrt(x * x + y * y);
	result.y = atan2(y, x);

	return result;
}






vector<cartesian_point> carts;
vector<cartesian_point> orbit_points(5);

vector<cartesian_point> orbit_velocities(5);



void idle_func(void)
{
	static size_t frame_count = 0;

	const double dt = 10000; // 10000 seconds == 2.77777 hours


	static bool calculated_ellipse = false;

	if (true == calculated_ellipse)
	{
		proceed_symplectic4(mercury_pos, mercury_vel, grav_constant, dt);
		positions.push_back(mercury_pos);
	}

	//if (calculated_ellipse == false && positions.size() != 0 && frame_count % 200 == 0)
	//	ellipse_positions.push_back(positions[positions.size() - 1]);

	if (false == calculated_ellipse)
	{
		calculated_ellipse = true;

		// Must have exactly 3 observations
		vector<timestamp_azimuth_data> measurements =
		{
			{hours_to_seconds(0),  deg_to_rad(360)},
			{hours_to_seconds(24), deg_to_rad(359)},
			{hours_to_seconds(48), deg_to_rad(357.98)}
				
			//{hours_to_seconds(0),  deg_to_rad(0) + pi / 2},
			//{hours_to_seconds(24), deg_to_rad(-1) + pi / 2},
			//{hours_to_seconds(52), deg_to_rad(-8) + pi / 2}

			//{hours_to_seconds(0),  deg_to_rad(0) + pi / 2},
			//{hours_to_seconds(24), deg_to_rad(-1) + pi / 2},
			//{hours_to_seconds(48), deg_to_rad(-2.01) + pi / 2}

			//{hours_to_seconds(0),   hours_to_radians(4.79) + pi / 2},
			//{hours_to_seconds(96),  hours_to_radians(4.78) + pi / 2},
			//{hours_to_seconds(192), hours_to_radians(4.77) + pi / 2}

		};





		// Produce 2 radii and velocities
		vector<radius_velocity_data> data_points(2);

		// Constant angular velocity, for example
		//double omega = 4.31e-8; // Ceres average angular velocity
		//double omega_min = 1.99e-7; // Earth average angular velocity

		for (size_t i = 0; i < measurements.size() - 1; i++)
		{
			// Variable angular velocity
			double omega = (measurements[i + 1].azimuth - measurements[i].azimuth) / (measurements[i + 1].timestamp - measurements[i].timestamp);
			double r = cbrt((grav_constant * sun_mass) / (omega * omega));
			double v = omega * r;

			data_points[i].angular_velocity = omega;
			data_points[i].radius = r;
			data_points[i].velocity = v;
		}

		// Produce input data

		double angle0 = measurements[0].azimuth;
		double omega = data_points[0].angular_velocity;
		double r = cbrt((grav_constant * sun_mass) / (omega * omega));
		double v = omega * r;
		radius_velocity_data data_point_0;
		data_point_0.angular_velocity = omega;
		data_point_0.radius = r;
		data_point_0.velocity = v;

		double angle1 = measurements[1].azimuth;
		double r1 = data_points[0].radius;
		double v1 = data_points[0].velocity;
		double a1 = data_points[0].angular_velocity;

		double angle2 = measurements[2].azimuth;
		double r2 = data_points[1].radius;
		double v2 = data_points[1].velocity;
		double a2 = data_points[1].angular_velocity;


		// Convert input data to Cartesian coordinates
		cartesian_point cart0 = to_cartesian(data_point_0.radius, angle0);
		cartesian_point cart1 = to_cartesian(r1, angle1);
		cartesian_point cart2 = to_cartesian(r2, angle2);

		//cartesian_point vel0;
		//vel0.x = (cart1.x - cart0.x) / (measurements[1].timestamp - measurements[0].timestamp);
		//vel0.y = (cart1.y - cart0.y) / (measurements[1].timestamp - measurements[0].timestamp);

		cartesian_point vel1;
		vel1.x = (cart2.x - cart1.x) / (measurements[2].timestamp - measurements[1].timestamp);
		vel1.y = (cart2.y - cart1.y) / (measurements[2].timestamp - measurements[1].timestamp);



		cartesian_point curr_pos = cart1;
		cartesian_point curr_vel = vel1;

		double dt_ = (measurements[2].timestamp - measurements[1].timestamp);

		orbit_points[0] = curr_pos;
		orbit_velocities[0] = curr_vel;

		for (size_t i = 1; i < 5; i++)
		{
			const cartesian_point grav_dir = cartesian_point(curr_pos);
			const double distance = grav_dir.length();

			cartesian_point accel;
			accel.x = -grav_dir.x / distance * (grav_constant * sun_mass / pow(distance, 2.0));
			accel.y = -grav_dir.y / distance * (grav_constant * sun_mass / pow(distance, 2.0));

			curr_vel.x += accel.x * dt_;
			curr_vel.y += accel.y * dt_;

			curr_pos.x += curr_vel.x * dt_;
			curr_pos.y += curr_vel.y * dt_;

			orbit_points[i] = curr_pos;
			orbit_velocities[i] = curr_vel;
		}



		//EllipseParameters ep = fitEllipse(orbit_points, cartesian_point(0, 0));






		//Eigen::MatrixXd A;
		//Eigen::VectorXd b;
		//gatherConstraints(orbit_points, orbit_velocities, A, b);

		//// Solve for ellipse coefficients
		//std::vector<double> coefficients = solveEllipseCoefficients(A, b);

		//// Extract ellipse parameters
		//EllipseParams_min params = extractEllipseParams(coefficients, cartesian_point(0, 0));























		//gauss_orbit_determination orbit_calculator;

		// Example observation times (in seconds)
		//double time_1 = 0;
		//double time_2 = (measurements[1].timestamp - measurements[0].timestamp);
		//double time_3 = (measurements[2].timestamp - measurements[1].timestamp);

		// Example observations in equatorial plane (x,y coordinates, z=0)
		//vector_3d pos_1_obs(orbit_points[0].x, orbit_points[0].y, 0);      // First observation
		//vector_3d pos_2_obs(orbit_points[1].x, orbit_points[1].y, 0);     // Second observation
		//vector_3d pos_3_obs(orbit_points[2].x, orbit_points[2].y, 0);      // Third observation

		// Calculate orbit
		//auto orbit_result = orbit_calculator.determine_orbit(
		//	time_1, time_2, time_3,
		//	pos_1_obs, pos_2_obs, pos_3_obs
		//);

		//// Calculate orbital elements
		//auto elements = orbit_calculator.calculate_orbit_elements(
		//	orbit_result.position,
		//	orbit_result.velocity
		//);

		// Print results
		//std::cout << "Semi-major axis: " << elements.semi_major_axis << " km\n";
		//std::cout << "Eccentricity: " << elements.eccentricity << "\n";
		//std::cout << "True anomaly: " << elements.true_anomaly << " degrees\n";






		//vector3d position = { orbit_points[0].x, orbit_points[0].y, 0.0 }; // Example position vector in km
		//vector3d velocity = { orbit_velocities[0].x, orbit_velocities[0].y, 0.0 }; // Example velocity vector in km/s

		//orbital_elements elements = gauss_method_orbital_elements(position, velocity);
		//print_orbital_elements(elements);





		//vector3d r1_ = { orbit_points[0].x, orbit_points[0].y, 0.0 }; // Example position vector 1 in km
		//vector3d r2_ = { orbit_points[1].x, orbit_points[1].y, 0.0 }; // Example position vector 2 in km
		//vector3d r3_ = { orbit_points[2].x, orbit_points[2].y, 0.0 }; // Example position vector 3 in km

		//double t1_ = dt; // Time at position 1 in seconds
		//double t2_ = 2*dt; // Time at position 2 in seconds
		//double t3_ = 3*dt; // Time at position 3 in seconds

		//orbital_elements elements = gauss_method_with_three_points(r1_, r2_, r3_, t1_, t2_, t3_);
		//print_orbital_elements(elements);




		//gauss_orbit_calculator calculator;

		//// Example times in seconds from epoch
		//double time_1 = 0;
		//double time_2 = (measurements[1].timestamp - measurements[0].timestamp);
		//double time_3 = (measurements[2].timestamp - measurements[1].timestamp);

		//// Example position vectors in equatorial plane (x,y,z) in kilometers
		//vector_3d pos_1(cart0.x, cart0.y, 0.0);
		//vector_3d pos_2(cart1.x, cart1.y, 0.0);
		//vector_3d pos_3(cart2.x, cart2.y, 0.0);

		//try {
		//	orbit_parameters orbit = calculator.determine_orbit(
		//		time_1, time_2, time_3,
		//		pos_1, pos_2, pos_3
		//	);

		//	std::cout << "Orbital Elements:\n";
		//	std::cout << "Semi-major axis (km): " << orbit.semi_major_axis << "\n";
		//	std::cout << "Eccentricity: " << orbit.eccentricity << "\n";
		//	std::cout << "True anomaly (rad): " << orbit.true_anomaly << "\n";
		//	std::cout << "Period (s): " << orbit.period << "\n";
		//	std::cout << "Mean motion (rad/s): " << orbit.mean_motion << "\n";
		//	std::cout << "Mean anomaly (rad): " << orbit.mean_anomaly << "\n";
		//	std::cout << "Eccentric anomaly (rad): " << orbit.eccentric_anomaly << "\n";

		//	std::cout << "\nAngular Momentum Vector (km²/s):\n";
		//	std::cout << "x: " << orbit.angular_momentum.x << "\n";
		//	std::cout << "y: " << orbit.angular_momentum.y << "\n";
		//	std::cout << "z: " << orbit.angular_momentum.z << "\n";

		//	std::cout << "\nEccentricity Vector:\n";
		//	std::cout << "x: " << orbit.eccentricity_vector.x << "\n";
		//	std::cout << "y: " << orbit.eccentricity_vector.y << "\n";
		//	std::cout << "z: " << orbit.eccentricity_vector.z << "\n";
		//}
		//catch (const std::runtime_error& e) {
		//	std::cerr << "Error: " << e.what() << std::endl;
		//	//return;
		//}


	//	std::vector<observation> observations = {
	//{0.0,  {orbit_points[0].x, orbit_points[0].y}},  // Example: Satellite at x=1000km, y=0 at t=0
	//{dt,   {orbit_points[1].x, orbit_points[1].y}}, // 1 minute later at x=2000km, y=0
	//{2*dt, {orbit_points[2].x, orbit_points[2].y}} // 2 minutes later at x=3000km, y=0
	//	};

	//	vector_2d velocity = gauss_method_orbit_determination(observations);
	//	double semi_major_axis, eccentricity;

	//	compute_orbit_parameters(observations[1].position, velocity, semi_major_axis, eccentricity);

	//	std::cout << "Semi-major axis: " << semi_major_axis << " km" << std::endl;
	//	std::cout << "Eccentricity: " << eccentricity << std::endl;




		//ns::Vector3 r1_ = { orbit_points[0].x, orbit_points[0].y, 0 }; // km
		//ns::Vector3 r2_ = { orbit_points[1].x, orbit_points[1].y, 0 }; // km
		//ns::Vector3 r3_ = { orbit_points[2].x, orbit_points[2].y, 0 }; // km
		 
		//double dt21_ = dt; // seconds
		//double dt31_ = 2*dt; // seconds

		//ns::Vector3 v2_ = solve_gauss(r1_, r2_, r3_, dt21_, dt31_);

		//ns::Vector3 r2_(orbit_points[0].x, orbit_points[0].y);
		//ns::Vector3 v2_(orbit_velocities[0].x, orbit_velocities[0].y);

		//OrbitElements orbit = calc_orbit_elements(r2_, v2_);

		//std::cout << "Semi-major Axis: " << orbit.semi_major_axis << " km" << std::endl;
		//std::cout << "Eccentricity: " << orbit.eccentricity << std::endl;
		//std::cout << "Inclination: " << orbit.inclination * 180 / pi << " degrees" << std::endl;
		//std::cout << "RAAN: " << orbit.raan * 180 / pi << " degrees" << std::endl;
		//std::cout << "Argument of Periapsis: " << orbit.arg_periapsis * 180 / pi << " degrees" << std::endl;
		//std::cout << "True Anomaly: " << orbit.true_anomaly * 180 / pi << " degrees" << std::endl;


		//global_ep.angle = orbit.true_anomaly;
		//global_ep.centerX = 0;
		//global_ep.centerY = 0;
		//global_ep.semiMajor = orbit.semi_major_axis;
		//global_ep.semiMinor = orbit.semi_major_axis * sqrt(1 - orbit.eccentricity*orbit.eccentricity);// ep.semiMinor;









		//cout << global_ep.angle << endl;
		//cout << global_ep.centerX << endl;
		//cout << global_ep.centerY << endl;
		//cout << global_ep.semiMajor << endl;
		//cout << global_ep.semiMinor << endl;


		//gauss_orbit_determination orbit_calculator;

		//// Example observation times (in seconds)
		//double time_1 = 0;
		//double time_2 = dt;    // 1 hour later
		//double time_3 = 2*dt;    // 2 hours later

		//// Example observations in equatorial plane (x,y coordinates, z=0)
		//vector_3d pos_1_obs(orbit_points[0].x, orbit_points[0].y, 0);      // First observation
		//vector_3d pos_2_obs(orbit_points[1].x, orbit_points[1].y, 0);     // Second observation
		//vector_3d pos_3_obs(orbit_points[2].x, orbit_points[2].y, 0);     // Third observation

		//// Calculate orbit
		//auto orbit_result = orbit_calculator.determine_orbit(
		//	time_1, time_2, time_3,
		//	pos_1_obs, pos_2_obs, pos_3_obs
		//);

		//// Calculate orbital elements
		//auto elements = orbit_calculator.calculate_orbit_elements(
		//	orbit_result.position,
		//	orbit_result.velocity
		//);

		//// Print results
		//std::cout << "Semi-major axis: " << elements.semi_major_axis << " km\n";
		//std::cout << "Eccentricity: " << elements.eccentricity << "\n";
		//std::cout << "True anomaly: " << elements.true_anomaly << " degrees\n";




		//VectorXd params = solveEllipseParameters(orbit_points, orbit_velocities, cartesian_point(0, 0));

		//double h = params[0], k = params[1], a = params[2], b = params[3], angle = params[4];

		//global_ep.angle = angle;
		//global_ep.centerX = 0;
		//global_ep.centerY = k;
		//global_ep.semiMajor = a;
		//global_ep.semiMinor = b;

		//cout << global_ep.angle << endl;
		//cout << global_ep.centerX << endl;
		//cout << global_ep.centerY << endl;
		//cout << global_ep.semiMajor << endl;
		//cout << global_ep.semiMinor << endl;

//
//vector_3d r1_ = { orbit_points[0].x, orbit_points[0].y, 0.0};  // Example position at t1
//vector_3d r2_ = { orbit_points[1].x, orbit_points[1].y, 0.0 };  // Example position at t2
//vector_3d r3_ = { orbit_points[2].x, orbit_points[2].y, 0.0 };  // Example position at t3
//
//double t1_ = 0.0;
//double t2_ = dt; 
//double t3_ = 2 * dt;




		//// Example usage with equatorial plane motion
		//vector3d r1_(orbit_points[0].x, orbit_points[0].y, 0.0); // First position vector
		//vector3d r2_(orbit_points[1].x, orbit_points[1].y, 0.0); // Second position vector
		//vector3d r3_(orbit_points[2].x, orbit_points[2].y, 0.0); // Third position vector

		//double t1 = dt; // Time of first observation
		//double t2 = dt*2.0; // Time of second observation
		//double t3 = dt*3.0; // Time of third observation
		//double mu = grav_constant*sun_mass; // Gravitational parameter (adjust according to the system)

		//determine_orbit_parameters(r1_, r2_, r3_, t1, t2, t3, mu);


		//calculate ellipse semi - major axis, eccentricity, angle, and center from 5 input points and 5 input velocities and a focus point. using c++. the ellipse must not be centered at the origin

//
//std::vector<observation> observations = {
//	{{cart0.x, cart0.y, 0}, 0},
//	{{cart1.x, cart1.y, 0}, (measurements[1].timestamp - measurements[0].timestamp)},
//	{{cart2.x, cart2.y, 0}, (measurements[2].timestamp - measurements[1].timestamp) + (measurements[1].timestamp - measurements[0].timestamp)}
//};



//double time_1 = 0;
//double time_2 = (measurements[1].timestamp - measurements[0].timestamp);
//double time_3 = (measurements[2].timestamp - measurements[1].timestamp);

//// Example position vectors in equatorial plane (x,y,z) in kilometers
//vector_3d pos_1(cart0.x, cart0.y, 0.0);
//vector_3d pos_2(cart1.x, cart1.y, 0.0);
//vector_3d pos_3(cart2.x, cart2.y, 0.0);

//
//orbital_elements elements = determine_orbit(observations);
//
//std::cout << "Semi-major axis: " << elements.semi_major_axis << " km\n";
//std::cout << "Eccentricity: " << elements.eccentricity << "\n";
//std::cout << "Inclination: " << elements.inclination << " rad\n";
//std::cout << "Longitude of ascending node: " << elements.longitude_of_ascending_node << " rad\n";
//std::cout << "Argument of periapsis: " << elements.argument_of_periapsis << " rad\n";
//std::cout << "Mean anomaly at epoch: " << elements.mean_anomaly_at_epoch << " rad\n";




//vector_3d r1_ = { cart0.x, cart0.y, 0.1 };
//vector_3d r2_ = { cart1.x, cart1.y, 0.2 };
//vector_3d r3_ = { cart2.x, cart2.y, 0.3 };
//
//// Define the time points
//double t1_ = 0.0;
//double t2_ = (measurements[1].timestamp - measurements[0].timestamp);
//double t3_ = (measurements[2].timestamp - measurements[1].timestamp) + (measurements[1].timestamp - measurements[0].timestamp);
//
//// Calculate the orbital elements
//OrbitalElements elements = gaussMethod(r1_, r2_, r3_, t1_, t2_, t3_);
//
//
//
//cout << "gauss: " << elements.semiMajorAxis << " " << elements.eccentricity << endl;

//
//global_ep.angle = elements.trueAnomaly;
//global_ep.centerX = 0;
//global_ep.centerY = 0;
//global_ep.semiMajor = elements.semiMajorAxis;
//global_ep.semiMinor = elements.semiMajorAxis * sqrt(1 - elements.eccentricity * elements.eccentricity);// ep.semiMinor;


//Vector3d pos1(orbit_points[0].x, orbit_points[0].y, 0.0);       // Position at time 1 (m)
//Vector3d pos2(orbit_points[1].x, orbit_points[1].y, 0.0);       // Position at time 2 (m)
//Vector3d pos3(orbit_points[2].x, orbit_points[2].y, 0.0);      // Position at time 3 (m)
//
//double t1 = 0.0;                        // Time 1 (seconds)
//double t2 = dt;// (measurements[1].timestamp - measurements[0].timestamp);       // Time 2 (2 days in seconds)
//double t3 = 2 * dt;// (measurements[2].timestamp - measurements[1].timestamp) + (measurements[1].timestamp - measurements[0].timestamp);       // Time 3 (4 days in seconds)
//
//GaussOrbitDetermination god(pos1, pos2, pos3, t1, t2, t3);
//god.solve();

//
//
//
//	// Example heliocentric positions in AU (converted to meters for consistency with GM)
//Vector3D r1_ = { cart0.x, cart0.y, 0.0 };
//Vector3D r2_ = { cart1.x, cart1.y, 0.0 };
//Vector3D r3_ = { cart2.x, cart2.y, 0.0 };
//
//// Times in seconds from some epoch (e.g., J2000)
//double t1_ = 0.0;
//double t2_ = (measurements[1].timestamp - measurements[0].timestamp);
//double t3_ = (measurements[2].timestamp - measurements[1].timestamp) + (measurements[1].timestamp - measurements[0].timestamp);
//
//gaussMethod(r1_, r2_, r3_, t1_, t2_, t3_);



//vector_3d r1_ = { cart0.x, cart0.y, 0.0 };
//vector_3d r2_ = { cart1.x, cart1.y, 0.0 };
//vector_3d r3_ = { cart2.x, cart2.y, 0.0 };
//
//double t1_ = 0.0;
//double t2_ = hours_to_seconds(24);
//double t3_ = hours_to_seconds(48);// 2 * dt;
//
//
//calculateOrbitalElements(r1_, r2_, r3_, t1_, t2_, t3_);



	// Define three position vectors (in meters) and their corresponding times (in seconds)
//vector<vector_3d> R = { { orbit_points[0].x, orbit_points[0].y, 0.0 }, { orbit_points[1].x, orbit_points[1].y, 0.0 }, { orbit_points[2].x, orbit_points[2].y, 0.0 } };
//vector<double> time = { 0,  dt,  2*dt };
//
//gaussMethod(R, time);
//
//
// 
// 
// 
// This almost works:
//vector_3d r1_ = { orbit_points[0].x, orbit_points[0].y, 0.0 };
//vector_3d r2_ = { orbit_points[1].x, orbit_points[1].y, 0.0 };
//vector_3d r3_ = { orbit_points[2].x, orbit_points[2].y, 0.0 };
//
//double t1_ = 0.0;
//double t2_ = dt_;// hours_to_seconds(24);
//double t3_ = 2 * dt_;// hours_to_seconds(48);// 2 * dt;
//
//gaussMethod(r1_, r2_, r3_, t1_, t2_, t3_);





vector_3d r1_ = { orbit_points[0].x, orbit_points[0].y, 0.0 };
vector_3d r2_ = { orbit_points[1].x, orbit_points[1].y, 0.0 };
vector_3d r3_ = { orbit_points[2].x, orbit_points[2].y, 0.0 };

double t1_ = 0.0;
double t2_ = dt_;// hours_to_seconds(24);
double t3_ = 2 * dt_;// hours_to_seconds(48);// 2 * dt;

gaussMethod(r1_, r2_, r3_, t1_, t2_, t3_);











		mercury_pos.x = orbit_points[0].x;
		mercury_pos.y = orbit_points[0].y;
		mercury_vel.x = orbit_velocities[0].x;
		mercury_vel.y = orbit_velocities[0].y;

		carts.push_back(cart1);
		carts.push_back(cart2);
		//carts.push_back(cart0);
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

	glColor3f(1.0, 1.0, 1.0);


	glVertex3d(sun_pos.x, sun_pos.y, sun_pos.z);

	if (carts.size() > 0)
	{
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(carts[0].x, carts[0].y, 0);

		glColor3f(0.0, 1.0, 0.0);
		glVertex3d(carts[1].x, carts[1].y, 0);
	}

	//glColor3f(1.0, 0.0, 1.0);


	for (size_t i = 0; i < orbit_points.size(); i++)
		glVertex3d(orbit_points[i].x, orbit_points[i].y, 0);


	for (size_t i = 0; i < positions.size(); i++)
		glVertex3d(positions[i].x, positions[i].y, 0);

	glEnd();







	glPushMatrix();

	glColor3f(1.0, 0.5, 0.0);

	glTranslated(global_ep.centerX, global_ep.centerY, 0);
	glRotated(global_ep.angle / (2 * pi) * 360.0, 0, 0, 1);
	glTranslated(-global_ep.centerX, -global_ep.centerY, 0);

	DrawEllipse(global_ep.centerX, global_ep.centerY, global_ep.semiMinor, global_ep.semiMajor, 100);




	glPopMatrix();






	//glBegin(GL_POINTS);

	//glColor3f(1.0, 0.0, 0.0);

	//for (size_t i = 0; i < double_check_ellipse_points.size(); i++)
	//	glVertex3d(double_check_ellipse_points[i][0], double_check_ellipse_points[i][1], double_check_ellipse_points[i][2]);

	//glEnd();



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

	main_camera.Set();
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}


