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
    double h = params[0], k = params[1], a = params[2], b = params[3];
    double error = 0;

	EllipseParameters2 ep;
	ep.angle = 0; 
	ep.center(0) = 0;
	ep.center(1) = k;
	ep.semi_major_axis = a;
	ep.eccentricity = sqrt(1 - (b * b) / (a * a));

	Eigen::Vector2d focus1;
	focus1(0) = focus.x;
	focus1(1) = focus.y;

	Eigen::Vector2d focus2 = calculateSecondFocus(ep, focus1);
 
	for(size_t i = 0; i < points.size(); i++)
	{
		cartesian_point p = points[i];
		cartesian_point v = velocities[i];

		double dist1 = distance(p.x, p.y,  focus1(0),  focus1(1));
		double dist2 = distance(p.x, p.y,  focus2(0),  focus2(1));

		error += pow(dist1 + dist2 - 2 * a, 2);

        // Since we're axis-aligned, we simplify velocity condition:
        // Velocity should be more in line with the axis of the ellipse
        double velError = 0;

        if (abs(v.x) > abs(v.y)) // Suggesting a is along x
			velError = pow(v.y / v.x - (k - p.y) / (h - p.x), 2); // Check alignment with y
		else
			velError = pow(v.x / v.y - (h - p.x) / (k - p.y), 2); // Check alignment with x

        error += velError;
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

	const double d = (mvec[4] - mvec[0]) / mvec[4];
	
	VectorXd params(4); // h, k, a, b

	cout << "d: " << d << endl;

	if(d < 0.1)
	params << 1, 1, m, m; // Initial guess
	else
	params << 1, 1, m*0.125, m*0.125; // Initial guess


	int iterations = 100000;
    double stepSize = 0.0001;

    for (int i = 0; i < iterations; i++) 
	{
        VectorXd gradient = VectorXd::Zero(4);

        for (int j = 0; j < 4; j++) 
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

	//const double dt = 10000; // 10000 seconds == 2.77777 hours


	static bool calculated_ellipse = false;

	//if (true == calculated_ellipse)
	//{
	//	proceed_symplectic4(mercury_pos, mercury_vel, grav_constant, dt);
	//	positions.push_back(mercury_pos);
	//}

	//if (calculated_ellipse == false && positions.size() != 0 && frame_count % 200 == 0)
	//	ellipse_positions.push_back(positions[positions.size() - 1]);

	if (false == calculated_ellipse)
	{
		calculated_ellipse = true;

		// Must have exactly 3 observations
		vector<timestamp_azimuth_data> measurements =
		{
			{hours_to_seconds(0),  deg_to_rad(360) + pi / 2},
			{hours_to_seconds(24), deg_to_rad(359) + pi / 2},
			{hours_to_seconds(48), deg_to_rad(357.95) + pi / 2}

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
		//cartesian_point cart0 = to_cartesian(data_point_0.radius, angle0);
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

		double dt = measurements[2].timestamp - measurements[1].timestamp;

		orbit_points[0] = curr_pos;
		orbit_velocities[0] = curr_vel;

		for (size_t i = 1; i < 5; i++)
		{
			const cartesian_point grav_dir = cartesian_point(curr_pos);
			const double distance = grav_dir.length();

			cartesian_point accel;
			accel.x = -grav_dir.x / distance * (grav_constant * sun_mass / pow(distance, 2.0));
			accel.y = -grav_dir.y / distance * (grav_constant * sun_mass / pow(distance, 2.0));

			curr_vel.x += accel.x * dt;
			curr_vel.y += accel.y * dt;

			curr_pos.x += curr_vel.x * dt;
			curr_pos.y += curr_vel.y * dt;

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











		VectorXd params = solveEllipseParameters(orbit_points, orbit_velocities, cartesian_point(0, 0));

		double h = params[0], k = params[1], a = params[2], b = params[3];

		global_ep.angle = 0;
		global_ep.centerX = 0;
		global_ep.centerY = k;
		global_ep.semiMajor = a;
		global_ep.semiMinor = b;

		cout << global_ep.angle << endl;
		cout << global_ep.centerX << endl;
		cout << global_ep.centerY << endl;
		cout << global_ep.semiMajor << endl;
		cout << global_ep.semiMinor << endl;

	










		//std::vector<Eigen::Vector2d> points = {
		//	Eigen::Vector2d(orbit_points[0].x, orbit_points[0].y),
		//	Eigen::Vector2d(orbit_points[1].x, orbit_points[1].y),
		//	Eigen::Vector2d(orbit_points[2].x, orbit_points[2].y),
		//	Eigen::Vector2d(orbit_points[3].x, orbit_points[3].y),
		//	Eigen::Vector2d(orbit_points[4].x, orbit_points[4].y)
		//};

		//Eigen::Vector2d focus(0, 0);

		//EllipseFitter fitter;

		//auto params1 = fitter.fitEllipse(points, focus, ErrorMetric::ALGEBRAIC);
		////auto params1 = fitter.fitEllipse(points, focus, ErrorMetric::GEOMETRIC);
		////auto params1 = fitter.fitEllipse(points, focus, ErrorMetric::FOCAL_SUM);
		////auto params4 = fitter.fitEllipse(points, focus); // Uses HYBRID by default


		//global_ep.angle = params1.angle;
		//global_ep.centerX = params1.center.x();
		//global_ep.centerY = params1.center.y();
		//global_ep.semiMajor = params1.semi_major_axis;
		//global_ep.semiMinor = params1.semi_major_axis * sqrt(1 - params1.eccentricity*params1.eccentricity);

		//cout << global_ep.angle << endl;
		//cout << global_ep.centerX << endl;
		//cout << global_ep.centerY << endl;
		//cout << global_ep.semiMajor << endl;
		//cout << global_ep.semiMinor << endl;









		//calculate ellipse semi - major axis, eccentricity, angle, and center from 5 input points and 5 input velocities and a focus point. using c++. the ellipse must not be centered at the origin



		//mercury_pos.x = cart1.x;
		//mercury_pos.y = cart1.y;
		//mercury_vel.x = vel1.x;
		//mercury_vel.y = vel1.y;

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

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}


