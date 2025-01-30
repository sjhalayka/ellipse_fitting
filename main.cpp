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


class EllipseParameters
{
public:
	double centerX = 0;
	double centerY = 0;
	double semiMajor = 0;
	double semiMinor = 0;
	double angle = 0;
};


// Used to draw things in the OpenGL code
EllipseParameters global_ep;









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


template<size_t N>
class Vector_nD 
{
public:
	std::array<double, N> components;

	// Helper function to get the sign of permutation
	static int permutationSign(const std::array<int, (N - 1) > &perm) 
	{
		int inversions = 0;
		for (int i = 0; i < (N - 2); i++) 
		{
			for (int j = i + 1; j < (N - 1); j++) 
			{
				if (perm[i] > perm[j]) inversions++;
			}
		}

		return (inversions % 2 == 0) ? 1 : -1;
	}

public:
	Vector_nD(const std::array<double, N>& comps) : components(comps) {}

	Vector_nD(void) 
	{
		for (size_t i = 0; i < N; i++)
			components[i] = 0;
	}


	double operator[](size_t index) const {
		if (index >= N) throw std::out_of_range("Index out of bounds");
		return components[index];
	}

	static Vector_nD cross_product(const std::vector<Vector_nD<N>>& vectors)
	{
		if (vectors.size() != (N - 1))
		{
			throw std::invalid_argument("nD cross product requires exactly (n - 1) vectors");
		}

		std::array<double, N> result;// = { 0, 0, 0, 0, 0, 0, 0 };

		for (size_t i = 0; i < N; i++)
			result[i] = 0.0;

		// These are the indices we'll use for each component calculation
		std::array<int, (N - 1)> baseIndices;// = { 0, 1, 2, 3, 4, 5 };

		for (int i = 0; i < N - 1; i++)
			baseIndices[i] = i;

		// For each component of the result vector
		for (int k = 0; k < N; k++) {
			// Skip k in our calculations - this is equivalent to removing the k-th column
			// For each permutation of the remaining 6 indices
			do {
				// Calculate sign of this term
				int sign = permutationSign(baseIndices);

				// Calculate the product for this permutation
				double product = 1.0;
				for (int i = 0; i < (N - 1); i++) {
					int col = baseIndices[i];
					// Adjust column index if it's past k
					int actualCol = (col >= k) ? col + 1 : col;
					product *= vectors[i][actualCol];
				}

				result[k] += sign * product;

			} while (std::next_permutation(baseIndices.begin(), baseIndices.end()));

			// Reset indices for next component
			for (int i = 0; i < N - 1; i++)
				baseIndices[i] = i;
		}

		return Vector_nD(result);
	}

	static double dot_product(const Vector_nD<N>& a, const Vector_nD<N>&b)
	{
		double dot_prod = 0;

		for (size_t i = 0; i < N; i++)
			dot_prod += a[i]*b[i];

		return dot_prod;
	}

	// Method to print vector (for debugging)
	void print() const {
		std::cout << "(";
		for (int i = 0; i < N; i++) {
			std::cout << components[i];
			if (i < (N - 1)) std::cout << ", ";
		}
		std::cout << ")" << std::endl;
	}
};


template <typename size_t N>
double determinant_nxn(const MatrixXd& m)
{
	if (m.cols() != m.rows())
	{
		cout << "Matrix must be square" << endl;
		return 0;
	}

	Vector_nD<N> a_nd;
	
	for (size_t i = 0; i < N; i++)
		a_nd.components[i] = m(0, i);

	std::vector<Vector_nD<N>> vectors;

	for (size_t i = 1; i < N; i++)
	{
		Vector_nD<N> non;

		for (size_t j = 0; j < N; j++)
			non.components[j] = m(i, j);

		vectors.push_back(non);
	}

	// Compute the cross product using (N - 1) vectors
	Vector_nD<N> result = Vector_nD<N>::cross_product(vectors);

	// Flip handedness
	for (size_t i = 0; i < result.components.size(); i++)
		if (i % 2 == 1)
			result.components[i] = -result.components[i];

	double d_dot = Vector_nD<N>::dot_product(a_nd, result);
	
	cout << d_dot << endl;
	cout << m.determinant() << endl << endl;

	return d_dot;
}





// Function to replace a column in a matrix
Matrix3d replace_column(
	const Matrix3d& m,
	const Vector3d& vector,
	int col)
{
	Matrix3d replacedMatrix = m;

	for (size_t i = 0; i < 3; i++)
		replacedMatrix(i, col) = vector(i);

	return replacedMatrix;
}

EllipseParameters extract_ellipse_parameters3x3(const Eigen::VectorXd& coefficients)
{
	//A(i, 0) = x * x;
	//A(i, 1) = 0;// x* y;
	//A(i, 2) = x * x;
	//A(i, 3) = 0;// x;
	//A(i, 4) = y;

	double a = coefficients(0);
	double b = 0;
	double c = coefficients(1);
	double d = 0;
	double e = coefficients(2);
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

	if (a2 > c2)
	{
		cout << "Rotating ellipse" << endl;
		std::swap(semiMajor, semiMinor);
		theta += pi / 2;
	}

	EllipseParameters params;
	params.centerX = centerX;
	params.centerY = centerY;
	params.semiMajor = semiMinor; // swapped
	params.semiMinor = semiMajor;
	params.angle = theta;

	return params;
}






EllipseParameters fit_ellipse3x3(const std::vector<cartesian_point>& points)
{
	if (points.size() != 3) {
		std::cerr << "Error: 3 points are required.\n";
		return EllipseParameters();
	}

	Eigen::MatrixXd A(3, 3);
	Eigen::VectorXd b(3);

	// Fill the matrix A and vector b with the equations from the points
	for (size_t i = 0; i < 3; i++)
	{
		double x = points[i].x;
		double y = points[i].y;

		A(i, 0) = x * x;
		A(i, 1) = y * y;
		A(i, 2) = y;

		b(i) = -1;
	}

	Eigen::VectorXd ellipse_coefficients_Jacobi = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

	//cout << ellipse_coefficients_Jacobi << endl;



	// Solve for the ellipse parameters
	double detMatrix = determinant_nxn<3>(A);
	EllipseParameters ep;

	if (detMatrix == 0)
	{
		cout << "no determinant" << endl;
	}
	else
	{
		double detMatrixa = determinant_nxn<3>(replace_column(A, b, 0));
		double detMatrixc = determinant_nxn<3>(replace_column(A, b, 1));
		double detMatrixe = determinant_nxn<3>(replace_column(A, b, 2));

		double a = detMatrixa / detMatrix;
		double c = detMatrixc / detMatrix;
		double e = detMatrixe / detMatrix;

		//std::cout << "Solution:" << std::endl;
		//std::cout << "a = " << a << std::endl;
		//std::cout << "c = " << c << std::endl;
		//std::cout << "e = " << e << std::endl;

		Eigen::Vector3d ellipse_coefficients;
		ellipse_coefficients(0) = a;
		ellipse_coefficients(1) = c;
		ellipse_coefficients(2) = e;

		ep = extract_ellipse_parameters3x3(ellipse_coefficients);
	}

	global_ep.angle = ep.angle;
	global_ep.centerX = ep.centerX;
	global_ep.centerY = ep.centerY;
	global_ep.semiMajor = ep.semiMajor;
	global_ep.semiMinor = ep.semiMinor;

	return ep;
}















// Structure to represent a 3D vector
class vector_3d {
public:
	double x, y, z;

	vector_3d operator+(const vector_3d& rhs) const
	{
		vector_3d ret;
		ret.x = x + rhs.x;
		ret.y = y + rhs.y;
		ret.z = z + rhs.z;

		return ret;
	}

	vector_3d operator*(const double& rhs) const
	{
		vector_3d ret;
		ret.x = x * rhs;
		ret.y = y * rhs;
		ret.z = z * rhs;

		return ret;
	}

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
void gauss_method(const vector_3d& r1, const vector_3d& r2, const vector_3d& r3, double t1, double t2, double t3)
{
	// Calculate average velocity
	const vector_3d v21 = { (r2.x - r1.x) / (t2 - t1), (r2.y - r1.y) / (t2 - t1), 0.0 };
	const vector_3d v32 = { (r3.x - r2.x) / (t3 - t2), (r3.y - r2.y) / (t3 - t2), 0.0 };
	const vector_3d v_avg = (v21 + v32) * 0.5;

	// Calculate average position
	const vector_3d r_avg = (r1 + r2 + r3) * (1.0 / 3.0);

	// Calculate the specific angular momentum vector
	const vector_3d h = crossProduct(r_avg, v_avg);

	// Calculate the semi-major axis
	double a_orbit = 1.0 / (2.0 / magnitude(r_avg) - dotProduct(v_avg, v_avg) / (sun_mass * grav_constant));

	// Calculate the eccentricity
	const double e = sqrt(1.0 - (magnitude(h) * magnitude(h) / ((sun_mass * grav_constant) * a_orbit)));

	// Calculate the semi-minor axis
	double b_orbit = a_orbit * sqrt(1.0 - e * e);



	// Calculate the center of the ellipse
	vector_3d center;
	center.x = 0;
	center.z = 0;

	const double dx = (magnitude(r2) - magnitude(r1)) + (magnitude(r3) - magnitude(r2));

	if (dx < 0)
	{
		//		swap(a_orbit, b_orbit);
		center.y = -a_orbit * e;
	}
	else
	{

		center.y = a_orbit * e;
	}

	global_ep.angle = 0;
	global_ep.centerX = center.x;
	global_ep.centerY = center.y;
	global_ep.semiMajor = a_orbit;
	global_ep.semiMinor = b_orbit;
}


int main(int argc, char** argv)
{
	cout << std::scientific << endl;
	//	cout << setprecision(20) << endl;

	//const size_t N = 10;

	//MatrixXd m(N, N);

	//for (size_t i = 0; i < N; i++)
	//	for (size_t j = 0; j < N; j++)
	//		m(i, j) = static_cast<double>(rand());

	//determinant_nxn<N>(m);

	//return 0;

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
	static double const cr2 = cbrt(2.0);

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




vector<cartesian_point> orbit_points;
vector<cartesian_point> orbit_velocities;




double calculate_orbit_radius(double omega, double omegaDot, double GM)
{
	// Check for valid angular velocity
	if (omega <= 0)
		return 0;

	// Calculate the radius using Newton's dynamics
	double radius = std::cbrt(GM / (omega * omega)); // Approximation for circular part

	// Adjust for angular velocity derivative
	double velocityAdjustment = -omegaDot / (2 * omega);
	double correctedRadius = radius + velocityAdjustment * radius; // Approximate correction

	return correctedRadius;
}










void idle_func(void)
{
	static size_t frame_count = 0;

	// For numerical integration
	// This variable is filled out upon finding the ellipse
	static double dt = 0;

	static bool calculated_ellipse = false;

	if (true == calculated_ellipse)
	{
		proceed_symplectic4(orbiter_pos, orbiter_vel, grav_constant, dt);
		positions.push_back(orbiter_pos);
	}

	if (false == calculated_ellipse)
	{
		calculated_ellipse = true;

		// Must have exactly 3 observations to produce the input data
		vector<timestamp_azimuth_data> measurements =
		{
			//{hours_to_seconds(0),  deg_to_rad(0)},
			//{hours_to_seconds(24), deg_to_rad(1)},
			//{hours_to_seconds(48), deg_to_rad(2.01)}

			{hours_to_seconds(0),  deg_to_rad(0)},
			{hours_to_seconds(24), deg_to_rad(1)},
			{hours_to_seconds(48), deg_to_rad(1.99)}
		};

		double prev_omega = (measurements[1].azimuth - measurements[0].azimuth) / (measurements[1].timestamp - measurements[0].timestamp);

		// Produce 3 angular accelerations. Only the last one is actually used,
		// so just use a quick dummy value of zero, for this one anyway
		vector<double> d_omega_data;
		d_omega_data.push_back(0);

		for (size_t i = 0; i < measurements.size() - 1; i++)
		{
			// Variable angular velocity
			const double omega = (measurements[i + 1].azimuth - measurements[i].azimuth) / (measurements[i + 1].timestamp - measurements[i].timestamp);

			// Variable angular acceleration
			const double d_omega = (omega - prev_omega) / (measurements[i + 1].timestamp - measurements[i].timestamp);
			prev_omega = omega;

			d_omega_data.push_back(d_omega);
		}

		// Use the last d_omega element to serve as the constant angular acceleration
		// We use constant angular acceleration because we only have 3 observations
		// If we wanted variable angular acceleration, we would need 4 observations instead
		const double constant_angular_acceleration = d_omega_data[2];

		// Produce 3 radii. The first one isn't actually used,
		// so just use a quick dummy value of zero
		vector<double> radii_data;
		radii_data.push_back(0);

		for (size_t i = 0; i < measurements.size() - 1; i++)
		{
			// Variable angular velocity, constant angular acceleration
			const double omega = (measurements[i + 1].azimuth - measurements[i].azimuth) / (measurements[i + 1].timestamp - measurements[i].timestamp);
			const double radius = calculate_orbit_radius(omega, constant_angular_acceleration, sun_mass * grav_constant);
			radii_data.push_back(radius);
		}

		// Produce input data
		const double angle1 = measurements[1].azimuth;
		const double r1 = radii_data[1];

		const double angle2 = measurements[2].azimuth;
		const double r2 = radii_data[2];

		// Convert input data to Cartesian coordinates
		cartesian_point cart1 = to_cartesian(r1, angle1);
		cartesian_point cart2 = to_cartesian(r2, angle2);

		cartesian_point vel1;
		vel1.x = (cart2.x - cart1.x) / (measurements[2].timestamp - measurements[1].timestamp);
		vel1.y = (cart2.y - cart1.y) / (measurements[2].timestamp - measurements[1].timestamp);

		// We now have an initial position and velocity
		cartesian_point curr_pos = cart1;
		cartesian_point curr_vel = vel1;

		// For the analytical method and the numerical method
		dt = (measurements[2].timestamp - measurements[1].timestamp);

		// Gauss' method uses 3 input points
		// System of linear equations uses 5 input points
		const size_t num_points_needed = 3;

		orbit_points.push_back(curr_pos);
		orbit_velocities.push_back(curr_vel);

		// Calculate (num_points_needed - 1) position and velocity data,
		// using 4th order symplectic integration
		for (size_t i = 1; i < num_points_needed; i++)
		{
			vector_3 curr_pos_3d(curr_pos.x, curr_pos.y, 0);
			vector_3 curr_vel_3d(curr_vel.x, curr_vel.y, 0);
			proceed_symplectic4(curr_pos_3d, curr_vel_3d, grav_constant, dt);
			curr_pos = cartesian_point(curr_pos_3d.x, curr_pos_3d.y);
			curr_vel = cartesian_point(curr_vel_3d.x, curr_vel_3d.y);

			orbit_points.push_back(curr_pos);
			orbit_velocities.push_back(curr_vel);
		}

		// This needs 3 points
		vector_3d r1_ = { orbit_points[0].x, orbit_points[0].y, 0.0 };
		vector_3d r2_ = { orbit_points[1].x, orbit_points[1].y, 0.0 };
		vector_3d r3_ = { orbit_points[2].x, orbit_points[2].y, 0.0 };

		double t1_ = 0.0;
		double t2_ = dt;
		double t3_ = 2 * dt;

		//gauss_method(r1_, r2_, r3_, t1_, t2_, t3_);




		EllipseParameters ep = fit_ellipse3x3(orbit_points);





		// Bootstrap the numerical integration,
		// to visualy double check the results
		orbiter_pos.x = orbit_points[0].x;
		orbiter_pos.y = orbit_points[0].y;
		orbiter_vel.x = orbit_velocities[0].x;
		orbiter_vel.y = orbit_velocities[0].y;
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
	glLineWidth(3.0f);



	glBegin(GL_POINTS);

	glColor3f(1.0, 1.0, 1.0);
	glVertex3d(sun_pos.x, sun_pos.y, sun_pos.z);

	//glColor3f(1.0, 0.0, 1.0);

	//for (size_t i = 0; i < orbit_points.size(); i++)
	//	glVertex3d(orbit_points[i].x, orbit_points[i].y, 0);

	glColor3f(0.0, 1.0, 0.0);

	for (size_t i = 0; i < positions.size(); i++)
		glVertex3d(positions[i].x, positions[i].y, 0);

	glEnd();


	glPushMatrix();

	glColor3f(1.0, 0.5, 0.0);

	glTranslated(global_ep.centerX, global_ep.centerY, 0);
	glRotated((global_ep.angle) * (2 * pi) * 360.0, 0, 0, 1);
	glTranslated(-global_ep.centerX, -global_ep.centerY, 0);

	DrawEllipse(global_ep.centerX, global_ep.centerY, global_ep.semiMinor, global_ep.semiMajor, 100);

	glPopMatrix();




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


