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
























// Structure to represent a 3D vector
class vector_3d {

public:
	double x, y, z;

	//vector_3d operator-(const vector_3d& rhs) const
	//{
	//	vector_3d ret;
	//	ret.x = x - rhs.x;
	//	ret.y = y - rhs.y;
	//	ret.z - z - rhs.z;
	//	return ret;
	//}

	//vector_3d operator+(const vector_3d& rhs) const
	//{
	//	vector_3d ret;
	//	ret.x = x + rhs.x;
	//	ret.y = y + rhs.y;
	//	ret.z - z + rhs.z;
	//	return ret;
	//}
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
	double a_orbit = 1.0 / (2.0 / magnitude(r2) - dotProduct(v2, v2) / (sun_mass * grav_constant));

	// Calculate the eccentricity
	double e = sqrt(1.0 - (magnitude(h) * magnitude(h) / ((sun_mass * grav_constant) * a_orbit)));

	// Calculate the semi-minor axis
	double b_orbit = a_orbit * sqrt(1.0 - e * e);

	// Calculate the center of the ellipse
	vector_3d center;
	center.x = 0;
	center.z = 0;

	const double dx = (magnitude(r2) - magnitude(r1)) + (magnitude(r3) - magnitude(r2));

	if (dx < 0)
		center.y = -a_orbit * e;
	else
		center.y = a_orbit * e;

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
vector<cartesian_point> orbit_points(6);
vector<cartesian_point> orbit_velocities(6);



double calculateOrbitRadius(double omega, double omegaDot, double GM) 
{
	// Check for valid angular velocity
	if (omega <= 0) 
	{
		throw std::invalid_argument("Angular velocity must be positive.");
	}

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
	const double dt = 10000; // 10000 seconds == 2.77777 hours

	static bool calculated_ellipse = false;

	if (true == calculated_ellipse)
	{
		proceed_symplectic4(mercury_pos, mercury_vel, grav_constant, dt);
		positions.push_back(mercury_pos);
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

		// Produce 2 radii
		vector<double> radii_data(2);

		// Constant angular velocity, for example
		//double omega = 4.31e-8; // Ceres average angular velocity
		//double omega_min = 1.99e-7; // Earth average angular velocity


		double prev_omega = (measurements[1].azimuth - measurements[0].azimuth) / (measurements[1].timestamp - measurements[0].timestamp);

		for (size_t i = 0; i < measurements.size() - 1; i++)
		{
			// Variable angular velocity
			const double omega = (measurements[i + 1].azimuth - measurements[i].azimuth) / (measurements[i + 1].timestamp - measurements[i].timestamp);

			const double d_omega = 0.0;// (omega - prev_omega) / (measurements[i + 1].timestamp - measurements[i].timestamp);
			prev_omega = omega;
	
			const double radius = calculateOrbitRadius(omega, d_omega, sun_mass*grav_constant);

			radii_data[i] = radius;
		}

		// Produce input data
		const double angle1 = measurements[1].azimuth;
		const double r1 = radii_data[0];

		const double angle2 = measurements[2].azimuth;
		const double r2 = radii_data[1];

		// Convert input data to Cartesian coordinates
		cartesian_point cart1 = to_cartesian(r1, angle1);
		cartesian_point cart2 = to_cartesian(r2, angle2);

		cartesian_point vel1;
		vel1.x = (cart2.x - cart1.x) / (measurements[2].timestamp - measurements[1].timestamp);
		vel1.y = (cart2.y - cart1.y) / (measurements[2].timestamp - measurements[1].timestamp);

		cartesian_point curr_pos = cart1;
		cartesian_point curr_vel = vel1;

		// For the analytical method
		double dt_ = (measurements[2].timestamp - measurements[1].timestamp);

		orbit_points[0] = curr_pos;
		orbit_velocities[0] = curr_vel;

		const size_t num_points_needed = 6;

		for (size_t i = 1; i < num_points_needed; i++)
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

		vector_3d r1_ = { orbit_points[0].x, orbit_points[0].y, 0.0 };
		vector_3d r2_ = { orbit_points[1].x, orbit_points[1].y, 0.0 };
		vector_3d r3_ = { orbit_points[2].x, orbit_points[2].y, 0.0 };

		double t1_ = 0.0;
		double t2_ = dt_;
		double t3_ = 2 * dt_;

		gaussMethod(r1_, r2_, r3_, t1_, t2_, t3_);






		// Bootstrap the numerical integration,
		// to double check the results
		mercury_pos.x = orbit_points[0].x;
		mercury_pos.y = orbit_points[0].y;
		mercury_vel.x = orbit_velocities[0].x;
		mercury_vel.y = orbit_velocities[0].y;
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

	//if (carts.size() > 0)
	//{
	//	glColor3f(1.0, 0.0, 0.0);
	//	glVertex3d(carts[0].x, carts[0].y, 0);

	//	glColor3f(0.0, 1.0, 0.0);
	//	glVertex3d(carts[1].x, carts[1].y, 0);
	//}

	//glColor3f(1.0, 0.0, 1.0);


	for (size_t i = 0; i < orbit_points.size(); i++)
		glVertex3d(orbit_points[i].x, orbit_points[i].y, 0);


	glColor3f(0.0, 1.0, 0.0);

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


