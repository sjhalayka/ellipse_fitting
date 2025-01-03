#include "main.h"


#include <Eigen/Dense>
using namespace Eigen;

// Function to fit an ellipse using Gauss method
VectorXd fitEllipse(const std::vector<double>& x, const std::vector<double>& y) {
	size_t n = x.size();
	if (n < 5) {
		std::cerr << "Need at least 5 points to fit an ellipse." << std::endl;
		return VectorXd(0);
	}

	// Step 1: Construct the design matrix D
	MatrixXd D(n, 6);
	for (size_t i = 0; i < n; ++i) {
		D(i, 0) = x[i] * x[i];
		D(i, 1) = x[i] * y[i];
		D(i, 2) = y[i] * y[i];
		D(i, 3) = x[i];
		D(i, 4) = y[i];
		D(i, 5) = 1.0;
	}

	// Step 2: Compute the scatter matrix S = D^T * D
	MatrixXd S = D.transpose() * D;

	// Step 3: Define the constraint matrix C
	MatrixXd C(6, 6);
	C << 0, 0, 2, 0, 0, 0,
		0, -1, 0, 0, 0, 0,
		2, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0;

	// Step 4: Solve the generalized eigenvalue problem S * v = ? * C * v
	EigenSolver<MatrixXd> solver(S);
	MatrixXd eigenvalues = solver.eigenvalues().real();
	MatrixXd eigenvectors = solver.eigenvectors().real();

	// Find the solution satisfying the ellipse constraint
	int bestIndex = -1;
	for (int i = 0; i < eigenvalues.size(); ++i) {
		VectorXd v = eigenvectors.col(i);
		if (v.transpose() * C * v > 0) {
			bestIndex = i;
			break;
		}
	}

	if (bestIndex == -1) {
		std::cerr << "No valid ellipse found." << std::endl;
		return VectorXd(0);
	}

	return eigenvectors.col(bestIndex).transpose();
}



class GaussEllipseFit {
public:
	static VectorXd fitEllipse(const std::vector<double>& x, const std::vector<double>& y) {
		if (x.size() != y.size() || x.size() < 5) {
			throw std::runtime_error("Need at least 5 points to fit an ellipse");
		}

		// Build design matrix D
		Eigen::MatrixXd D(x.size(), 6);
		for (size_t i = 0; i < x.size(); i++) {
			D(i, 0) = x[i] * x[i];
			D(i, 1) = x[i] * y[i];
			D(i, 2) = y[i] * y[i];
			D(i, 3) = x[i];
			D(i, 4) = y[i];
			D(i, 5) = 1.0;
		}

		// Build scatter matrix S
		Eigen::MatrixXd S = D.transpose() * D;

		// Build constraint matrix C
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6);
		C(0, 2) = 2;
		C(2, 0) = 2;
		C(1, 1) = -1;

		// Solve generalized eigensystem
		Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(S, C);

		// Find the positive eigenvalue
		Eigen::VectorXd eigenvalues = ges.eigenvalues().real();
		int positiveIndex = -1;
		for (int i = 0; i < eigenvalues.size(); i++) {
			if (eigenvalues(i) > 0) {
				positiveIndex = i;
				break;
			}
		}

		if (positiveIndex == -1) {
			throw std::runtime_error("No valid ellipse found");
		}

		// Get the coefficients
		Eigen::VectorXd coefficients = ges.eigenvectors().real().col(positiveIndex);

		// Extract ellipse parameters
		return coefficients;// extractEllipseParameters(coefficients);
	}

	EllipseParameters extractEllipseParameters(const Eigen::VectorXd& coefficients) 
	{
		double a = coefficients(0);
		double b = coefficients(1) / 2;
		double c = coefficients(2);
		double d = coefficients(3);
		double e = coefficients(4);
		double f = coefficients(5);

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
		params.rotation = theta;

		return params;
	}
};


struct EllipseFoci {
	double focus1X;
	double focus1Y;
	double focus2X;
	double focus2Y;
};




EllipseFoci calculateFoci(const EllipseParameters& params) {
	// c� = a� - b�
	double c = sqrt(params.semiMajor * params.semiMajor -
		params.semiMinor * params.semiMinor);

	// Calculate foci before rotation
	double fx = c;
	double fy = 0;

	// Rotate foci
	double fx1 = fx * cos(params.rotation) - fy * sin(params.rotation);
	double fy1 = fx * sin(params.rotation) + fy * cos(params.rotation);

	// Return both foci (translated to ellipse center)
	EllipseFoci foci;
	foci.focus1X = params.centerX + fx1;
	foci.focus1Y = params.centerY + fy1;
	foci.focus2X = params.centerX - fx1;
	foci.focus2Y = params.centerY - fy1;

	return foci;
}


void DrawEllipse(float cx, float cy, float rx, float ry, int num_segments)
{
	float theta = 2 * 3.1415926f / float(num_segments);
	float c = cosf(theta);//precalculate the sine and cosine
	float s = sinf(theta);
	float t;

	float x = 1;//we start at angle = 0 
	float y = 0;

	glBegin(GL_LINE_LOOP);
	for (int ii = 0; ii < num_segments; ii++)
	{
		//apply radius and offset
		glVertex2f(x * rx + cx, y * ry + cy);//output vertex 

		//apply the rotation matrix
		t = x;
		x = c * x - s * y;
		y = s * t + c * y;
	}
	glEnd();
}




int main(int argc, char **argv)
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

custom_math::vector_3 grav_acceleration(const custom_math::vector_3 &pos, const custom_math::vector_3 &vel, const double G)
{
	custom_math::vector_3 grav_dir = sun_pos - pos;

	double distance = grav_dir.length();

	grav_dir.normalize();
	custom_math::vector_3 accel = grav_dir * (G*sun_mass / pow(distance, 2.0));

	return accel;
}



void proceed_Euler(custom_math::vector_3 &pos, custom_math::vector_3 &vel, const double G, const double dt)
{
	custom_math::vector_3 accel = grav_acceleration(pos, vel, G);

	vel += accel * dt;
	pos += vel * dt;
}

void proceed_RK4(custom_math::vector_3 &pos, custom_math::vector_3 &vel, const double G, const double dt)
{
	static const double one_sixth = 1.0 / 6.0;

	custom_math::vector_3 k1_velocity = vel;
	custom_math::vector_3 k1_acceleration = grav_acceleration(pos, k1_velocity, G);
	custom_math::vector_3 k2_velocity = vel + k1_acceleration * dt*0.5;
	custom_math::vector_3 k2_acceleration = grav_acceleration(pos + k1_velocity * dt*0.5, k2_velocity, G);
	custom_math::vector_3 k3_velocity = vel + k2_acceleration * dt*0.5;
	custom_math::vector_3 k3_acceleration = grav_acceleration(pos + k2_velocity * dt*0.5, k3_velocity, G);
	custom_math::vector_3 k4_velocity = vel + k3_acceleration * dt;
	custom_math::vector_3 k4_acceleration = grav_acceleration(pos + k3_velocity * dt, k4_velocity, G);

	vel += (k1_acceleration + (k2_acceleration + k3_acceleration)*2.0 + k4_acceleration)*one_sixth*dt;
	pos += (k1_velocity + (k2_velocity + k3_velocity)*2.0 + k4_velocity)*one_sixth*dt;
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




void idle_func(void)
{
	static size_t frame_count = 0;

	frame_count++;

	const double dt = 10000; // 10000 seconds == 2.77777 hours

	// Pick an integrator:

	//proceed_Euler(mercury_pos, mercury_vel, grav_constant, dt);
	//proceed_RK4(mercury_pos, mercury_vel, grav_constant, dt);
	proceed_symplectic4(mercury_pos, mercury_vel, grav_constant, dt);

    positions.push_back(mercury_pos);


	static bool calculated_ellipse = false;


	if (calculated_ellipse == false && frame_count % 123 == 0)
		ellipse_positions.push_back(positions[positions.size() - 1]);




	if (false == calculated_ellipse && ellipse_positions.size() == 20)
	{
		calculated_ellipse = true;



		//vector<Point> pointsf;

		//for (size_t i = 0; i < ellipse_positions.size(); i++)
		//{
		//	pointsf.push_back(Point(ellipse_positions[i].x, ellipse_positions[i].y));
		//}

		//RotatedRect box = cv::fitEllipse(pointsf);


		//Mat cimage = Mat::zeros(800, 600, CV_8UC3);

		//ellipse(cimage, box, Scalar(0, 0, 255), 1, LINE_AA);

		////ellipse(cimage, box, Scalar(0, 0, 255), 1, LINE_AA);
		//
		//imshow("result", cimage);

		double largest_distance = 0;
		double avg_distance = 0;
		vector<double> distances(ellipse_positions.size(), 0);

		for (size_t i = 0; i < ellipse_positions.size(); i++)
		{
			if (ellipse_positions[i].length() > largest_distance)
				largest_distance = ellipse_positions[i].length();

			avg_distance += ellipse_positions[i].length();
			distances[i] = ellipse_positions[i].length();
		}

		avg_distance /= ellipse_positions.size();

		global_ep.centerX = 0;
		global_ep.centerY = 0;
		global_ep.semiMajor = largest_distance;
		global_ep.semiMinor = largest_distance;
		global_ep.rotation = 0;




		//for (size_t i = 0; i < global_errors.size(); i++)
		//{
		//	const double val = (ellipse_positions[i].x * ellipse_positions[i].x / (global_ep.semiMinor * global_ep.semiMinor)) + (ellipse_positions[i].y * ellipse_positions[i].y / (global_ep.semiMajor * global_ep.semiMajor)) - 1.0;

		//	global_errors[i] = abs(val);

		//	global_total_error += global_errors[i];
		//}
		



		for (size_t k = 0; k < 1000; k++)
		{
			vector<double> global_errors(ellipse_positions.size(), 0);
			double global_total_error = 0;

			global_total_error = 0;

			for (size_t i = 0; i < global_errors.size(); i++)
			{
				const double val = (ellipse_positions[i].x * ellipse_positions[i].x / (global_ep.semiMinor * global_ep.semiMinor)) + (ellipse_positions[i].y * ellipse_positions[i].y / (global_ep.semiMajor * global_ep.semiMajor)) - 1.0;

				global_errors[i] = abs(val);

				global_total_error += global_errors[i];
			}



			const size_t refinement_count = 1;
















			if(k % 2 == 0)
			for (size_t i = 0; i < refinement_count; i++)
			{
				EllipseParameters local_ep;
				local_ep.centerX = global_ep.centerX;
				local_ep.centerY = global_ep.centerY;
				local_ep.semiMajor = global_ep.semiMajor;
				local_ep.semiMinor = global_ep.semiMinor * 0.9;
				local_ep.rotation = 0;

				vector<double> local_errors(global_errors.size(), 0);
				double local_total_error = 0;

				EllipseFoci ef = calculateFoci(local_ep);

				for (size_t j = 0; j < local_errors.size(); j++)
				{
					const double val = (ellipse_positions[j].x * ellipse_positions[j].x / (local_ep.semiMinor * local_ep.semiMinor)) + (ellipse_positions[j].y * ellipse_positions[j].y / (local_ep.semiMajor * local_ep.semiMajor)) - 1.0;
					local_errors[j] = abs(val);

					local_total_error += local_errors[j];
				}

				static double last_local_total_error = DBL_MAX;

				if (local_total_error < global_total_error)
				{
					if (last_local_total_error < local_total_error)
						continue;

					last_local_total_error = local_total_error;

					EllipseFoci ef = calculateFoci(local_ep);

					//cout << local_total_error << " " << global_total_error << endl;

					global_ep.centerX = local_ep.centerX;
					global_ep.centerY = ef.focus1X;// local_ep.centerY;
					global_ep.semiMajor = local_ep.semiMajor;
					global_ep.semiMinor = local_ep.semiMinor;
					global_ep.rotation = 0;

					//continue;
					//break;
				}

				//global_total_error = local_total_error;
			}
			else
			for (size_t i = 0; i < refinement_count; i++)
			{
				EllipseParameters local_ep;
				local_ep.centerX = global_ep.centerX;
				local_ep.centerY = global_ep.centerY;
				local_ep.semiMajor = global_ep.semiMajor * 0.9;
				local_ep.semiMinor = global_ep.semiMinor;
				local_ep.rotation = 0;

				vector<double> local_errors(global_errors.size(), 0);
				double local_total_error = 0;
				
				if (local_ep.semiMajor <= local_ep.semiMinor)
				{
					//EllipseFoci ef = calculateFoci(local_ep);

					//global_ep.centerX = local_ep.centerX;
					//global_ep.centerY = local_ep.centerY;
					//global_ep.semiMajor = local_ep.semiMajor;
					//global_ep.semiMinor = local_ep.semiMinor;
					//global_ep.rotation = 0;

					//break;
				}

				for (size_t j = 0; j < local_errors.size(); j++)
				{
					const double val = (ellipse_positions[j].x * ellipse_positions[j].x / (local_ep.semiMinor * local_ep.semiMinor)) + (ellipse_positions[j].y * ellipse_positions[j].y / (local_ep.semiMajor * local_ep.semiMajor)) - 1.0;
					local_errors[j] = abs(val);

					local_total_error += local_errors[j];
				}

				static double last_local_total_error = DBL_MAX;

				if (local_total_error < global_total_error)
				{
					if (last_local_total_error < local_total_error)
						continue;

					last_local_total_error = local_total_error;

					EllipseFoci ef = calculateFoci(local_ep);

					cout << local_total_error << " " << global_total_error << endl;
					cout << "lala" << endl;

					global_ep.centerX = local_ep.centerX;
					global_ep.centerY = ef.focus1X;// local_ep.centerY;
					global_ep.semiMajor = local_ep.semiMajor;
					global_ep.semiMinor = local_ep.semiMinor;
					global_ep.rotation = 0;

					//global_total_error = local_total_error;
					//continue;// break;

				}
			}



			//for (size_t i = 0; i < refinement_count; i++)
			//{
			//	EllipseParameters local_ep;
			//	local_ep.centerX = global_ep.centerX;
			//	local_ep.centerY = global_ep.centerY;
			//	local_ep.semiMajor = global_ep.semiMajor * 0.99;
			//	local_ep.semiMinor = global_ep.semiMinor;
			//	local_ep.rotation = 0;

			//	//if (local_ep.semiMajor <= local_ep.semiMinor)
			//	//{
			//	//	k = 1000;
			//	//	//EllipseFoci ef = calculateFoci(local_ep);

			//	//	//global_ep.centerX = local_ep.centerX;
			//	//	//global_ep.centerY = local_ep.centerY;
			//	//	//global_ep.semiMajor = local_ep.semiMajor;
			//	//	//global_ep.semiMinor = local_ep.semiMinor;
			//	//	//global_ep.rotation = 0;

			//	//	break;
			//	//}

			//	vector<double> local_errors(global_errors.size(), 0);
			//	double local_total_error = 0;



			//	for (size_t j = 0; j < local_errors.size(); j++)
			//	{
			//		const double val = (ellipse_positions[j].x * ellipse_positions[j].x / (local_ep.semiMinor * local_ep.semiMinor)) + (ellipse_positions[j].y * ellipse_positions[j].y / (local_ep.semiMajor * local_ep.semiMajor)) - 1.0;
			//		local_errors[j] = abs(val);

			//		local_total_error += local_errors[j];
			//	}




			//	static double last_local_total_error = DBL_MAX;

			//	if (local_total_error < global_total_error)
			//	{
			//		if (last_local_total_error < local_total_error)
			//			break;

			//		last_local_total_error = local_total_error;

			//		EllipseFoci ef = calculateFoci(local_ep);

			//		cout << local_total_error << " " << global_total_error << endl;

			//		global_ep.centerX = local_ep.centerX;
			//		global_ep.centerY = ef.focus1X;// local_ep.centerY;
			//		global_ep.semiMajor = local_ep.semiMajor;
			//		global_ep.semiMinor = local_ep.semiMinor;
			//		global_ep.rotation = 0;

			//		//global_total_error = local_total_error;

			//		break;
			//	}
			//}






		}




	}


    glutPostRedisplay();
}

void init_opengl(const int &width, const int &height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH);
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

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutSetWindow(win_id);
	glutReshapeWindow(win_x, win_y);
	glViewport(0, 0, win_x, win_y);

	main_camera.Set(main_camera.u, main_camera.v, main_camera.w, main_camera.fov, win_x, win_y, camera_near, camera_far);
}

// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
void render_string(int x, const int y, void *font, const string &text)
{
	for(size_t i = 0; i < text.length(); i++)
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
  

    glPointSize(2.0);
    glLineWidth(2.0f);

    
    glBegin(GL_POINTS);
    glVertex3d(sun_pos.x, sun_pos.y, sun_pos.z);
    
    glColor3f(1.0, 1.0, 1.0);
    
    for(size_t i = 0; i < ellipse_positions.size(); i++)
		glVertex3d(ellipse_positions[i].x, ellipse_positions[i].y, ellipse_positions[i].z);
    
    glEnd();
    

	glColor3f(1.0, 0.5, 0.0);
	//global_ep.centerX = local_ep.centerX;
	//global_ep.centerY = local_ep.centerY;
	//global_ep.semiMajor = local_ep.semiMajor;
	//global_ep.semiMinor = local_ep.semiMinor;
	//global_ep.rotation = 0;


	EllipseFoci ef = calculateFoci(global_ep);

	//cout << ef.focus1X << " " << ef.focus1Y << endl;
	//cout << ef.focus2X << " " << ef.focus2Y << endl;
	//cout << endl;

	//glTranslated(0, ef.focus1Y, 0);

	DrawEllipse(global_ep.centerX, global_ep.centerY, global_ep.semiMinor, global_ep.semiMajor, 100);


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
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	// Draw the model's components using OpenGL/GLUT primitives.
	draw_objects();

	if(true == draw_control_list)
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
		render_string(10, static_cast<int>(start + 1*break_size), GLUT_BITMAP_HELVETICA_18, string("  LMB + drag: Rotate camera"));
		render_string(10, static_cast<int>(start + 2*break_size), GLUT_BITMAP_HELVETICA_18, string("  RMB + drag: Zoom camera"));

		render_string(10, static_cast<int>(start + 4*break_size), GLUT_BITMAP_HELVETICA_18, string("Keyboard controls:"));
        render_string(10, static_cast<int>(start + 5*break_size), GLUT_BITMAP_HELVETICA_18, string("  w: Draw axis"));
		render_string(10, static_cast<int>(start + 6*break_size), GLUT_BITMAP_HELVETICA_18, string("  e: Draw text"));
		render_string(10, static_cast<int>(start + 7*break_size), GLUT_BITMAP_HELVETICA_18, string("  u: Rotate camera +u"));
		render_string(10, static_cast<int>(start + 8*break_size), GLUT_BITMAP_HELVETICA_18, string("  i: Rotate camera -u"));
		render_string(10, static_cast<int>(start + 9*break_size), GLUT_BITMAP_HELVETICA_18, string("  o: Rotate camera +v"));
		render_string(10, static_cast<int>(start + 10*break_size), GLUT_BITMAP_HELVETICA_18, string("  p: Rotate camera -v"));


		
        custom_math::vector_3 eye = main_camera.eye;
		custom_math::vector_3 eye_norm = eye;
		eye_norm.normalize();

		oss.clear();
		oss.str("");		
		oss << "Camera position: " << eye.x << ' ' << eye.y << ' ' << eye.z;
		render_string(10, static_cast<int>(win_y - 2*break_size), GLUT_BITMAP_HELVETICA_18, oss.str());

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
	switch(tolower(key))
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
	if(GLUT_LEFT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			lmb_down = true;
		else
			lmb_down = false;
	}
	else if(GLUT_MIDDLE_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			mmb_down = true;
		else
			mmb_down = false;
	}
	else if(GLUT_RIGHT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
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

	if(true == lmb_down && (0 != mouse_delta_x || 0 != mouse_delta_y))
	{
		main_camera.u -= static_cast<float>(mouse_delta_y)*u_spacer;
		main_camera.v += static_cast<float>(mouse_delta_x)*v_spacer;
	}
	else if(true == rmb_down && (0 != mouse_delta_y))
	{
		main_camera.w -= static_cast<float>(mouse_delta_y)*w_spacer;

		if(main_camera.w < 1.1f)
			main_camera.w = 1.1f;

	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}




