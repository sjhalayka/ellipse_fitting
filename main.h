#ifndef main_H
#define main_H

#include "uv_camera.h"
#include "custom_math.h"
using custom_math::vector_3;



//using namespace custom_math;

#include <cstdlib>
#include <GL/glut.h>       //GLUT Library
#include <stdexcept>
#include <optional>
using namespace std;

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;
using std::istringstream;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <set>
using std::set;

#include <map>
using std::map;

#include <utility>
using std::pair;

#include <random>


void idle_func(void);
void init_opengl(const int &width, const int &height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);

void render_string(int x, const int y, void *font, const string &text);
void draw_objects(void);

const double speed_of_light = 299792458;
const double grav_constant = 6.673e-11;
const double sun_mass = 1.989e30;
const double pi = 4.0 * atan(1.0);

custom_math::vector_3 sun_pos(0, 0, 0);

custom_math::vector_3 orbiter_pos(0, 0, 0);
custom_math::vector_3 orbiter_vel(0, 0, 0);

vector<custom_math::vector_3> positions;

custom_math::vector_3 background_colour(0.0f, 0.0f, 0.0f);
custom_math::vector_3 control_list_colour(1.0f, 1.0f, 1.0f);

bool draw_axis = true;
bool draw_control_list = true;

uv_camera main_camera;

GLint win_id = 0;
GLint win_x = 800, win_y = 600;
float camera_w = 5e11f;

float camera_fov = 45;
float camera_x_transform = 0;
float camera_y_transform = 0;
float u_spacer = 0.01f;
float v_spacer = 0.5f*u_spacer;
float w_spacer = 0.1f;
float camera_near = 0.01f;
float camera_far = 1e15f;

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;


#endif
