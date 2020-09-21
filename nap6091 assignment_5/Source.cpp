// Nicholas Pyfrom-Day		nap6091
// CMPS 415
// Assignment 5 source code
/* image "stars.ppm" obtained from http://photojournal.jpl.nasa.gov/jpeg/PIA15256.jpg */

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <gmtl/gmtl.h> 
#include <stdlib.h>
#include <stdio.h>

#include "text.h"

#pragma comment (lib, "opengl32.lib")
#pragma comment (lib, "glew32.lib")
#pragma comment (lib, "glfw3.lib")

using namespace std;

const int num_birds = 6;
const float pi = M_PI;
const float r_sph = 0.4;								// radius of sphere
const float r_cyl = 0.005;								// radius of cylinder
const float flight_radius = r_sph + 0.15;
const float threshold = (11 * pi) / 30;
const float time = 0.05;
const int nm = 10;										// number of meridians is 10
const int np = 10;										// number of parallels is 10
const int sphere_size = (np + 2)*(2 * (nm + 1) + 1) + 2 * (nm + 1);
const int cyl_size = 2 * ((np - 1) * (nm + 1)) + np - 2;	// size of index list for sphere and cylinder
const int num_cyl_vert = np*(nm + 1);						// number of vertices
const int num_sphere_vert = (np + 4)*(nm + 1);
const float phi_inc = (2 * pi) / nm;					// initialize the increment value for phi
const float theta_inc = pi / (np - 1);					// initialize the increment value for theta
float omega[num_birds];
float distance[num_birds][num_birds];
float wing1_angle = pi / 6;								// current angle of joint 1
float wing1_angle_inc = pi / 6;							// angle to rotate wing 1
float wing2_angle = -pi / 6;							// current angle of joint 2
float wing2_angle_inc = -pi / 6;						// angle to rotate wing 2
double cam_previous_x_pos = 0.0;
double cam_previous_y_pos = 0.0;
double lon = 0.0;
double lat = 0.0;
double azi = 0.0;
double ele = 0.0;
bool cam_1_active = true;
bool draw_axes = false;
bool light1 = true;
bool light2 = true;
bool update = false;

/* -------------------------Model Matrices-------------------------------------- */
gmtl::Matrix44f M_sph_to_world;		// sphere model matrix
gmtl::Matrix44f M_cyl_to_world;		// cylinder model matrix
gmtl::Matrix44f M_bird[num_birds];				// bird model matrix
gmtl::Matrix44f M_wing1[num_birds];			// wing 1 model matrix
gmtl::Matrix44f M_wing2[num_birds];			// wing 2 model matrix
											/* --------------------------Transforms----------------------------------------- */
gmtl::Matrix44f T_cyl;					// cyl object wrt W, R, or O frame
gmtl::Matrix44f T_bird_obj_wrt_B;		// bird object wrt bird base frame
gmtl::Matrix44f T_wing1_obj_wrt_J1;		// wing object 1 wrt joint 1 frame
gmtl::Matrix44f T_wing2_obj_wrt_J2;		// wing object 2 wrt joint 2 frame
gmtl::Matrix44f T_J1_wrt_B;				// {J1} wrt {bird base}
gmtl::Matrix44f T_J2_wrt_B;				// {J2} wrt {bird base}
gmtl::Matrix44f T_B_wrt_O;				// {bird base} wrt {O}
gmtl::Matrix44f T_O_wrt_R;				// {O} wrt {R}
gmtl::Matrix44f Q[num_birds];						// {R} wrt {W}
gmtl::Matrix44f cam_1_wrt_W;			// camera pose describing {camera 1} wrt {W}
gmtl::Matrix44f Cam_2;					// camera pose describing {camera 2} wrt {O}
gmtl::Matrix44f Cam_2_wrt_W;			// camera pose describing {camera 2} wrt {W}
										/* --------------------------View Matrices-------------------------------------- */
gmtl::Matrix44f V_cam_1;			// view matrix for camera 1
gmtl::Matrix44f V_cam_2;			// view matrix for camera 2
gmtl::Matrix44f V_current;			// the active view matrix
gmtl::Matrix44f P;
/* --------------------------Rotation Matrices---------------------------------- */
gmtl::Matrix44f x_rot;		// transform matrix for positive X rotation
gmtl::Matrix44f z_rot;		// transform matrix for positive Z rotation
gmtl::Matrix44f y_rot;
gmtl::Matrix44f z_translation;

gmtl::Vec4f light_position;
gmtl::Vec4f world_origin;
gmtl::Vec3f light_color_1;
gmtl::Vec3f light_color_2;

char* filetobuf(char *file)
{
	FILE *fptr;
	FILE **ptr = &fptr;
	long length;
	char *buf;

	fopen_s(ptr, file, "rb");
	if (!fptr)
		return NULL;
	fseek(fptr, 0, SEEK_END);
	length = ftell(fptr);
	buf = (char*)malloc(length + 1);
	fseek(fptr, 0, SEEK_SET);
	fread(buf, length, 1, fptr);
	fclose(fptr);
	buf[length] = 0;

	return buf;
}

void error_callback(int error, const char* description)
{
	fprintf(stderr, "Error: %s\n", description);
}

GLuint setupShaderProgram()
{
	GLuint vertex_shader, fragment_shader, shader_program;
	int IsCompiled_VS, IsCompiled_FS, IsLinked, max_length;
	char *vertex_shader_log;
	char *fragment_shader_log;
	char *shader_program_log;

	char* vertex_source = filetobuf("OpenGL_Example.vert");
	char* fragment_source = filetobuf("OpenGL_Example.frag");

	vertex_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex_shader, 1, &vertex_source, NULL);
	glCompileShader(vertex_shader);

	glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &IsCompiled_VS);
	if (IsCompiled_VS == GL_FALSE)
	{
		glGetShaderiv(vertex_shader, GL_INFO_LOG_LENGTH, &max_length);

		vertex_shader_log = (char *)malloc(max_length);

		glGetShaderInfoLog(vertex_shader, max_length, &max_length, vertex_shader_log);
		printf("Error: %s", vertex_shader_log);
		free(vertex_shader_log);
		free(vertex_source);
		return 0;
	}

	fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment_shader, 1, &fragment_source, NULL);
	glCompileShader(fragment_shader);
	glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &IsCompiled_FS);
	if (IsCompiled_FS == GL_FALSE)
	{
		glGetShaderiv(fragment_shader, GL_INFO_LOG_LENGTH, &max_length);

		fragment_shader_log = (char *)malloc(max_length);

		glGetShaderInfoLog(fragment_shader, max_length, &max_length, fragment_shader_log);
		printf("Error: %s", fragment_shader_log);
		free(fragment_shader_log);
		free(vertex_source);
		free(fragment_source);
		return 0;
	}

	shader_program = glCreateProgram();

	glAttachShader(shader_program, vertex_shader);
	glAttachShader(shader_program, fragment_shader);

	glLinkProgram(shader_program);

	glGetProgramiv(shader_program, GL_LINK_STATUS, (int *)&IsLinked);
	if (IsLinked == GL_FALSE)
	{
		glGetProgramiv(shader_program, GL_INFO_LOG_LENGTH, &max_length);

		shader_program_log = (char *)malloc(max_length);

		glGetProgramInfoLog(shader_program, max_length, &max_length, shader_program_log);
		printf("Error: %s", shader_program_log);
		free(shader_program_log);
		free(vertex_source);
		free(fragment_source);
		return 0;
	}

	free(vertex_source);
	free(fragment_source);

	return shader_program;

}

void init(GLuint shader_program, GLuint &pos_loc_out, GLuint &color_loc_out, GLuint &normals_loc_out, GLuint &uv_out,
	GLuint &m_loc_out, GLuint &mvp_loc_out, GLuint &normMat_loc_out,
	GLuint &light_loc_out, GLuint &light_color_1_loc_out, GLuint &light_color_2_loc_out, GLuint &world_origin_loc_out, GLuint &texture_loc_out)
{
	pos_loc_out = glGetAttribLocation(shader_program, "in_position");
	color_loc_out = glGetAttribLocation(shader_program, "in_color");
	normals_loc_out = glGetAttribLocation(shader_program, "in_normals");
	uv_out = glGetAttribLocation(shader_program, "in_uv");
	m_loc_out = glGetUniformLocation(shader_program, "M");
	mvp_loc_out = glGetUniformLocation(shader_program, "MVP");
	normMat_loc_out = glGetUniformLocation(shader_program, "normalMatrix");
	light_loc_out = glGetUniformLocation(shader_program, "light_position");
	light_color_1_loc_out = glGetUniformLocation(shader_program, "light_color_1");
	light_color_2_loc_out = glGetUniformLocation(shader_program, "light_color_2");
	world_origin_loc_out = glGetUniformLocation(shader_program, "world_origin");
	texture_loc_out = glGetUniformLocation(shader_program, "texture_Colors");
}

void set_x_rot(double angle) {
	gmtl::identity(x_rot);
	x_rot[1][1] = x_rot[2][2] = cos(angle);
	x_rot[1][2] = -sin(angle);
	x_rot[2][1] = sin(angle);
	x_rot.setState(gmtl::Matrix44f::AFFINE);
}

void set_y_rot(double angle) {
	gmtl::identity(y_rot);
	y_rot[2][2] = y_rot[0][0] = cos(angle);
	y_rot[0][2] = sin(angle);
	y_rot[2][0] = -sin(angle);
	y_rot.setState(gmtl::Matrix44f::AFFINE);
}

void set_z_rot(double angle) {
	gmtl::identity(z_rot);
	z_rot[1][1] = z_rot[0][0] = cos(angle);
	z_rot[0][1] = -sin(angle);
	z_rot[1][0] = sin(angle);
	z_rot.setState(gmtl::Matrix44f::AFFINE);
}

void set_z_translation(double delta_z) {
	z_translation[2][3] = z_translation[2][3] + delta_z;
}

void set_T_cyl(int frame, int axis) {
	if (frame == 0) {						// frame: 0 for cyl_wrt_W ; 1 for cyl_wrt_O ; 2 for cyl_wrt_R
		gmtl::identity(T_cyl);
		T_cyl.setState(gmtl::Matrix44f::AFFINE);
		if (axis == 1) {					// axis: 0 for y-axis ; 1 for x-axis ; 2 for z-axis 
			set_z_rot(pi / 2);
			T_cyl = z_rot * T_cyl;
		}
		else if (axis == 2) {
			set_x_rot(pi / 2);
			T_cyl = x_rot * T_cyl;
		}
		M_cyl_to_world = T_cyl;
	}
	if (frame == 1) {
		gmtl::identity(T_cyl);
		T_cyl[0][0] = T_cyl[1][1] = T_cyl[2][2] = 0.33f;
		T_cyl.setState(gmtl::Matrix44f::AFFINE);
		if (axis == 1) {
			set_z_rot(pi / 2);
			T_cyl = z_rot * T_cyl;
		}
		else if (axis == 2) {
			set_x_rot(pi / 2);
			T_cyl = x_rot * T_cyl;
		}
		M_cyl_to_world = Q[0] * T_O_wrt_R * T_cyl;
	}
	if (frame == 2) {
		gmtl::identity(T_cyl);
		T_cyl[0][0] = T_cyl[2][2] = 1.5f;
		T_cyl[1][1] = 0.5f;
		T_cyl.setState(gmtl::Matrix44f::AFFINE);
		if (axis == 1) {
			set_z_rot(pi / 2);
			T_cyl = z_rot * T_cyl;
		}
		else if (axis == 2) {
			set_x_rot(pi / 2);
			T_cyl = x_rot * T_cyl;
		}
		M_cyl_to_world = Q[0] * T_cyl;
	}
}

void setup_projection() {
	float near, far, left, right, top, bottom;
	near = 1.0;
	far = 5.0;
	left = -1.5;
	right = 1.5;
	top = 1.0;
	bottom = -1.0;

	gmtl::zero(P);
	P[0][0] = (2 * near) / (right - left);
	P[0][2] = (right + left) / (right - left);
	P[1][1] = (2 * near) / (top - bottom);
	P[1][2] = (top + bottom) / (top - bottom);
	P[2][2] = -(far + near) / (far - near);
	P[2][3] = -(2 * far*near) / (far - near);
	P[3][2] = -1.0;
	P.setState(gmtl::Matrix44f::AFFINE);
}

void setup_transform_matrices() {
	gmtl::identity(z_translation);
	z_translation[2][3] = 2.0;
	z_translation.setState(gmtl::Matrix44f::AFFINE);

	// {J1} wrt {bird base}
	gmtl::identity(T_J1_wrt_B);
	T_J1_wrt_B[0][3] = -0.025;
	T_J1_wrt_B[1][3] = -0.1;
	T_J1_wrt_B[2][3] = 0.025;
	T_J1_wrt_B.setState(gmtl::Matrix44f::AFFINE);

	// {J2} wrt {bird base}
	gmtl::identity(T_J2_wrt_B);
	T_J2_wrt_B[0][3] = 0.025;
	T_J2_wrt_B[1][3] = -0.1;
	T_J2_wrt_B[2][3] = 0.025;
	T_J2_wrt_B.setState(gmtl::Matrix44f::AFFINE);

	// {bird base} wrt {O}
	gmtl::identity(T_B_wrt_O);
	set_y_rot(pi / 2);
	T_B_wrt_O = y_rot * T_B_wrt_O;
	T_B_wrt_O.setState(gmtl::Matrix44f::AFFINE);

	// {O} wrt {R}
	gmtl::identity(T_O_wrt_R);
	T_O_wrt_R[0][3] = flight_radius;
	T_O_wrt_R.setState(gmtl::Matrix44f::AFFINE);

	// Q
	for (int i = 0; i < num_birds; i++) {
		set_y_rot(i*(2 * pi) / num_birds);
		set_x_rot(i*(pi) / 4);
		gmtl::identity(Q[i]);
		Q[i] = Q[i] * y_rot * x_rot;
		Q[i].setState(gmtl::Matrix44f::AFFINE);
	}

	// {camera 1} wrt {W}
	gmtl::identity(cam_1_wrt_W);
	cam_1_wrt_W = z_translation;
	cam_1_wrt_W.setState(gmtl::Matrix44f::AFFINE);
	V_cam_1 = gmtl::invert(cam_1_wrt_W);

	// {camera 2} wrt {W}
	gmtl::identity(Cam_2);
	set_z_rot(-pi / 2);
	Cam_2 = z_rot * z_translation;
	Cam_2.setState(gmtl::Matrix44f::AFFINE);
	Cam_2_wrt_W = Q[0] * T_O_wrt_R * Cam_2;
	V_cam_2 = gmtl::invert(Cam_2_wrt_W);

	// scale bird object to bird base frame
	gmtl::identity(T_bird_obj_wrt_B);
	T_bird_obj_wrt_B[0][0] = T_bird_obj_wrt_B[1][1] = T_bird_obj_wrt_B[2][2] = 0.1f;
	T_bird_obj_wrt_B.setState(gmtl::Matrix44f::AFFINE);

	// scale wing 1 to joint 1 frame
	gmtl::identity(T_wing1_obj_wrt_J1);
	T_wing1_obj_wrt_J1[0][0] = T_wing1_obj_wrt_J1[1][1] = T_wing1_obj_wrt_J1[2][2] = 0.1f;
	T_wing1_obj_wrt_J1.setState(gmtl::Matrix44f::AFFINE);

	// scale wing 2 to joint 2 frame
	gmtl::identity(T_wing2_obj_wrt_J2);
	T_wing2_obj_wrt_J2[0][0] = T_wing2_obj_wrt_J2[1][1] = T_wing2_obj_wrt_J2[2][2] = 0.1f;
	T_wing2_obj_wrt_J2.setState(gmtl::Matrix44f::AFFINE);

	// bird model matrix
	for (int i = 0; i < num_birds; i++) {
		gmtl::identity(M_bird[i]);
		M_bird[i].setState(gmtl::Matrix44f::AFFINE);
		M_bird[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_bird_obj_wrt_B;
	}

	// wing 1 model matrix
	for (int i = 0; i < num_birds; i++) {
		gmtl::identity(M_wing1[i]);
		M_wing1[i].setState(gmtl::Matrix44f::AFFINE);
		M_wing1[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_J1_wrt_B * T_wing1_obj_wrt_J1;
	}


	// wing 2 model matrix
	for (int i = 0; i < num_birds; i++) {
		gmtl::identity(M_wing2[i]);
		M_wing2[i].setState(gmtl::Matrix44f::AFFINE);
		M_wing2[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_J2_wrt_B * T_wing2_obj_wrt_J2;
	}

	// sphere model matrix
	gmtl::identity(M_sph_to_world);
	M_sph_to_world.setState(gmtl::Matrix44f::AFFINE);

	// cylinder model matrix 
	gmtl::identity(M_cyl_to_world);
	M_cyl_to_world.setState(gmtl::Matrix44f::AFFINE);

	for (int i = 0; i < num_birds; i++) {
		omega[i] = (i + 1);
	}
}

void roll_bird(double gamma) {
	set_y_rot(gamma);
	T_B_wrt_O = y_rot * T_B_wrt_O;
}

void flap_wing(gmtl::Matrix44f &T, float &alpha, float &alpha_inc) {
	T[2][2] = T[0][0] = cos(alpha);
	T[0][2] = sin(alpha);
	T[2][0] = -sin(alpha);
	T.setState(gmtl::Matrix44f::AFFINE);
	if ((alpha == pi / 2) || (alpha == -pi / 2)) {
		alpha_inc = -alpha_inc;
	}
	alpha = alpha + alpha_inc;
}

void update_cam(double delta_x_pos, double delta_y_pos) {
	if (cam_1_active == true) {
		lon = lon + delta_x_pos;
		lat = lat + delta_y_pos;
		set_y_rot(lon);
		set_x_rot(lat);
		cam_1_wrt_W = y_rot * x_rot *z_translation;
		V_cam_1 = gmtl::invert(cam_1_wrt_W);
	}
	else {
		ele = ele + delta_x_pos;
		azi = azi + delta_y_pos;
		set_x_rot(azi);
		set_y_rot(ele);
		set_z_rot(-pi / 2);
		Cam_2 = x_rot * y_rot * z_rot *z_translation;
		Cam_2_wrt_W = Q[0] * T_O_wrt_R * Cam_2;
		V_cam_2 = gmtl::invert(Cam_2_wrt_W);
	}
}

/* most new code is in this function*/
GLuint update_birds(GLuint force_VAO, GLuint position_loc, GLuint color_loc) {
	float angle[num_birds][num_birds];
	float distance[num_birds][num_birds];
	float influence[2][num_birds][num_birds];
	gmtl::Vec3f X_axis[num_birds];
	gmtl::Vec3f Y_axis[num_birds];
	gmtl::Vec3f Z_axis;
	gmtl::Vec3f U_axis[num_birds][num_birds];
	gmtl::Vec3f dispersion[num_birds][num_birds];
	gmtl::Vec3f centering[num_birds][num_birds];
	gmtl::Vec3f velocity_matching[num_birds][num_birds];
	gmtl::Vec3f total_dispersion[num_birds];
	gmtl::Vec3f total_centering[num_birds];
	gmtl::Vec3f total_velocity_matching[num_birds];
	gmtl::Vec3f world_velocity[num_birds];
	gmtl::Vec3f temp_velocity;
	gmtl::Vec3f acceleration[num_birds];
	gmtl::Matrix33f vel;

	for (int i = 0; i < num_birds; i++) {
		X_axis[i].set(Q[i][0][0], Q[i][1][0], Q[i][2][0]);
		Y_axis[i].set(Q[i][0][1], Q[i][1][1], Q[i][2][1]);
		world_velocity[i] = Y_axis[i] * flight_radius * omega[i];
		gmtl::normalize(X_axis[i]);
		gmtl::normalize(Y_axis[i]);
	}

	for (int i = 0; i < num_birds; i++) {
		for (int j = 0; j < num_birds; j++) {
			cout << i << j << endl;
			if (i == j) {
				angle[i][j] = distance[i][j] = influence[0][i][j] = influence[1][i][j] = 0;
			}
			else {
				angle[i][j] = acos(gmtl::dot(X_axis[i], X_axis[j]));
				distance[i][j] = angle[i][j] * flight_radius;
				if (distance[i][j] > threshold) {
					influence[0][i][j] = 0;
					influence[1][i][j] = 0;
				}
				else if (distance[i][j] < .01) {
					influence[0][i][j] = 1;
					influence[1][i][j] = 0;
				}
				else {
					influence[0][i][j] = 1 / (1 + 4 * distance[i][j] + 5 * (distance[i][j]) * (distance[i][j]));
					influence[1][i][j] = (4 * distance[i][j] + 5 * (distance[i][j]) * (distance[i][j])) / (1 + 4 * distance[i][j] + 5 * (distance[i][j]) * (distance[i][j]));
					U_axis[i][j] = gmtl::makeCross(X_axis[i], X_axis[j]);
					gmtl::normalize(U_axis[i][j]);

					temp_velocity = world_velocity[j];
					vel = gmtl::make< gmtl::Matrix33f >(gmtl::AxisAngle<float>(-angle[i][j], U_axis[i][j]));
					temp_velocity = vel * temp_velocity;
					temp_velocity = temp_velocity - world_velocity[i];

					dispersion[i][j] = -influence[0][i][j] * gmtl::makeCross(U_axis[i][j], X_axis[i]);
					centering[i][j] = influence[1][i][j] * gmtl::makeCross(U_axis[i][j], X_axis[i]);
					velocity_matching[i][j] = influence[0][i][j] * temp_velocity;
					cout << "disp = " << dispersion[i][j] << endl;
					cout << "cent = " << centering[i][j] << endl;
					cout << "vel = " << velocity_matching[i][j] << endl;
				}
			}
		}
		for (int j = 0; j < num_birds; j++) {

			total_dispersion[i] = total_dispersion[i] + dispersion[i][j];
			total_centering[i] = total_centering[i] + centering[i][j];
			total_velocity_matching[i] = total_velocity_matching[i] + velocity_matching[i][j];
		}
		const float t = 1 / (num_birds - 1);
		total_dispersion[i] = 0.2f * total_dispersion[i];
		total_centering[i] = 0.2f * total_centering[i];
		total_velocity_matching[i] = 0.2f * total_velocity_matching[i];
		acceleration[i] = total_dispersion[i] + total_centering[i] + total_velocity_matching[i];
		cout << "tot disp = " << total_dispersion[i] << endl;
		cout << "tot cent = " << total_centering[i] << endl;
		cout << "tot vel = " << total_velocity_matching[i] << endl;
		cout << "accel = " << acceleration[i] << endl;
	}

	for (int i = 0; i < num_birds; i++) {
		world_velocity[i] = world_velocity[i] + (time * acceleration[i]);
	}

	for (int i = 0; i < num_birds; i++) {
		for (int j = 0; j < 3; j++) {
			Q[i][j][2] = gmtl::makeNormal(gmtl::makeCross(X_axis[i], world_velocity[i]))[j];
		}
		Z_axis.set(Q[i][0][2], Q[i][1][2], Q[i][2][2]);
		for (int j = 0; j < 3; j++) {
			Q[i][j][1] = gmtl::makeCross(Z_axis, X_axis[i])[j];
			Y_axis[i].set(Q[i][0][1], Q[i][1][1], Q[i][2][1]);
		}
		omega[i] = gmtl::dot(Y_axis[i], world_velocity[i]) / flight_radius;
	}

	for (int i = 0; i < num_birds; i++) {
		set_z_rot(omega[i] * time);
		Q[i] = Q[i] * z_rot;
		M_bird[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_bird_obj_wrt_B;
		M_wing1[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_J1_wrt_B * T_wing1_obj_wrt_J1;
		M_wing2[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_J2_wrt_B * T_wing2_obj_wrt_J2;
	}

	Cam_2_wrt_W = Q[0] * T_O_wrt_R * Cam_2;
	V_cam_2 = gmtl::invert(Cam_2_wrt_W);

	/*------------update force vao------------------*/
	GLfloat bird_forces[6 * num_birds][3];
	GLfloat force_colors[6 * num_birds][3];
	GLfloat force_index_list[6 * num_birds];
	for (int i = 0; i < num_birds; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				bird_forces[(6 * i) + (2 * j)][k] = T_O_wrt_R[k][3];
				if (j == 0)
					bird_forces[(6 * i) + (2 * j) + 1][k] = total_dispersion[i][k];
				else if (j == 1)
					bird_forces[(6 * i) + (2 * j) + 1][k] = total_centering[i][k];
				else
					bird_forces[(6 * i) + (2 * j) + 1][k] = total_velocity_matching[i][k];
			}
		}
	}

	for (int i = 0; i < num_birds; i++) {
		for (int j = 0; j < 3; j++) {
			if (j == 0) {
				force_colors[(6 * i) + (2 * j)][0] = force_colors[(6 * i) + (2 * j) + 1][0] = 1.0;
				force_colors[(6 * i) + (2 * j)][1] = force_colors[(6 * i) + (2 * j) + 1][1] = 0.0;
				force_colors[(6 * i) + (2 * j)][2] = force_colors[(6 * i) + (2 * j) + 1][2] = 0.0;
			}
			else if (j == 1) {
				force_colors[(6 * i) + (2 * j)][0] = force_colors[(6 * i) + (2 * j) + 1][0] = 0.0;
				force_colors[(6 * i) + (2 * j)][1] = force_colors[(6 * i) + (2 * j) + 1][1] = 1.0;
				force_colors[(6 * i) + (2 * j)][2] = force_colors[(6 * i) + (2 * j) + 1][2] = 0.0;
			}
			else {
				force_colors[(6 * i) + (2 * j)][0] = force_colors[(6 * i) + (2 * j) + 1][0] = 0.0;
				force_colors[(6 * i) + (2 * j)][1] = force_colors[(6 * i) + (2 * j) + 1][1] = 0.0;
				force_colors[(6 * i) + (2 * j)][2] = force_colors[(6 * i) + (2 * j) + 1][2] = 1.0;
			}
		}
	}

	for (int i = 0; i < 9 * num_birds; i = i + 3) {
		force_index_list[i] = i;
		force_index_list[i+1] = i + 1;
		force_index_list[i + 2] = 0xFFFF;
	}

	GLuint VBO[3];;
	glGenBuffers(3, VBO);
	glBindVertexArray(force_VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, 6 * num_birds * 3 * sizeof(GLfloat), bird_forces, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, 6 * num_birds * 3 * sizeof(GLfloat), force_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 9 * num_birds * sizeof(GLuint), force_index_list, GL_STATIC_DRAW);

	//glDeleteBuffers(2, VBO);

	return force_VAO;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	double gamma = 0.0;
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	else if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
		set_z_rot(pi / 6);
		Q[0] = Q[0] * z_rot;
	}
	else if (key == GLFW_KEY_Y && action == GLFW_PRESS) {
		set_x_rot(pi / 6);
		Q[0] = Q[0] * x_rot;
	}
	else if (key == GLFW_KEY_U && action == GLFW_PRESS) {
		set_x_rot(-pi / 6);
		Q[0] = Q[0] * x_rot;
	}
	else if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		gamma = pi / 6;
		roll_bird(gamma);
	}
	else if (key == GLFW_KEY_E && action == GLFW_PRESS) {
		gamma = -pi / 6;
		roll_bird(gamma);
	}
	else if (key == GLFW_KEY_F && action == GLFW_PRESS) {
		flap_wing(T_J1_wrt_B, wing1_angle, wing1_angle_inc);
	}
	else if (key == GLFW_KEY_G && action == GLFW_PRESS) {
		flap_wing(T_J2_wrt_B, wing2_angle, wing2_angle_inc);
	}
	else if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
		cam_1_active = true;
	}
	else if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
		cam_1_active = false;
	}
	else if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_REPEAT) {
		draw_axes = true;
	}
	else if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE) {
		draw_axes = false;
	}
	else if (key == GLFW_KEY_O && action == GLFW_RELEASE) {
		light1 = !light1;
	}
	else if (key == GLFW_KEY_L && action == GLFW_RELEASE) {
		light2 = !light2;
	}
	else if (key == GLFW_KEY_P && action == GLFW_PRESS) {
		update = !update;
	}
	for (int i = 0; i < num_birds; i++) {
		M_bird[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_bird_obj_wrt_B;
		M_wing1[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_J1_wrt_B * T_wing1_obj_wrt_J1;
		M_wing2[i] = Q[i] * T_O_wrt_R * T_B_wrt_O * T_J2_wrt_B * T_wing2_obj_wrt_J2;
	}
	Cam_2_wrt_W = Q[0] * T_O_wrt_R * Cam_2;
	V_cam_2 = gmtl::invert(Cam_2_wrt_W);
}

void scroll_callback(GLFWwindow* window, double xoffset, double zoffset) {
	set_z_translation(-zoffset / 10);
	update_cam(0.0, 0.0);
}

void cursor_callback(GLFWwindow* window, double x_pos, double y_pos) {
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
		double delta_long = 0.0;
		double delta_lat = 0.0;
		delta_long = -(x_pos - cam_previous_x_pos) / (500 / (2 * pi));
		delta_lat = -(y_pos - cam_previous_y_pos) / (500 / (2 * pi));
		cam_previous_x_pos = x_pos;
		cam_previous_y_pos = y_pos;
		update_cam(delta_long, delta_lat);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) {
		cam_previous_x_pos = x_pos;
		cam_previous_y_pos = y_pos;
	}
}

GLuint* setupDrawnObjects(GLuint position_loc, GLuint color_loc, GLuint normals_loc, GLuint uv_loc)
{
	float phi = 0.0;	// initialize phi
	float theta = 0.0;	// initialize theta
	int ind = 0;		// initialize vertex list index counter
	float r_par = 0;	// radius of sphere at different parallels

						/***-----------------------------Setup sphere object----------------------------------***/
	GLfloat sph_vertices[num_sphere_vert][3];	// initialize array of sphere vertices
	GLfloat sph_normals[num_sphere_vert][4];	// initialize array of sphere vertices
	GLfloat uv_coords[num_sphere_vert][2];
	GLfloat sph_colors[num_sphere_vert][3];
	GLuint sph_index_list[sphere_size];
	float magnitude = 0.0;						// matnitude of normal vector
												/*----------Sphere bottom hemisphere-------------*/
	for (int i = 0; i < np / 2; i++) {
		r_par = sin(theta);
		for (int j = 0; j < nm + 1; j++) {
			sph_vertices[ind][0] = (r_par * cos(phi)) * r_sph;
			sph_vertices[ind][1] = -cos(theta) * r_sph;
			sph_vertices[ind][2] = (r_par * -sin(phi)) * r_sph;
			// determine magnitude and normals for sphere
			magnitude = sqrt((sph_vertices[ind][0])*(sph_vertices[ind][0]) + (sph_vertices[ind][1])*(sph_vertices[ind][1]) + (sph_vertices[ind][2])*(sph_vertices[ind][2]));
			sph_normals[ind][0] = sph_vertices[ind][0] / magnitude;
			sph_normals[ind][1] = sph_vertices[ind][1] / magnitude;
			sph_normals[ind][2] = sph_vertices[ind][2] / magnitude;
			sph_normals[ind][3] = 0.0;
			// texture coordinates
			uv_coords[ind][0] = phi / (2 * pi);
			uv_coords[ind][1] = theta / pi;
			// increment index and phi
			ind = ind + 1;
			phi = phi + phi_inc;
		}
		phi = 0.0;
		theta = theta + theta_inc;
	}
	theta = theta - theta_inc;
	for (int j = 0; j < nm + 1; j++) {
		r_par = sin(theta);
		sph_vertices[ind][0] = (r_par * cos(phi)) * r_sph;
		sph_vertices[ind][1] = -cos(theta) * r_sph;
		sph_vertices[ind][2] = (r_par * -sin(phi)) * r_sph;
		sph_normals[ind][0] = 0.0;
		sph_normals[ind][1] = 1.0;
		sph_normals[ind][2] = 0.0;
		sph_normals[ind][3] = 0.0;
		// texture coordinates
		uv_coords[ind][0] = phi / (2 * pi);
		uv_coords[ind][1] = theta / pi;
		// increment index and phi
		ind = ind + 1;
		phi = phi + phi_inc;
	}
	phi = 0.0;
	for (int i = 0; i < 2; i++) {
		r_par = sin(theta);
		for (int j = 0; j < nm + 1; j++) {
			sph_vertices[ind][0] = (r_par * cos(phi)) * 0.75* r_sph;
			sph_vertices[ind][1] = -cos(theta) * r_sph;
			sph_vertices[ind][2] = (r_par * -sin(phi)) * 0.75 * r_sph;
			magnitude = sqrt((sph_vertices[ind][0])*(sph_vertices[ind][0]) + (sph_vertices[ind][2])*(sph_vertices[ind][2]));
			sph_normals[ind][0] = sph_vertices[ind][0] / magnitude;
			sph_normals[ind][1] = 0.0;
			sph_normals[ind][2] = sph_vertices[ind][2] / magnitude;
			sph_normals[ind][3] = 0.0;
			// texture coordinates
			uv_coords[ind][0] = phi / (2 * pi);
			uv_coords[ind][1] = (theta + theta_inc) / pi;
			// increment index and phi
			ind = ind + 1;
			phi = phi + phi_inc;
		}
		phi = 0.0;
		theta = theta + theta_inc;
	}
	theta = theta - theta_inc;
	for (int j = 0; j < nm + 1; j++) {
		r_par = sin(theta);
		sph_vertices[ind][0] = (r_par * cos(phi)) * r_sph;
		sph_vertices[ind][1] = -cos(theta) * r_sph;
		sph_vertices[ind][2] = (r_par * -sin(phi)) * r_sph;
		sph_normals[ind][0] = 0.0;
		sph_normals[ind][1] = -1.0;
		sph_normals[ind][2] = 0.0;
		sph_normals[ind][3] = 0.0;
		// texture coordinates
		uv_coords[ind][0] = phi / (2 * pi);
		uv_coords[ind][1] = (theta + theta_inc) / pi;
		// increment index and phi
		ind = ind + 1;
		phi = phi + phi_inc;
	}
	phi = 0.0;
	for (int i = 0; i < np / 2; i++) {
		r_par = sin(theta);
		for (int j = 0; j < nm + 1; j++) {
			sph_vertices[ind][0] = (r_par * cos(phi)) * r_sph;
			sph_vertices[ind][1] = -cos(theta) * r_sph;
			sph_vertices[ind][2] = (r_par * -sin(phi)) * r_sph;
			// determine magnitude and normals for sphere
			magnitude = sqrt((sph_vertices[ind][0])*(sph_vertices[ind][0]) + (sph_vertices[ind][1])*(sph_vertices[ind][1]) + (sph_vertices[ind][2])*(sph_vertices[ind][2]));
			sph_normals[ind][0] = sph_vertices[ind][0] / magnitude;
			sph_normals[ind][1] = sph_vertices[ind][1] / magnitude;
			sph_normals[ind][2] = sph_vertices[ind][2] / magnitude;
			sph_normals[ind][3] = 0.0;
			// texture coordinates
			uv_coords[ind][0] = phi / (2 * pi);
			uv_coords[ind][1] = theta / pi;
			// increment index and phi
			ind = ind + 1;
			phi = phi + phi_inc;
		}
		phi = 0.0;
		theta = theta + theta_inc;
	}

	/*------------sphere colors-------------*/
	float r = 1.0;
	float g = 0.0;
	float b = 0.0;
	float delta = 1.0 / num_sphere_vert;
	for (int i = 0; i < num_sphere_vert; i++) {
		sph_colors[i][0] = r;
		sph_colors[i][1] = g;
		sph_colors[i][2] = b;
		r = r - delta;
		g = g + delta;
	}
	for (int i = 55; i < 66; i++) {
		sph_colors[i][0] = 0.0;
		sph_colors[i][1] = 1.0;
		sph_colors[i][2] = 0.0;
	}
	for (int i = 66; i < 88; i++) {
		sph_colors[i][0] = 0.0;
		sph_colors[i][1] = 0.0;
		sph_colors[i][2] = 0.5;
	}
	for (int i = 88; i < 99; i++) {
		sph_colors[i][0] = 1.0;
		sph_colors[i][1] = 0.0;
		sph_colors[i][2] = 0.0;
	}

	/*----------sphere index list-------------*/
	sph_index_list[0] = nm + 1;
	sph_index_list[1] = 0;
	for (int i = 2; i < sphere_size; i++) {
		if (i % (2 * (nm + 1) + 1) == (2 * (nm + 1))) {
			sph_index_list[i] = 0xFFFF;
			sph_index_list[i + 1] = sph_index_list[i - 2] + 1;
			sph_index_list[i + 2] = sph_index_list[i - 1] + 1;
			i = i + 2;
		}
		else {
			sph_index_list[i] = sph_index_list[i - 2] + 1;
		}
	}

	/*----------cylinder verts and normals-----------------*/
	phi = 0.0;		// initialize phi
	theta = 0.0;	// initialize theta
	ind = 0;		// initialize vertex list index counter
	GLfloat cyl_vertices[num_cyl_vert][3];
	GLfloat cyl_normals[num_cyl_vert][4];
	for (int i = 0; i < np; i++) {
		for (int j = 0; j < nm + 1; j++) {
			cyl_vertices[ind][0] = r_cyl * cos(phi);	// x coordinate
			cyl_vertices[ind][1] = -cos(theta);			// y coordinate
			cyl_vertices[ind][2] = r_cyl * -sin(phi);	// x coordinate
			magnitude = sqrt((cyl_vertices[ind][0])*(cyl_vertices[ind][0]) + (cyl_vertices[ind][2])*(cyl_vertices[ind][2]));
			cyl_normals[ind][0] = cyl_vertices[ind][0] / magnitude;
			cyl_normals[ind][1] = 0.0;
			cyl_normals[ind][2] = cyl_vertices[ind][2] / magnitude;
			cyl_normals[ind][3] = 0.0;
			// increment index and phi
			ind = ind + 1;
			phi = phi + phi_inc;
		}
		phi = 0.0;
		theta = theta + theta_inc;
	}

	/*--------cylinder index list-------------------------*/
	GLuint cyl_index_list[cyl_size];
	cyl_index_list[0] = nm + 1;
	cyl_index_list[0] = 0;
	for (int i = 2; i < cyl_size; i++) {
		if (i % (2 * (nm + 1) + 1) == (2 * (nm + 1))) {
			cyl_index_list[i] = 0xFFFF;
			cyl_index_list[i + 1] = cyl_index_list[i - 2] + 1;
			cyl_index_list[i + 2] = cyl_index_list[i - 1] + 1;
			i = i + 2;
		}
		else {
			cyl_index_list[i] = cyl_index_list[i - 2] + 1;
		}
	}

	/* ------------Cylinder Colors----------------- */
	GLfloat cyl_white[num_cyl_vert][3];
	GLfloat cyl_blue[num_cyl_vert][3];
	for (int i = 0; i < num_cyl_vert; i++) {
		cyl_white[i][0] = 1.0;
		cyl_white[i][1] = 1.0;
		cyl_white[i][2] = 1.0;
		cyl_blue[i][0] = 0.0;
		cyl_blue[i][1] = 0.0;
		cyl_blue[i][2] = 1.0;
	}



	/* --------------Setup Bird Object------------------------ */
	GLfloat bird_vertices[][3] = {
		{ 0, 1.0, 0 },
		{ -0.25, 0.0, 0.25 },
		{ 0, 0, 0 },
		{ 0.25, 0, 0.25 },
		{ -0.25, -1.0, 0.25 },
		{ 0.0, -1.0, 0.0 },
		{ 0.25, -1.0, 0.25 }
	};
	GLuint bird_index[13] = {
		0, 1, 0, 2, 0, 3, 0xFFFF,
		1, 4, 2, 5, 3, 6
	};
	GLfloat bird_colors[][3] = {
		{ 1.0, 1.0, 1.0 },
		{ 0.0, 0.0, 1.0 },
		{ 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, 1.0 },
		{ 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0 }
	};

	GLfloat bird_normals[][4] = {
		{ 0.0, 1.0, 0.0, 0.0 },
		{ sqrt(2) / 2, 1.0, sqrt(2) / 2, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ -sqrt(2) / 2, 1.0, sqrt(2) / 2, 0.0 },
		{ sqrt(2) / 2, 1.0, sqrt(2) / 2, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ -sqrt(2) / 2, 1.0, sqrt(2) / 2, 0.0 }
	};

	GLfloat bird_uv_coords[][2] = {
		{ 0.5, 1.0 },
		{ 0.25, 0.5 },
		{ 0.5, 0.5 },
		{ 0.75, 0.5 },
		{ 0.25, 0.0 },
		{ 0.5, 0.0 },
		{ 0.75, 0.0 }
	};

	/* -----------------Setup wing 1---------------------------------- */
	GLfloat wing1_vertices[][3] = {
		{ 0.0, 1.0, 0.0 },
		{ -0.25, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 }
	};
	GLuint wing1_index_list[3] = {
		0, 1, 2
	};
	GLfloat wing1_colors[][3] = {
		{ 1.0, 0.0, 1.0 },
		{ 1.0, 0.0, 1.0 },
		{ 1.0, 0.0, 1.0 }
	};

	GLfloat wing1_normals[][4] = {
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 }
	};

	GLfloat wing1_uv_coords[][2] = {
		{ 0.25, 0.5 },
		{ 0.0, 0.0 },
		{ 0.25, 0.0 }
	};

	/* -----------------Setup wing 2---------------------------------- */
	GLfloat wing2_vertices[][3] = {
		{ 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, 0.0 },
		{ 0.25, 0.0, 0.0 }
	};
	GLuint wing2_index_list[3] = {
		0, 1, 2
	};
	GLfloat wing2_colors[][3] = {
		{ 1.0, 1.0, 0.0 },
		{ 1.0, 1.0, 0.0 },
		{ 1.0, 1.0, 0.0 }
	};

	GLfloat wing2_normals[][4] = {
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 }
	};

	GLfloat wing2_uv_coords[][2] = {
		{ 0.75, 0.5 },
		{ 0.75, 0.0 },
		{ 1.0, 0.0 }
	};

	/* -------------------------------Setup VAO's----------------------------------------------- */
	GLuint VAO_out[6], VBO[22], EAB[6];

	glGenVertexArrays(6, VAO_out);
	glGenBuffers(22, VBO);
	glGenBuffers(6, EAB);

	/* -----------------------sphere-------------------------------- */
	glBindVertexArray(VAO_out[0]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, num_sphere_vert * 3 * sizeof(GLfloat), sph_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, num_sphere_vert * 3 * sizeof(GLfloat), sph_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ARRAY_BUFFER, num_sphere_vert * 4 * sizeof(GLfloat), sph_normals, GL_STATIC_DRAW);
	glVertexAttribPointer(normals_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(normals_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[3]);
	glBufferData(GL_ARRAY_BUFFER, num_sphere_vert * 2 * sizeof(GLfloat), uv_coords, GL_STATIC_DRAW);
	glEnableVertexAttribArray(uv_loc);
	glVertexAttribPointer(uv_loc, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sphere_size * sizeof(GLuint), sph_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------white cyliner--------------------------- */
	glBindVertexArray(VAO_out[1]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[4]);
	glBufferData(GL_ARRAY_BUFFER, num_cyl_vert * 3 * sizeof(GLfloat), cyl_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[5]);
	glBufferData(GL_ARRAY_BUFFER, num_cyl_vert * 3 * sizeof(GLfloat), cyl_white, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[6]);
	glBufferData(GL_ARRAY_BUFFER, num_cyl_vert * 4 * sizeof(GLfloat), cyl_normals, GL_STATIC_DRAW);
	glVertexAttribPointer(normals_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(normals_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[1]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, cyl_size * sizeof(GLuint), cyl_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------blue cyliner---------------------------- */
	glBindVertexArray(VAO_out[2]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[7]);
	glBufferData(GL_ARRAY_BUFFER, num_cyl_vert * 3 * sizeof(GLfloat), cyl_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[8]);
	glBufferData(GL_ARRAY_BUFFER, num_cyl_vert * 3 * sizeof(GLfloat), cyl_blue, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[9]);
	glBufferData(GL_ARRAY_BUFFER, num_cyl_vert * 4 * sizeof(GLfloat), cyl_normals, GL_STATIC_DRAW);
	glVertexAttribPointer(normals_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(normals_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[2]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, cyl_size * sizeof(GLuint), cyl_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------bird----------------------------------- */
	glBindVertexArray(VAO_out[3]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[10]);
	glBufferData(GL_ARRAY_BUFFER, 7 * 3 * sizeof(GLfloat), bird_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[11]);
	glBufferData(GL_ARRAY_BUFFER, 7 * 3 * sizeof(GLfloat), bird_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[12]);
	glBufferData(GL_ARRAY_BUFFER, 7 * 4 * sizeof(GLfloat), bird_normals, GL_STATIC_DRAW);
	glVertexAttribPointer(normals_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(normals_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[13]);
	glBufferData(GL_ARRAY_BUFFER, 7 * 2 * sizeof(GLfloat), bird_uv_coords, GL_STATIC_DRAW);
	glEnableVertexAttribArray(uv_loc);
	glVertexAttribPointer(uv_loc, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[3]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 13 * sizeof(GLuint), bird_index, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------wing 1------------------------------ */
	glBindVertexArray(VAO_out[4]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[14]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing1_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[15]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing1_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[16]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 4 * sizeof(GLfloat), wing1_normals, GL_STATIC_DRAW);
	glVertexAttribPointer(normals_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(normals_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[17]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 2 * sizeof(GLfloat), wing1_uv_coords, GL_STATIC_DRAW);
	glEnableVertexAttribArray(uv_loc);
	glVertexAttribPointer(uv_loc, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[4]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * sizeof(GLuint), wing1_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* --------------------wing 2----------------------------- */
	glBindVertexArray(VAO_out[5]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[18]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing2_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[19]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing2_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[20]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 4 * sizeof(GLfloat), wing2_normals, GL_STATIC_DRAW);
	glVertexAttribPointer(normals_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(normals_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[21]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 2 * sizeof(GLfloat), wing2_uv_coords, GL_STATIC_DRAW);
	glEnableVertexAttribArray(uv_loc);
	glVertexAttribPointer(uv_loc, 2, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[5]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * sizeof(GLuint), wing2_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	//glDeleteBuffers(22, VBO);
	//glDeleteBuffers(6, EAB);

	return VAO_out;
}

GLFWwindow* setupWindow()
{
	GLFWwindow *window;
	glfwSetErrorCallback(error_callback);

	if (!glfwInit())
		exit(EXIT_FAILURE);

	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); //Update later
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	window = glfwCreateWindow(1600, 900, "Assignment 5", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwSetKeyCallback(window, key_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwMakeContextCurrent(window);
	//gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
	glfwSwapInterval(1);

	if (glewInit() != GLEW_OK)
		exit(EXIT_FAILURE);

	glEnable(GL_DEPTH_TEST);

	return window;
}

void display(GLFWwindow* window, GLuint shader_program, GLuint m_location, GLuint mvp_location, GLuint norm_mat_location,
	GLuint light_location, GLuint light_color_1_location, GLuint light_color_2_location, GLuint world_origin_location, GLuint texture_location,
	unsigned int tex1_width, unsigned int tex1_height, unsigned char *imagedata1, unsigned int tex2_width, unsigned int tex2_height, unsigned char *imagedata2,
	unsigned int tex3_width, unsigned int tex3_height, unsigned char *imagedata3, GLuint *VAO, GLuint bird_VAO[][3], gmtl::Matrix44f mat[][3])
{
	gmtl::Matrix33f normMat[num_birds + 1][3];
	if (cam_1_active == true)
		V_current = V_cam_1;
	else
		V_current = V_cam_2;
	for (int i = 0; i < num_birds + 1; i++) {
		for (int j = 0; j < 3; j++) {
			mat[i][j] = V_current * mat[i][j];
			normMat[i][j][0][0] = mat[i][j][0][0];
			normMat[i][j][0][1] = mat[i][j][0][1];
			normMat[i][j][0][2] = mat[i][j][0][2];
			normMat[i][j][1][0] = mat[i][j][1][0];
			normMat[i][j][1][1] = mat[i][j][1][1];
			normMat[i][j][1][2] = mat[i][j][1][2];
			normMat[i][j][2][0] = mat[i][j][2][0];
			normMat[i][j][2][1] = mat[i][j][2][1];
			normMat[i][j][2][2] = mat[i][j][2][2];
			normMat[i][j] = gmtl::invert(normMat[i][j]);
		}
	}

	light_position.set(0.0, 0.75, 1.0, 1.0);
	light_position = V_current * light_position;
	world_origin.set(V_current[0][3], V_current[1][3], V_current[2][3], V_current[3][3]);
	light_color_1.set(0.0, 0.0, 0.0);
	light_color_2.set(0.0, 0.0, 0.0);
	if (light1 == true)
		light_color_1.set(1.0, 1.0, 1.0);
	if (light2 == true)
		light_color_2.set(1.0, 1.0, 0.0);



	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(shader_program);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(0xFFFF);

	glUniform4f(light_location, light_position[0], light_position[1], light_position[2], light_position[3]);
	glUniform4f(world_origin_location, world_origin[0], world_origin[1], world_origin[2], world_origin[3]);
	glUniform3f(light_color_1_location, light_color_1[0], light_color_1[1], light_color_1[2]);
	glUniform3f(light_color_2_location, light_color_2[0], light_color_2[1], light_color_2[2]);

	/* --------------------Draw the cylinders-------------------------------- */
	if (draw_axes == true) {
		glBindVertexArray(VAO[1]);

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 3; j++) {
				set_T_cyl(i, j);
				glUniformMatrix4fv(m_location, 1, GL_FALSE, (V_current * M_cyl_to_world).mData);
				glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (P*V_current * M_cyl_to_world).mData);
				glUniformMatrix4fv(norm_mat_location, 1, GL_FALSE, (V_current * M_cyl_to_world).mData);
				glDrawElements(GL_TRIANGLE_STRIP, cyl_size, GL_UNSIGNED_INT, (void*)0);
			}
		}

		glBindVertexArray(VAO[2]);
		for (int i = 0; i < 3; i++) {
			set_T_cyl(2, i);
			glUniformMatrix4fv(m_location, 1, GL_FALSE, (V_current * M_cyl_to_world).mData);
			glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (P*V_current * M_cyl_to_world).mData);
			glUniformMatrix4fv(norm_mat_location, 1, GL_FALSE, (V_current * M_cyl_to_world).mData);
			glDrawElements(GL_TRIANGLE_STRIP, cyl_size, GL_UNSIGNED_INT, (void*)0);
		}
	}

	glUniform1i(texture_location, 0);

	glActiveTexture(GL_TEXTURE0);

	glBindTexture(GL_TEXTURE_2D, texture_location);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex1_height, tex1_width, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	bool bird_text = false; // change to false to render without texture
	if (bird_text == true) {
		glUniform1i(texture_location, 0);

		glActiveTexture(GL_TEXTURE0);

		glBindTexture(GL_TEXTURE_2D, texture_location);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex2_height, tex2_width, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata2);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	}

	for (int i = 0; i < num_birds; i++) {
		/* -------------------------Draw the bird----------------------------- */
		glBindVertexArray(bird_VAO[i][0]);
		glUniformMatrix4fv(m_location, 1, GL_FALSE, (mat[i][0]).mData);
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (P*mat[i][0]).mData);
		glUniformMatrix3fv(norm_mat_location, 1, GL_FALSE, normMat[i][0].mData);
		glDrawElements(GL_TRIANGLE_STRIP, 13, GL_UNSIGNED_INT, (void*)0);

		/* --------------------------Draw wing 1------------------------------ */
		glBindVertexArray(bird_VAO[i][1]);
		glUniformMatrix4fv(m_location, 1, GL_FALSE, (mat[i][1]).mData);
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (P*mat[i][1]).mData);
		glUniformMatrix3fv(norm_mat_location, 1, GL_FALSE, normMat[i][1].mData);
		glDrawElements(GL_TRIANGLE_STRIP, 3, GL_UNSIGNED_INT, (void*)0);

		/* ---------------------------Draw wing 2---------------------------- */
		glBindVertexArray(bird_VAO[i][2]);
		glUniformMatrix4fv(m_location, 1, GL_FALSE, (mat[i][2]).mData);
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (P*mat[i][2]).mData);
		glUniformMatrix3fv(norm_mat_location, 1, GL_FALSE, normMat[i][2].mData);
		glDrawElements(GL_TRIANGLE_STRIP, 3, GL_UNSIGNED_INT, (void*)0);
	}

	/* -------------------------Draw the sphere--------------------------- */
	bool sphere_text = false;  // change to false to render without texture
	if (sphere_text == true) {
		glUniform1i(texture_location, 0);

		glActiveTexture(GL_TEXTURE0);

		glBindTexture(GL_TEXTURE_2D, texture_location);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex3_height, tex3_width, 0, GL_RGB, GL_UNSIGNED_BYTE, imagedata3);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	}

	glBindVertexArray(VAO[0]);
	glUniformMatrix4fv(m_location, 1, GL_FALSE, mat[num_birds][0].mData);
	glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (P*mat[num_birds][0]).mData);
	glUniformMatrix3fv(norm_mat_location, 1, GL_FALSE, normMat[num_birds][0].mData);
	glDrawElements(GL_TRIANGLE_STRIP, sphere_size, GL_UNSIGNED_INT, (void*)0);

	glBindVertexArray(VAO[3]);
	glDrawElements(GL_LINES, 6 * num_birds, GL_UNSIGNED_INT, (void*)0);
	//glDrawArrays(GL_LINES, 0, 6 * num_birds);

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void load_images(unsigned int &texwidth, unsigned int &texheight, unsigned char *&imagedata, char* texture) {
	LoadPPM(texture, &texwidth, &texheight, &imagedata, 1);
}

int main(int argc, char *argv[])
{
	GLFWwindow* mainwindow = NULL;
	GLuint program = NULL, VAO[4], bird_VAO[num_birds][3];
	GLuint pos_location = NULL, color_location = NULL, normals_location = NULL, uv_location = NULL;
	GLuint m_location = NULL, mvp_location = NULL, norm_mat_location = NULL;
	GLuint light_location = NULL, light_color_1_location = NULL, light_color_2_location = NULL, world_origin_location = NULL, texture_location = NULL;
	gmtl::Matrix44f matrix_list[num_birds + 1][3];
	unsigned int tex1_width, tex1_height, tex2_width, tex2_height, tex3_width, tex3_height;
	unsigned char *imagedata1 = NULL, *imagedata2 = NULL, *imagedata3 = NULL;

	/* This function defines the model and transform matrices for the sphere and cylinder */
	setup_transform_matrices();
	setup_projection();

	mainwindow = setupWindow();

	program = setupShaderProgram();

	load_images(tex1_width, tex1_height, imagedata1, "white.ppm");
	load_images(tex2_width, tex2_height, imagedata2, "stars.ppm");
	load_images(tex3_width, tex3_height, imagedata3, "marbles2.ppm");

	init(program, pos_location, color_location, normals_location, uv_location, m_location, mvp_location, norm_mat_location, light_location, light_color_1_location, light_color_2_location, world_origin_location, texture_location);

	for (int i = 0; i < 3; i++) {
		VAO[i] = *(setupDrawnObjects(pos_location, color_location, normals_location, uv_location) + i);
	}
	for (int i = 0; i < num_birds; i++) {
		for (int j = 0; j < 3; j++) {
			bird_VAO[i][j] = *(setupDrawnObjects(pos_location, color_location, normals_location, uv_location) + 3 + j);
		}
	}
	//VAO[3] = *(setupDrawnObjects(pos_location, color_location, normals_location, uv_location) + 6);
	VAO[3] = update_birds(VAO[3], pos_location, color_location);


	while (!glfwWindowShouldClose(mainwindow)) {
		for (int i = 0; i < num_birds; i++) {
			matrix_list[i][0] = M_bird[i];
			matrix_list[i][1] = M_wing1[i];
			matrix_list[i][2] = M_wing2[i];
		}
		matrix_list[num_birds][0] = M_sph_to_world;

		display(mainwindow, program, m_location, mvp_location, norm_mat_location,
			light_location, light_color_1_location, light_color_2_location, world_origin_location, texture_location, tex1_width, tex1_height, imagedata1, tex2_width, tex2_height, imagedata2, tex3_width, tex3_height, imagedata3, VAO, bird_VAO, matrix_list);
		
		glfwSwapBuffers(mainwindow);
		glfwPollEvents();
		if (update) {
			VAO[3] = update_birds(VAO[3], pos_location, color_location);
			//update_birds();
		}
	}
	glfwDestroyWindow(mainwindow);

	glUseProgram(NULL);

	glDisableVertexAttribArray(pos_location);
	glDisableVertexAttribArray(color_location);

	glDeleteProgram(program);
	glDeleteVertexArrays(4, VAO);

	return 0;
}
