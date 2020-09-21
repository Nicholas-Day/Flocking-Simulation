// Nicholas Pyfrom-Day		nap6091
// CMPS 415
// Assignment 3 source code

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <gmtl/gmtl.h> 
#include <stdlib.h>
#include <stdio.h>

#pragma comment (lib, "opengl32.lib")
#pragma comment (lib, "glew32.lib")
#pragma comment (lib, "glfw3.lib")

const float pi = M_PI;
const float r_sph = 0.4;								// radius of sphere
const float r_cyl = 0.005;								// radius of cylinder
const int nm = 100;										// number of meridians is 100
const int np = 100;										// number of parallels is 100
const int size = 2 * ((np - 1) * (nm + 1)) + np - 2;	// size of index list for sphere and cylinder
const int num_vert = np*(nm + 1);						// number of vertices
const float phi_inc = (2 * pi) / nm;					// initialize the increment value for phi
const float theta_inc = pi / (np - 1);					// initialize the increment value for theta
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

/* -------------------------Model Matrices-------------------------------------- */
gmtl::Matrix44f M_sph_to_world;		// sphere model matrix
gmtl::Matrix44f M_cyl_to_world;		// cylinder model matrix
gmtl::Matrix44f M_bird;				// bird model matrix
gmtl::Matrix44f M_wing1;			// wing 1 model matrix
gmtl::Matrix44f M_wing2;			// wing 2 model matrix
/* --------------------------Transforms----------------------------------------- */
gmtl::Matrix44f T_cyl;					// cyl object wrt W, R, or O frame
gmtl::Matrix44f T_bird_obj_wrt_B;		// bird object wrt bird base frame
gmtl::Matrix44f T_wing1_obj_wrt_J1;		// wing object 1 wrt joint 1 frame
gmtl::Matrix44f T_wing2_obj_wrt_J2;		// wing object 2 wrt joint 2 frame
gmtl::Matrix44f T_J1_wrt_B;				// {J1} wrt {bird base}
gmtl::Matrix44f T_J2_wrt_B;				// {J2} wrt {bird base}
gmtl::Matrix44f T_B_wrt_O;				// {bird base} wrt {O}
gmtl::Matrix44f T_O_wrt_R;				// {O} wrt {R}
gmtl::Matrix44f Q;						// {R} wrt {W}
gmtl::Matrix44f cam_1_wrt_W;			// camera pose describing {camera 1} wrt {W}
gmtl::Matrix44f Cam_2;					// camera pose describing {camera 2} wrt {O}
gmtl::Matrix44f Cam_2_wrt_W;			// camera pose describing {camera 2} wrt {W}
/* --------------------------View Matrices-------------------------------------- */
gmtl::Matrix44f V_cam_1;			// view matrix for camera 1
gmtl::Matrix44f V_cam_2;			// view matrix for camera 2
gmtl::Matrix44f V_current;			// the active view matrix
/* --------------------------Rotation Matrices---------------------------------- */
gmtl::Matrix44f x_rot;		// transform matrix for positive X rotation
gmtl::Matrix44f z_rot;		// transform matrix for positive Z rotation
gmtl::Matrix44f y_rot;

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

void init(GLuint shader_program, GLuint &pos_loc_out, GLuint &color_loc_out, GLuint &m_loc_out)
{
	pos_loc_out = glGetAttribLocation(shader_program, "in_position");
	color_loc_out = glGetAttribLocation(shader_program, "in_color");
	m_loc_out = glGetUniformLocation(shader_program, "M");
}

void setup_transform_matrices() {
	// {J1} wrt {bird base}
	gmtl::identity(T_J1_wrt_B);
	T_J1_wrt_B[0][3] = -0.05;
	T_J1_wrt_B[1][3] = -0.2;
	T_J1_wrt_B[2][3] = 0.05;
	T_J1_wrt_B.setState(gmtl::Matrix44f::AFFINE);

	// {J2} wrt {bird base}
	gmtl::identity(T_J2_wrt_B);
	T_J2_wrt_B[0][3] = 0.05;
	T_J2_wrt_B[1][3] = -0.2;
	T_J2_wrt_B[2][3] = 0.05;
	T_J2_wrt_B.setState(gmtl::Matrix44f::AFFINE);

	// {bird base} wrt {O}
	gmtl::identity(T_B_wrt_O);
	T_B_wrt_O.setState(gmtl::Matrix44f::AFFINE);

	// {O} wrt {R}
	gmtl::identity(T_O_wrt_R);
	T_O_wrt_R.setState(gmtl::Matrix44f::AFFINE);
	T_O_wrt_R[0][3] = r_sph + 0.15;

	// {R} wrt {W}
	gmtl::identity(Q);
	Q.setState(gmtl::Matrix44f::AFFINE);

	// {camera 1} wrt {W}
	gmtl::identity(cam_1_wrt_W);
	cam_1_wrt_W.setState(gmtl::Matrix44f::AFFINE);

	// {camera 2} wrt {W}
	gmtl::identity(Cam_2);
	Cam_2.setState(gmtl::Matrix44f::AFFINE);

	// scale bird object to bird base frame
	gmtl::identity(T_bird_obj_wrt_B);
	T_bird_obj_wrt_B[0][0] = T_bird_obj_wrt_B[1][1] = T_bird_obj_wrt_B[2][2] = 0.2f;
	T_bird_obj_wrt_B.setState(gmtl::Matrix44f::AFFINE);

	// scale wing 1 to joint 1 frame
	gmtl::identity(T_wing1_obj_wrt_J1);
	T_wing1_obj_wrt_J1[0][0] = T_wing1_obj_wrt_J1[1][1] = T_wing1_obj_wrt_J1[2][2] = 0.2f;
	T_wing1_obj_wrt_J1.setState(gmtl::Matrix44f::AFFINE);

	// scale wing 2 to joint 2 frame
	gmtl::identity(T_wing2_obj_wrt_J2);
	T_wing2_obj_wrt_J2[0][0] = T_wing2_obj_wrt_J2[1][1] = T_wing2_obj_wrt_J2[2][2] = 0.2f;
	T_wing2_obj_wrt_J2.setState(gmtl::Matrix44f::AFFINE);

	// bird model matrix
	gmtl::identity(M_bird);
	M_bird.setState(gmtl::Matrix44f::AFFINE);
	M_bird = T_O_wrt_R * T_bird_obj_wrt_B;

	// wing 1 model matrix
	gmtl::identity(M_wing1);
	M_wing1.setState(gmtl::Matrix44f::AFFINE);
	M_wing1 = Q * T_O_wrt_R * T_B_wrt_O * T_J1_wrt_B * T_wing1_obj_wrt_J1;

	// wing 2 model matrix
	gmtl::identity(M_wing2);
	M_wing2.setState(gmtl::Matrix44f::AFFINE);
	M_wing2 = Q * T_O_wrt_R * T_B_wrt_O * T_J2_wrt_B * T_wing2_obj_wrt_J2;

	// sphere model matrix
	gmtl::identity(M_sph_to_world);
	M_sph_to_world.setState(gmtl::Matrix44f::AFFINE);

	// cylinder model matrix 
	gmtl::identity(M_cyl_to_world);
	M_cyl_to_world.setState(gmtl::Matrix44f::AFFINE);
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
		M_cyl_to_world = Q * T_O_wrt_R * T_cyl;
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
		M_cyl_to_world = Q * T_cyl;
	}
}

void roll_bird(double gamma) {
	set_y_rot(gamma);
	T_B_wrt_O = y_rot * T_B_wrt_O;
}

void flap_wing(gmtl::Matrix44f &T,float &alpha, float &alpha_inc) {
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
		cam_1_wrt_W = y_rot * x_rot;
		V_cam_1 = gmtl::invert(cam_1_wrt_W);
	}
	else {
		ele = ele + delta_x_pos;
		azi = azi + delta_y_pos;
		set_x_rot(azi);
		set_y_rot(ele);
		set_z_rot(-pi / 2);
		Cam_2 = x_rot * y_rot * z_rot;
		Cam_2_wrt_W = Q * T_O_wrt_R * Cam_2;
		V_cam_2 = gmtl::invert(Cam_2_wrt_W);
	}
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	double gamma = 0.0;
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	else if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
		set_z_rot(pi / 6);
		Q = Q*z_rot;
	}
	else if (key == GLFW_KEY_Y && action == GLFW_PRESS) {
		set_x_rot(pi / 6);
		Q = Q*x_rot;
	}
	else if (key == GLFW_KEY_U && action == GLFW_PRESS) {
		set_x_rot(-pi / 6);
		Q = Q * x_rot;
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
	M_bird = Q * T_O_wrt_R * T_B_wrt_O * T_bird_obj_wrt_B;
	M_wing1 = Q * T_O_wrt_R * T_B_wrt_O * T_J1_wrt_B * T_wing1_obj_wrt_J1;
	M_wing2 = Q * T_O_wrt_R * T_B_wrt_O * T_J2_wrt_B * T_wing2_obj_wrt_J2;
	Cam_2_wrt_W = Q * T_O_wrt_R * Cam_2;
	V_cam_2 = gmtl::invert(Cam_2_wrt_W);
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

GLuint* setupDrawnObjects(GLuint position_loc, GLuint color_loc) 
{
	float phi = 0.0;	// initialize phi
	float theta = 0.0;	// initialize theta
	int ind = 0;		// initialize vertex list index counter
	float r_par = 0;	// radius of sphere at different parallels

	/* ----------Sphere and Cylinder vertex lists---------------------------------- */
	GLfloat sph_vertices[num_vert][3];	// initialize array of sphere vertices
	GLfloat cyl_vertices[num_vert][3];	// initialize array of cylinder vertices
	for (int i = 0; i < np; i++) {
		r_par = sin(theta);
		for (int j = 0; j < nm + 1; j++) {
				sph_vertices[ind][0] = (r_par * cos(phi)) * r_sph;	// x coordinate
				sph_vertices[ind][1] = -cos(theta) * r_sph;			// y coordinate
				sph_vertices[ind][2] = (r_par * -sin(phi)) * r_sph;	// z coordinate

				cyl_vertices[ind][0] = r_cyl * cos(phi);	// x coordinate
				cyl_vertices[ind][1] = -cos(theta);			// y coordinate
				cyl_vertices[ind][2] = r_cyl * -sin(phi);	// x coordinate

			// increment index and phi
			ind = ind + 1;
			phi = phi + phi_inc;
		}
		phi = 0.0;
		theta = theta + theta_inc;
	}

	/* ---------Sphere and Cylinder index lists-------------------------- */
	GLuint sph_index_list[size];
	GLuint cyl_index_list[size];
	cyl_index_list[0] = sph_index_list[0] = nm + 1;
	cyl_index_list[0] = sph_index_list[1] = 0;
	for (int i = 2; i < size; i++) {
		if (i % (2 * (nm + 1) + 1) == (2 * (nm + 1))) {
			cyl_index_list[i] = sph_index_list[i] = 0xFFFF;
			cyl_index_list[i+1] = sph_index_list[i + 1] = sph_index_list[i - 2] + 1;
			cyl_index_list[i+2] = sph_index_list[i + 2] = sph_index_list[i - 1] + 1;
			i = i + 2;
		}
		else {
			cyl_index_list[i] = sph_index_list[i] = sph_index_list[i - 2] + 1;
		}
	}
	
	/* ------------Sphere and Cylinder Colors------------------------- */
	GLfloat sph_colors[num_vert][3];
	GLfloat cyl_white[num_vert][3];
	GLfloat cyl_blue[num_vert][3];
	float r = 1.0;
	float g = 0.0;
	float b = 0.0;
	float delta = 1.0 / (np*(nm + 1));
	for (int i = 0; i < num_vert; i++) {
		sph_colors[i][0] = r;
		sph_colors[i][1] = g;
		sph_colors[i][2] = b;
		r = r - delta;
		g = g + delta;
		cyl_white[i][0] = 1.0;
		cyl_white[i][1] = 1.0;
		cyl_white[i][2] = 1.0;
		cyl_blue[i][0] = 0.0;
		cyl_blue[i][1] = 0.0;
		cyl_blue[i][2] = 1.0;
	}

	/* ---------Setup Bird Vertices, color, and index list--------------- */
	GLfloat bird_vertices[][3] = {
		{ 0, 1.0, 0 },
		{ -0.25, 0.0, 0.25 },
		{ 0, 0, 0 },
		{ 0.25, 0, 0.25 },
		{ -0.25, -1.0, 0.25},
		{ 0.0, -1.0, 0.0 },
		{ 0.25, -1.0, 0.25}
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

	/* -----------------Setup wing 1---------------------------------- */
	GLfloat wing1_vertices[][3] = {
		{ 0.0, 1.0, 0.0},
		{ -0.25, 0.0, 0.0},
		{ 0.0, 0.0, 0.0}
	};
	GLuint wing1_index_list[3] = {
		0, 1, 2
	};
	GLfloat wing1_colors[][3] = {
		{ 1.0, 0.0, 1.0 },
		{ 1.0, 0.0, 1.0 },
		{ 1.0, 0.0, 1.0 }
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

	/* -------------------------------Setup VAO's----------------------------------------------- */
	GLuint VAO_out[6], VBO[12], EAB[6]; 

	glGenVertexArrays(6, VAO_out);
	glGenBuffers(12, VBO);
	glGenBuffers(6, EAB);

	/* -----------------------sphere-------------------------------- */
	glBindVertexArray(VAO_out[0]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), sph_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), sph_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(GLuint), sph_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------white cyliner--------------------------- */
	glBindVertexArray(VAO_out[1]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), cyl_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[3]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), cyl_white, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[1]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(GLuint), cyl_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------blue cyliner---------------------------- */
	glBindVertexArray(VAO_out[2]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[4]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), cyl_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[5]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), cyl_blue, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[2]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(GLuint), cyl_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------bird----------------------------------- */
	glBindVertexArray(VAO_out[3]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[6]);
	glBufferData(GL_ARRAY_BUFFER, 7*3 * sizeof(GLfloat), bird_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[7]);
	glBufferData(GL_ARRAY_BUFFER, 7*3 * sizeof(GLfloat), bird_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[3]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 13 * sizeof(GLuint), bird_index, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* -------------------wing 1------------------------------ */
	glBindVertexArray(VAO_out[4]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[8]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing1_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[9]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing1_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[4]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * sizeof(GLuint), wing1_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/* --------------------wing 2----------------------------- */
	glBindVertexArray(VAO_out[5]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[10]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing2_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[11]);
	glBufferData(GL_ARRAY_BUFFER, 3 * 3 * sizeof(GLfloat), wing2_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[5]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * sizeof(GLuint), wing2_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// glDeleteBuffers(10, VBO);
	// glDeleteBuffers(5, EAB);

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
	window = glfwCreateWindow(500, 500, "Assignment 2", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwMakeContextCurrent(window);
	//gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
	glfwSwapInterval(1);

	if (glewInit() != GLEW_OK)
		exit(EXIT_FAILURE);

	glEnable(GL_DEPTH_TEST);

	return window;
}

void display(GLFWwindow* window, GLuint shader_program, GLuint m_location, GLuint *VAO, gmtl::Matrix44f *mat)
{
	if (cam_1_active == true)
		V_current = V_cam_1;
	else
		V_current = V_cam_2;
	for (int i = 0; i < 4; i++) {
		mat[i] = V_current * mat[i];
	}
	
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(shader_program);

	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(0xFFFF);

	/* --------------------Draw the cylinders-------------------------------- */
	if (draw_axes == true) {
		glBindVertexArray(VAO[1]);

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 3; j++) {
				set_T_cyl(i, j);
				glUniformMatrix4fv(m_location, 1, GL_FALSE, (V_current * M_cyl_to_world).mData);
				glDrawElements(GL_TRIANGLE_STRIP, size, GL_UNSIGNED_INT, (void*)0);
			}
		}

		glBindVertexArray(VAO[2]);
		for (int i = 0; i < 3; i++) {
			set_T_cyl(2, i);
			glUniformMatrix4fv(m_location, 1, GL_FALSE, (V_current * M_cyl_to_world).mData);
			glDrawElements(GL_TRIANGLE_STRIP, size, GL_UNSIGNED_INT, (void*)0);
		}
	}

	/* -------------------------Draw the sphere--------------------------- */
	glBindVertexArray(VAO[0]);
	glUniformMatrix4fv(m_location, 1, GL_FALSE, mat[0].mData);
	glDrawElements(GL_TRIANGLE_STRIP, size, GL_UNSIGNED_INT, (void*)0);

	/* -------------------------Draw the bird----------------------------- */
	glBindVertexArray(VAO[3]);
	glUniformMatrix4fv(m_location, 1, GL_FALSE, mat[1].mData);
	glDrawElements(GL_TRIANGLE_STRIP, 13, GL_UNSIGNED_INT, (void*)0);

	/* --------------------------Draw wing 1------------------------------ */
	glBindVertexArray(VAO[4]);
	glUniformMatrix4fv(m_location, 1, GL_FALSE, mat[2].mData);
	glDrawElements(GL_TRIANGLE_STRIP,3,GL_UNSIGNED_INT,(void*)0);

	/* ---------------------------Draw wing 2---------------------------- */
	glBindVertexArray(VAO[5]);
	glUniformMatrix4fv(m_location, 1, GL_FALSE, mat[3].mData);
	glDrawElements(GL_TRIANGLE_STRIP, 3, GL_UNSIGNED_INT, (void*)0);

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

int main(int argc, char *argv[])
{
	GLFWwindow* mainwindow = NULL;
	GLuint program = NULL, VAO[6];
	GLuint pos_location = NULL, color_location = NULL, m_location = NULL;
	gmtl::Matrix44f matrix_list[4];

	/* This function defines the model and transform matrices for the sphere and cylinder */
	setup_transform_matrices();

	mainwindow = setupWindow();

	program = setupShaderProgram();

	init(program, pos_location, color_location, m_location);

	for (int i = 0; i < 6; i++) {
		VAO[i] = *(setupDrawnObjects(pos_location, color_location) + i);
	}

	while (!glfwWindowShouldClose(mainwindow)){
		matrix_list[0] = M_sph_to_world;
		matrix_list[1] = M_bird;
		matrix_list[2] = M_wing1;
		matrix_list[3] = M_wing2;

		display(mainwindow, program, m_location, VAO, matrix_list);
		
		glfwSwapBuffers(mainwindow);
		glfwPollEvents();
	}
	glfwDestroyWindow(mainwindow);

	glUseProgram(NULL);

	glDisableVertexAttribArray(pos_location);
	glDisableVertexAttribArray(color_location);

	glDeleteProgram(program);
	glDeleteVertexArrays(5, VAO);

	return 0;
}
