// Nicholas Pyfrom-Day
// nap6091

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
const float r_sph = 0.4;	// radius of sphere
const float r_cyl = 0.005;	// radius of cylinder
const int nm = 100;			// number of meridians is 10
const int np = 100;			// number of parallels is 10
const int size = 2 * ((np - 1) * (nm + 1)) + np - 2; // size of index list
const float phi_inc = (2 * pi) / nm; // initialize the increment value for phi
const float theta_inc = pi / (np - 1); // initialize the increment value for theta
const int num_vert = np*(nm + 1); // number of vertices

gmtl::Matrix44f M_sph_to_world;	// sphere model matrix
gmtl::Matrix44f M_cyl_to_world;	// cylinder model matrix
gmtl::Matrix44f x_trans_pos;	// transform matrix for positive X translation
gmtl::Matrix44f x_trans_neg;	// transform matrix for negative X translation
gmtl::Matrix44f x_rot_pos;		// transform matrix for positive X rotation
gmtl::Matrix44f x_rot_neg;		// transform matrix for negative X rotation
gmtl::Matrix44f z_rot_pos;		// transform matrix for positive Z rotation
gmtl::Matrix44f z_rot_neg;		// transform matrix for negative Z rotation
gmtl::Matrix44f cyl_x_rot;		// transform matrix for cylinder X rotation
gmtl::Matrix44f cyl_z_rot;		// transform matrix for cylinder Z rotation

char* filetobuf(char *file)
{
	FILE *fptr;
	long length;
	char *buf;

	fptr = fopen(file, "rb"); /* Open file for reading */
	if (!fptr) /* Return NULL on failure */
		return NULL;
	fseek(fptr, 0, SEEK_END); /* Seek to the end of the file */
	length = ftell(fptr); /* Find out how many bytes into the file we are */
	buf = (char*)malloc(length + 1); /* Allocate a buffer for the entire length of the file and a null terminator */
	fseek(fptr, 0, SEEK_SET); /* Go back to the beginning of the file */
	fread(buf, length, 1, fptr); /* Read the contents of the file in to the buffer */
	fclose(fptr); /* Close the file */
	buf[length] = 0; /* Null terminator */

	return buf; /* Return the buffer */
}

void error_callback(int error, const char* description)
{
	fprintf(stderr, "Error: %s\n", description);
}

/* This function defines the model and transform matrices for the sphere and cylinder */
void setup_transform_matrices() {
	// set sphere model matrix to identity
	gmtl::identity(M_sph_to_world);
	M_sph_to_world.setState(gmtl::Matrix44f::AFFINE);

	// set cylinder model matrix to identity
	gmtl::identity(M_cyl_to_world);
	M_cyl_to_world.setState(gmtl::Matrix44f::AFFINE);

	// matrix for positive X translation
	gmtl::identity(x_trans_pos);
	x_trans_pos[0][3] = 0.1;
	x_rot_pos.setState(gmtl::Matrix44f::AFFINE);

	// matrix for negative X translation
	gmtl::identity(x_trans_neg);
	x_trans_neg[0][3] = -0.1;
	x_rot_neg.setState(gmtl::Matrix44f::AFFINE);

	// matrix for positive Z rotation
	gmtl::identity(z_rot_pos);
	z_rot_pos[1][1] = z_rot_pos[0][0] = cos(pi / 6);
	z_rot_pos[0][1] = -sin(pi / 6);
	z_rot_pos[1][0] = sin(pi / 6);
	z_rot_pos.setState(gmtl::Matrix44f::AFFINE);

	// matrix for negative Z rotation
	gmtl::identity(z_rot_neg);
	z_rot_neg[1][1] = z_rot_neg[0][0] = cos(-pi / 6);
	z_rot_neg[0][1] = -sin(-pi / 6);
	z_rot_neg[1][0] = sin(-pi / 6);
	z_rot_neg.setState(gmtl::Matrix44f::AFFINE);

	// matrix for positive X rotation
	gmtl::identity(x_rot_pos);
	x_rot_pos[1][1] = x_rot_pos[2][2] = cos(pi / 6);
	x_rot_pos[1][2] = -sin(pi / 6);
	x_rot_pos[2][1] = sin(pi / 6);
	x_rot_pos.setState(gmtl::Matrix44f::AFFINE);

	// matrix for negative X rotation
	gmtl::identity(x_rot_neg);
	x_rot_neg[1][1] = x_rot_neg[2][2] = cos(-pi / 6);
	x_rot_neg[1][2] = -sin(-pi / 6);
	x_rot_neg[2][1] = sin(-pi / 6);
	x_rot_neg.setState(gmtl::Matrix44f::AFFINE);

	// cylinder Y to Z axis (x-rotation)
	gmtl::identity(cyl_x_rot);
	cyl_x_rot[1][1] = cyl_x_rot[2][2] = cos(pi / 2);
	cyl_x_rot[1][2] = -sin(pi / 2);
	cyl_x_rot[2][1] = sin(pi / 2);
	cyl_x_rot.setState(gmtl::Matrix44f::AFFINE);

	// cylinder Y to X axis (z-rotation)
	gmtl::identity(cyl_z_rot);
	cyl_z_rot[1][1] = cyl_z_rot[2][2] = cos(pi / 2);
	cyl_z_rot[0][1] = -sin(pi / 2);
	cyl_z_rot[1][0] = sin(pi / 2);
	cyl_z_rot.setState(gmtl::Matrix44f::AFFINE);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	// positive world-X translation of sphere
	else if (key == GLFW_KEY_X && action == GLFW_PRESS) {
		M_sph_to_world = x_trans_pos * M_sph_to_world;
	}
	// negative world-X translation of sphere
	else if (key == GLFW_KEY_Z && action == GLFW_PRESS) {
		M_sph_to_world = x_trans_neg * M_sph_to_world;
	}
	// positive local-X translation of sphere
	else if (key == GLFW_KEY_V && action == GLFW_PRESS) {
		M_sph_to_world = M_sph_to_world * x_trans_pos;
	}
	// negative local-X translation of sphere
	else if (key == GLFW_KEY_C && action == GLFW_PRESS) {
		M_sph_to_world = M_sph_to_world * x_trans_neg;
	}
	// positive world-Z rotation of sphere
	else if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
		M_sph_to_world = z_rot_pos * M_sph_to_world;
	}
	// negative world-Z rotation of sphere
	else if (key == GLFW_KEY_E && action == GLFW_PRESS) {
		M_sph_to_world = z_rot_neg * M_sph_to_world;
	}
	// positive local-X rotation of sphere
	else if (key == GLFW_KEY_S && action == GLFW_PRESS) {
		M_sph_to_world = M_sph_to_world * x_rot_pos;
	}
	// negative local-X rotation of sphere
	else if (key == GLFW_KEY_W && action == GLFW_PRESS) {
		M_sph_to_world = M_sph_to_world * x_rot_neg;
	}
	// positive local-Z rotation of sphere
	else if (key == GLFW_KEY_A && action == GLFW_PRESS) {
		M_sph_to_world = M_sph_to_world * z_rot_pos;
	}
	// negative local-Z rotation of sphere
	else if (key == GLFW_KEY_D && action == GLFW_PRESS) {
		M_sph_to_world = M_sph_to_world * z_rot_neg;
	}
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
	glfwMakeContextCurrent(window);
	//gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
	glfwSwapInterval(1);

	if (glewInit() != GLEW_OK)
		exit(EXIT_FAILURE);

	glEnable(GL_DEPTH_TEST);

	return window;
}

GLuint setupShaderProgram()
{
	GLuint vertex_shader, fragment_shader, shader_program;
	int IsCompiled_VS, IsCompiled_FS, IsLinked, max_length;
	char *vertex_shader_log;
	char *fragment_shader_log;
	char *shader_program_log;

	/* Read our shaders into the appropriate buffers */
	char* vertex_source = filetobuf("OpenGL_Example.vert");
	char* fragment_source = filetobuf("OpenGL_Example.frag");

	vertex_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex_shader, 1, &vertex_source, NULL);
	glCompileShader(vertex_shader);

	glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &IsCompiled_VS);
	if (IsCompiled_VS == GL_FALSE)
	{
		glGetShaderiv(vertex_shader, GL_INFO_LOG_LENGTH, &max_length);

		/* The max_length includes the NULL character */
		vertex_shader_log = (char *)malloc(max_length);

		glGetShaderInfoLog(vertex_shader, max_length, &max_length, vertex_shader_log);
		printf("Error: %s", vertex_shader_log);
		/* Handle the error in an appropriate way such as displaying a message or writing to a log file. */
		/* In this simple program, we'll just leave */
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

		/* The max_length includes the NULL character */
		fragment_shader_log = (char *)malloc(max_length);

		glGetShaderInfoLog(fragment_shader, max_length, &max_length, fragment_shader_log);
		printf("Error: %s", fragment_shader_log);
		/* Handle the error in an appropriate way such as displaying a message or writing to a log file. */
		/* In this simple program, we'll just leave */
		free(fragment_shader_log);
		free(vertex_source);
		free(fragment_source);
		return 0;
	}

	/* If we reached this point it means the vertex and fragment shaders compiled and are syntax error free. */
	/* We must link them together to make a GL shader program */
	/* GL shader programs are monolithic. It is a single piece made of 1 vertex shader and 1 fragment shader. */
	/* Assign our program handle a "name" */
	shader_program = glCreateProgram();

	/* Attach our shaders to our program */
	glAttachShader(shader_program, vertex_shader);
	glAttachShader(shader_program, fragment_shader);

	/* Link our program */
	/* At this stage, the vertex and fragment programs are inspected, optimized and a binary code is generated for the shader. */
	/* The binary code is uploaded to the GPU, if there is no error. */
	glLinkProgram(shader_program);

	/* Again, we must check and make sure that it linked. If it fails, it would mean either there is a mismatch between the vertex */
	/* and fragment shaders. It might be that you have surpassed your GPU's abilities. Perhaps too many ALU operations or */
	/* too many texel fetch instructions or too many interpolators or dynamic loops. */

	glGetProgramiv(shader_program, GL_LINK_STATUS, (int *)&IsLinked);
	if (IsLinked == GL_FALSE)
	{
		/* Noticed that glGetProgramiv is used to get the length for a shader program, not glGetShaderiv. */
		glGetProgramiv(shader_program, GL_INFO_LOG_LENGTH, &max_length);

		/* The max_length includes the NULL character */
		shader_program_log = (char *)malloc(max_length);

		/* Notice that glGetProgramInfoLog, not glGetShaderInfoLog. */
		glGetProgramInfoLog(shader_program, max_length, &max_length, shader_program_log);
		printf("Error: %s", shader_program_log);
		/* Handle the error in an appropriate way such as displaying a message or writing to a log file. */
		/* In this simple program, we'll just leave */
		free(shader_program_log);
		free(vertex_source);
		free(fragment_source);
		return 0;
	}
	//	glBindAttribLocation(shader_program, SHADER_POSITION_INDEX, "in_position");
	//	glBindAttribLocation(shader_program, SHADER_COLOR_INDEX, "in_color");

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

GLuint* setupDrawnObjects(GLuint position_loc, GLuint color_loc)
{
	float phi = 0.0; // initialize phi
	float theta = 0.0; // initialize theta
	int ind = 0; // initialize vertex list index counter
	float r_par = 0;  // radius of sphere at different parallels

	GLfloat sph_vertices[num_vert][3]; // initialize array of sphere vertices
	GLfloat cyl_vertices[num_vert][3]; // initialize array of cylinder vertices
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
		// reset phi and increment theta
		phi = 0.0;
		theta = theta + theta_inc;
	}

	GLuint sph_index_list[size]; // index array for spher
	GLuint cyl_index_list[size]; // index array for cylinder
	cyl_index_list[0] = sph_index_list[0] = nm + 1;
	cyl_index_list[0] = sph_index_list[1] = 0;
	for (int i = 2; i < size; i++) {
		if (i % (2 * (nm + 1) + 1) == (2 * (nm + 1))) {
			cyl_index_list[i] = sph_index_list[i] = 0xFFFF;
			cyl_index_list[i + 1] = sph_index_list[i + 1] = sph_index_list[i - 2] + 1;
			cyl_index_list[i + 2] = sph_index_list[i + 2] = sph_index_list[i - 1] + 1;
			i = i + 2;
		}
		else {
			cyl_index_list[i] = sph_index_list[i] = sph_index_list[i - 2] + 1;
		}
	}

	GLfloat sph_colors[num_vert][3]; // color array for sphere
	GLfloat cyl_colors[num_vert][3]; // color array for cylinder
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
		cyl_colors[i][0] = 1.0;
		cyl_colors[i][1] = 1.0;
		cyl_colors[i][2] = 1.0;
	}

	GLuint VAO_out[2], VBO[4], EAB[2]; 
	glGenVertexArrays(2, VAO_out);
	glGenBuffers(4, VBO);
	glGenBuffers(2, EAB);

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

	glBindVertexArray(VAO_out[1]);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), cyl_vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(position_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(position_loc);

	glBindBuffer(GL_ARRAY_BUFFER, VBO[3]);
	glBufferData(GL_ARRAY_BUFFER, num_vert * 3 * sizeof(GLfloat), cyl_colors, GL_STATIC_DRAW);
	glVertexAttribPointer(color_loc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(color_loc);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EAB[1]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(GLuint), cyl_index_list, GL_STATIC_DRAW);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	/*
	These calls just remove the names/IDs from use. The buffered data are still associated
	with the vertex array object. Since they are only scoped to this function, however,
	we would normally remove them here or GL will never use them again within the application.
	This probably wouldn't cause errors for the assignment, so we will omit it here.
	*/
	// glDeleteBuffers(2, VBO);
	// glDeleteBuffers(1, &EAB);

	return VAO_out;
}

void display(GLFWwindow* window, GLuint shader_program, GLuint m_location, GLuint VAO_sph, GLuint VAO_cyl, gmtl::Matrix44f mat_sph, gmtl::Matrix44f mat_cyl, gmtl::Matrix44f mat_cyl_x, gmtl::Matrix44f mat_cyl_z)
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(shader_program);

	// Enable primitive restart
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(0xFFFF);

	// bind cylinder vao
	glBindVertexArray(VAO_cyl);

	/* Draw cylinder for Y-axis */
	glUniformMatrix4fv(
		m_location,		
		1,				
		GL_FALSE,		
		mat_cyl.mData		
	);
	glDrawElements(
		GL_TRIANGLE_STRIP,		
		size,					
		GL_UNSIGNED_INT,		
		(void*)0				
	);

	/* Draw cylinder for X-axis */
	glUniformMatrix4fv(
		m_location,		
		1,				
		GL_FALSE,		
		mat_cyl_z.mData		
	);
	glDrawElements(
		GL_TRIANGLE_STRIP,		
		size,					
		GL_UNSIGNED_INT,		
		(void*)0				
	);

	/* Draw cylinder for Z-axis */
	glUniformMatrix4fv(
		m_location,		
		1,			
		GL_FALSE,		
		mat_cyl_x.mData		
	);
	glDrawElements(
		GL_TRIANGLE_STRIP,		
		size,					
		GL_UNSIGNED_INT,		
		(void*)0			
	);

	/* Draw the sphere */
	glBindVertexArray(VAO_sph);
	glUniformMatrix4fv(
		m_location,		
		1,				
		GL_FALSE,		
		mat_sph.mData		
	);
	glDrawElements(
		GL_TRIANGLE_STRIP,	
		size,					
		GL_UNSIGNED_INT,		
		(void*)0				
	);

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

int main(int argc, char *argv[])
{
	GLFWwindow* mainwindow = NULL;
	GLuint program = NULL, VAO[2] = { NULL, NULL };
	GLuint pos_location = NULL, color_location = NULL, m_location = NULL;

	setup_transform_matrices();

	mainwindow = setupWindow();

	program = setupShaderProgram();

	init(program, pos_location, color_location, m_location);

	VAO[0] = *(setupDrawnObjects(pos_location, color_location));		// sphere VAO
	VAO[1] = *(setupDrawnObjects(pos_location, color_location) + 1);	// cylinder VAO


	while (!glfwWindowShouldClose(mainwindow))
	{	display(mainwindow, program, m_location, VAO[0], VAO[1], M_sph_to_world, M_cyl_to_world, cyl_x_rot, cyl_z_rot);

		glfwSwapBuffers(mainwindow);
		glfwPollEvents();
	}
	glfwDestroyWindow(mainwindow);

	glUseProgram(NULL);

	glDisableVertexAttribArray(pos_location);
	glDisableVertexAttribArray(color_location);

	glDeleteProgram(program);
	glDeleteVertexArrays(2, VAO);

	return 0;
}