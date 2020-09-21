// Nicholas Pyfrom-Day    nap6091
// CMPS 415

#include <GLFW/glfw3.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#pragma comment (lib, "opengl32.lib")
#pragma comment (lib, "glfw3.lib")

#define WIDTH 400		// width of window (also frame buffer's width)
#define HEIGHT 300		// height of window (also frame buffer's height)
static GLubyte frame_buffer[HEIGHT][WIDTH][3];

int click_num = 1;    // variable to distinguish between first and second mouse click
double xs_pos = 0.0;  // starting x position
double ys_pos = 0.0;  // starting y position
double xe_pos = 0.0;  // ending x position
double ye_pos = 0.0;  // ending y position

void error_callback(int error, const char* description)
{
	fprintf(stderr, "Error: %s\n", description);
}

void draw_line(double xi_pos, double yi_pos, double xf_pos, double yf_pos) {
	//set x and y positions
	double x_pos = xi_pos;
	double y_pos = yi_pos;

	// initialize slope, decrement value and decision variable
	double m = (yf_pos - yi_pos) / (xf_pos - xi_pos);
	double w = 1.0 - m;
	double d = 0.5;

	// initialize RGB for endpoints
	double Rs = 255, Re = 0, R = Rs; // Rs is starting value for red, Re is ending value, R is assigned to a pixel and incremented
	double Gs = 0, Ge = 255, G = Gs;
	double Bs = 0, Be = 0, B = Bs;

	// calculate RGB increments
	double dR = (Re - Rs) / (xf_pos - xi_pos); // change in R
	double dG = (Ge - Gs) / (xf_pos - xi_pos); // change in G
	double dB = (Be - Bs) / (xf_pos - xi_pos); // change in B

	//plot initial pixel
	frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = Rs;

	// for vertical line
	if (xf_pos == xi_pos) {
		// calculate RGB increments
		dR = (Re - Rs) / (yf_pos - yi_pos);
		dG = (Ge - Gs) / (yf_pos - yi_pos);
		dB = (Be - Bs) / (yf_pos - yi_pos);
		// vertical line down
		if (yf_pos > yi_pos) {
			for (int i = 1; i < (yf_pos - yi_pos); i++) {
				// increment RGB
				R = R + dR;
				G = G + dG;
				B = B + dB;
				// set color for current pixel
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos][2] = B;
				// antialiasing left of line
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos - 1][0] = 0.5*R + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) - i][(int)x_pos - 1][0]);
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos - 1][1] = 0.5*G + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) - i][(int)x_pos - 1][1]);
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos - 1][2] = 0.5*B + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) - i][(int)x_pos - 1][2]);
				// antialiasing left of line
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos + 1][0] = 0.5*R + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) - i][(int)x_pos + 1][0]);
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos + 1][1] = 0.5*G + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) - i][(int)x_pos + 1][1]);
				frame_buffer[HEIGHT - (int)y_pos - i][(int)x_pos + 1][2] = 0.5*B + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) - i][(int)x_pos + 1][2]);
			}
		}
		else {
			for (int i = 1; i < (yi_pos - yf_pos); i++) {
				// increment RGB
				R = R - dR;
				G = G - dG;
				B = B - dB;
				// set color for current pixel
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos][2] = B;
				// antialiasing left of line
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos - 1][0] = 0.5*R + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) + i][(int)x_pos - 1][0]);
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos - 1][1] = 0.5*G + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) + i][(int)x_pos - 1][1]);
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos - 1][2] = 0.5*B + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) + i][(int)x_pos - 1][2]);
				// antialiasing left of line
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos + 1][0] = 0.5*R + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) + i][(int)x_pos + 1][0]);
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos + 1][1] = 0.5*G + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) + i][(int)x_pos + 1][1]);
				frame_buffer[HEIGHT - (int)y_pos + i][(int)x_pos + 1][2] = 0.5*B + 0.5*(frame_buffer[(HEIGHT - (int)y_pos) + i][(int)x_pos + 1][2]);
			}
		}
	}

	// for horizontal line
	else if (yf_pos == yi_pos) {
		if (xi_pos < xf_pos) {
			for (int i = 1; i < (xf_pos - xi_pos); i++) {
				// increment RGB
				R = R + dR;
				G = G + dG;
				B = B + dB;
				// color current pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + i][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + i][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + i][2] = B;
				//antialiasing above
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos + i][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos + i][0];
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos + i][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos + i][1];
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos + i][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos + i][2];
				//antialiasing below
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos + i][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos + i][0];
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos + i][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos + i][1];
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos + i][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos + i][2];
			}
		}
		else {
			for (int i = 1; i < (xi_pos - xf_pos); i++) {
				// increment RGB
				R = R - dR;
				G = G - dG;
				B = B - dB;
				// color current pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - i][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - i][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - i][2] = B;
				//antialiasing above
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos - i][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos - i][0];
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos - i][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos - i][1];
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos - i][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos - i][2];
				//antialiasing below
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos - i][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos - i][0];
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos - i][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos - i][1];
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos - i][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos - i][2];
			}
		}
	}

	// for diagonal line
	else if (fabs(yf_pos - yi_pos) == fabs(xf_pos - xi_pos)) {
		if (xf_pos > xi_pos) {
			if (yf_pos > yi_pos) {
				while (y_pos < yf_pos) {
					// increment x and y position
					x_pos++;
					y_pos++;
					// increment color channels
					R = R + dR;
					G = G + dG;
					B = B + dB;
					// color current pixel
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;
					// antialiasing above line
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2];
					//antialiasing below
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2];
				}
			}
			else {
				while (y_pos > yf_pos) {
					// increment x and y position
					x_pos++;
					y_pos--;
					// increment color channels
					R = R + dR;
					G = G + dG;
					B = B + dB;
					// color current pixel
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;
					// antialiasing above line
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2];
					//antialiasing below
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2];
				}
			}
		}
		else {
			if (yf_pos > yi_pos) {
				while (y_pos < yf_pos) {
					// increment x and y position
					x_pos--;
					y_pos++;
					// increment color channels
					R = R - dR;
					G = G - dG;
					B = B - dB;
					// color current pixel
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;
					// antialiasing above line
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2];
					//antialiasing below
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2];
				}
			}
			else {
				while (y_pos > yf_pos) {
					// increment x and y position
					x_pos--;
					y_pos--;
					// increment color channels
					R = R - dR;
					G = G - dG;
					B = B - dB;
					// color current pixel
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
					frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;
					// antialiasing above line
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2];
					//antialiasing below
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = 0.5*R + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = 0.5*G + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1];
					frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = 0.5*B + 0.5*frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2];
				}
			}
		}
	}

	// for 0 < slope < 1
	else if (m > 0 && m < 1) {
		if (xi_pos > xf_pos) {
			while (x_pos > xf_pos) {
				x_pos--;
				R = R - dR;
				G = G - dG;
				B = B - dB;
				if (d < w) {
					d = d + m;
				}
				else {
					y_pos--;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0]);
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1]);
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0]);
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1]);
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2]);
			}
		}
		else {
			while (x_pos < xf_pos) {
				x_pos++;
				R = R + dR;
				G = G + dG;
				B = B + dB;
				if (d < w) {
					d = d + m;
				}
				else {
					y_pos++;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0]);
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1]);
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0]);
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1]);
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2]);
			}
		}
	}

	// for -1 < slope < 0
	else if (m < 0 && m > -1) {
		m = -m;
		w = 1.0 - m;
		if (xi_pos > xf_pos) {
			while (x_pos > xf_pos) {
				x_pos--;
				R = R - dR;
				G = G - dG;
				B = B - dB;
				if (d < w) {
					d = d + m;
				}
				else {
					y_pos++;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0]);
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1]);
				frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0]);
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1]);
				frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2]);
			}
		}
	else {
		while (x_pos < xf_pos) {
			x_pos++;
			R = R + dR;
			G = G + dG;
			B = B + dB;
			if (d < w) {
				d = d + m;
			}
			else {
				y_pos--;
				d = d - w;
			}
			//color the pixel
			frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
			frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
			frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

			//anitaliasing below pixel
			frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][0]);
			frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][1]);
			frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)(y_pos + 1)][(int)x_pos][2]);

			//antialiasing above pixel
			frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][0]);
			frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][1]);
			frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)(y_pos - 1)][(int)x_pos][2]);
		}
	}
	}

	// for slope > 1
	else if (m > 1) {
		m = (xf_pos - xi_pos) / (yf_pos - yi_pos);
		w = 1.0 - m;
		// calculate RGB increments
		dR = (Re - Rs) / (yf_pos - yi_pos);
		dG = (Ge - Gs) / (yf_pos - yi_pos);
		dB = (Be - Bs) / (yf_pos - yi_pos);
		if (yi_pos > yf_pos) {
			while (y_pos > yf_pos) {
				y_pos--;
				R = R - dR;
				G = G - dG;
				B = B - dB;
				if (d < w) {
					d = d + m;
				}
				else {
					x_pos--;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2]);
			}
		}
		else {
			while (y_pos < yf_pos) {
				y_pos++;
				R = R + dR;
				G = G + dG;
				B = B + dB;
				if (d < w) {
					d = d + m;
				}
				else {
					x_pos++;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2]);
			}
		}
	}

	// for slope < -1
	else if (m < -1) {
		m = (xi_pos - xf_pos) / (yf_pos - yi_pos);
		w = 1.0 - m;
		// calculate RGB increments
		dR = (Re - Rs) / (yf_pos - yi_pos);
		dG = (Ge - Gs) / (yf_pos - yi_pos);
		dB = (Be - Bs) / (yf_pos - yi_pos);
		if (yi_pos > yf_pos) {
			while (y_pos > yf_pos) {
				y_pos--;
				R = R - dR;
				G = G - dG;
				B = B - dB;
				if (d < w) {
					d = d + m;
				}
				else {
					x_pos++;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2]);
			}
		}
		else {
			while (y_pos < yf_pos) {
				y_pos++;
				R = R + dR;
				G = G + dG;
				B = B + dB;
				if (d < w) {
					d = d + m;
				}
				else {
					x_pos--;
					d = d - w;
				}
				//color the pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][0] = R;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][1] = G;
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos][2] = B;

				//anitaliasing below pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0] = d*R + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1] = d*G + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2] = d*B + (1 - d)*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos - 1][2]);

				//antialiasing above pixel
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0] = (1 - d)*R + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][0]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1] = (1 - d)*G + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][1]);
				frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2] = (1 - d)*B + d*(frame_buffer[HEIGHT - (int)y_pos][(int)x_pos + 1][2]);
			}
		}
	}
	
	//plot final pixel
	frame_buffer[HEIGHT - (int)yf_pos][(int)xf_pos][1] = Ge;
}

/* Called when a mouse button is pressed: */
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		// store poition of first click
		if (click_num == 1) {
			glfwGetCursorPos(window, &xs_pos, &ys_pos);
			click_num = 2;
		}
		// store position of second click and call function to plot the line
		else if (click_num == 2) {
			glfwGetCursorPos(window, &xe_pos, &ye_pos);
			click_num = 1;
			draw_line(xs_pos, ys_pos, xe_pos, ye_pos);
		}
	}
}

/* Called when a keyboard button is pressed: */
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void display(void)
{
	glViewport(0, 0, WIDTH, HEIGHT);
	glClear(GL_COLOR_BUFFER_BIT);

	// Set the raster position to the lower-left corner to avoid a problem
	// (with glDrawPixels) when the window is resized to smaller dimensions.
	glRasterPos2i(-1, -1);

	// Write the information stored in "frame_buffer" to the color buffer
	glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, frame_buffer);
}

int main(void)
{
	GLFWwindow* window;
	glfwSetErrorCallback(error_callback);
	if (!glfwInit())
		exit(EXIT_FAILURE);

	window = glfwCreateWindow(WIDTH, HEIGHT, "Assignment 1", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	glfwMakeContextCurrent(window);

	while (!glfwWindowShouldClose(window))
	{
		display();

		glfwSwapBuffers(window);
		//* For continuous rendering
		glfwPollEvents(); 
		//* For rendering only in response to events
		//glfwWaitEvents();
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}